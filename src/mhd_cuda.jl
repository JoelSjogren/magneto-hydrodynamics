# CUDA backend for the MHD solver (mhd.jl), FP64. Not included by the
# FractalToroid module: run scripts include this file only when the "gpu"
# flag is passed, so CPU-only runs never load CUDA. The step mirrors
# mhd_step! operation-for-operation; diagnostics, rendering, ICs and
# checkpoints stay on the CPU — the evolving fields live on the device and
# are copied back with download! at frame/checkpoint times.
#
# The CPU flux routine accumulates each face's flux into both neighbor
# cells while sweeping a pencil serially — a scatter that races between
# GPU threads. Here axis fluxes go to face arrays first (one thread per
# face), then a second kernel differences them into K.

module MHDCuda

using CUDA
using FractalToroid
import FractalToroid: MRHO, MMX, MMY, MMZ, MBX, MBY, MBZ, MPSI,
                      MHDSim, cellvol, ssprk3!, _lincomb!, _minmod

export to_gpu, download!, gpu_step!

struct GPUState{A<:AbstractArray{Float64,3}}
    S::Vector{A}
    S1::Vector{A}
    S2::Vector{A}
    K::Vector{A}
    F::Vector{A}                    # face fluxes for the current axis
    tmp::A                          # scratch: signal speeds, mask^dt
    mask::Union{Nothing,A}
end

function to_gpu(sim::MHDSim)
    CUDA.functional() || error("CUDA is not functional on this machine")
    n = sim.box.n
    z8() = [CUDA.zeros(Float64, n, n, n) for _ in 1:8]
    GPUState([CuArray(A) for A in sim.S], z8(), z8(), z8(), z8(),
             CUDA.zeros(Float64, n, n, n),
             sim.mask === nothing ? nothing : CuArray(sim.mask))
end

"Copy the evolving fields back into the CPU sim (diagnostics/checkpoint)."
function download!(sim::MHDSim, G::GPUState)
    for q in 1:8
        copyto!(sim.S[q], G.S[q])
    end
    sim
end

# ssprk3! works on the GPU state vectors once _lincomb! knows CuArrays
function FractalToroid._lincomb!(Y::Vector{<:CuArray}, a, A, b, B)
    for f in eachindex(Y)
        @. Y[f] = a * A[f] + b * B[f]
    end
end

@inline _wp(i, n) = i == n ? 1 : i + 1
@inline _wm(i, n) = i == 1 ? n : i - 1

@inline function _neighbors(i, j, k, ax, n)
    if ax == 1
        (CartesianIndex(_wm(i, n), j, k), CartesianIndex(i, j, k),
         CartesianIndex(_wp(i, n), j, k), CartesianIndex(_wp(_wp(i, n), n), j, k))
    elseif ax == 2
        (CartesianIndex(i, _wm(j, n), k), CartesianIndex(i, j, k),
         CartesianIndex(i, _wp(j, n), k), CartesianIndex(i, _wp(_wp(j, n), n), k))
    else
        (CartesianIndex(i, j, _wm(k, n)), CartesianIndex(i, j, k),
         CartesianIndex(i, j, _wp(k, n)), CartesianIndex(i, j, _wp(_wp(k, n), n)))
    end
end

"Tuple version of _flux8!: GLM-MHD flux of state U along axis ax."
@inline function _flux8(U::NTuple{8,Float64}, ax, cs2, ch2, fl)
    ρ = max(U[1], fl)
    vd = U[1+ax] / ρ
    Bd = U[4+ax]
    b2 = U[5]^2 + U[6]^2 + U[7]^2
    pt = cs2 * ρ + 0.5 * b2
    F = (U[1+ax],
         U[2] * vd - U[5] * Bd + (ax == 1 ? pt : 0.0),
         U[3] * vd - U[6] * Bd + (ax == 2 ? pt : 0.0),
         U[4] * vd - U[7] * Bd + (ax == 3 ? pt : 0.0),
         ax == 1 ? U[8] : U[5] * vd - (U[2] / ρ) * Bd,
         ax == 2 ? U[8] : U[6] * vd - (U[3] / ρ) * Bd,
         ax == 3 ? U[8] : U[7] * vd - (U[4] / ρ) * Bd,
         ch2 * Bd)
    F, abs(vd) + sqrt(cs2 + b2 / ρ)
end

@inline function _decode(t, n)
    i = (t - 1) % n + 1
    j = ((t - 1) ÷ n) % n + 1
    k = (t - 1) ÷ (n * n) + 1
    i, j, k
end

# MUSCL + Rusanov flux through the face between cell c1 and its +ax
# neighbor c2, stored at c1's index.
function k_faceflux!(F, S, ax, n, cs2, ch, ch2, fl)
    t = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    t > n * n * n && return
    i, j, k = _decode(t, n)
    c0, c1, c2, c3 = _neighbors(i, j, k, ax, n)
    UL = ntuple(Val(8)) do q
        @inbounds u0 = S[q][c0]; @inbounds u1 = S[q][c1]
        @inbounds u2 = S[q][c2]
        u1 + 0.5 * _minmod(u1 - u0, u2 - u1)
    end
    UR = ntuple(Val(8)) do q
        @inbounds u1 = S[q][c1]; @inbounds u2 = S[q][c2]
        @inbounds u3 = S[q][c3]
        u2 - 0.5 * _minmod(u2 - u1, u3 - u2)
    end
    FL, sL = _flux8(UL, ax, cs2, ch2, fl)
    FR, sR = _flux8(UR, ax, cs2, ch2, fl)
    α = max(sL, sR, ch)
    for q in 1:8
        @inbounds F[q][c1] = 0.5 * (FL[q] + FR[q]) - 0.5 * α * (UR[q] - UL[q])
    end
    return
end

# K[cell] += (flux in through the −ax face − flux out through the +ax face)/dx
function k_facediff!(K, F, ax, n, idx)
    t = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    t > n * n * n && return
    i, j, k = _decode(t, n)
    c = CartesianIndex(i, j, k)
    cm = ax == 1 ? CartesianIndex(_wm(i, n), j, k) :
         ax == 2 ? CartesianIndex(i, _wm(j, n), k) :
                   CartesianIndex(i, j, _wm(k, n))
    for q in 1:8
        @inbounds K[q][c] += (F[q][cm] - F[q][c]) * idx
    end
    return
end

"η∇²B for one component."
function k_laplacian!(dB, B, n, ηi2)
    t = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    t > n * n * n && return
    i, j, k = _decode(t, n)
    @inbounds dB[i, j, k] += ηi2 *
        (B[_wp(i, n), j, k] + B[_wm(i, n), j, k] +
         B[i, _wp(j, n), k] + B[i, _wm(j, n), k] +
         B[i, j, _wp(k, n)] + B[i, j, _wm(k, n)] - 6 * B[i, j, k])
    return
end

@inline _nthreads() = 256
@inline _nblocks(n) = cld(n * n * n, _nthreads())

function gpu_rhs!(K, S, G::GPUState, sim::MHDSim)
    n = sim.box.n
    idx = 1 / sim.box.dx
    cs2 = sim.cs^2
    ch = sim.ch
    fl = sim.rho_floor
    for q in 1:8
        fill!(K[q], 0.0)
    end
    Kt = ntuple(q -> K[q], Val(8))
    St = ntuple(q -> S[q], Val(8))
    Ft = ntuple(q -> G.F[q], Val(8))
    for ax in 1:3
        @cuda threads = _nthreads() blocks = _nblocks(n) k_faceflux!(
            Ft, St, ax, n, cs2, ch, ch * ch, fl)
        @cuda threads = _nthreads() blocks = _nblocks(n) k_facediff!(
            Kt, Ft, ax, n, idx)
    end
    if sim.eta != 0
        ηi2 = sim.eta / sim.box.dx^2
        for q in MBX:MBZ
            @cuda threads = _nthreads() blocks = _nblocks(n) k_laplacian!(
                K[q], S[q], n, ηi2)
        end
    end
    K
end

"Same signal speed max as max_speed(sim, S), reduced on the device."
function gpu_max_speed(G::GPUState, sim::MHDSim)
    S = G.S
    cs2 = sim.cs^2
    fl = sim.rho_floor
    @. G.tmp = max(abs(S[MMX]), abs(S[MMY]), abs(S[MMZ])) / max(S[MRHO], fl) +
               sqrt(cs2 + (S[MBX]^2 + S[MBY]^2 + S[MBZ]^2) / max(S[MRHO], fl))
    maximum(G.tmp)
end

"GPU mirror of mhd_step!: scalars stay in `sim`, fields stay on the device."
function gpu_step!(sim::MHDSim, G::GPUState)
    smax = max(1e-10, gpu_max_speed(G, sim))
    sim.ch = smax
    dt = sim.cfl * sim.box.dx / smax
    sim.dt = dt
    rhs! = (K, S, t) -> gpu_rhs!(K, S, G, sim)
    ssprk3!(G.S, G.S1, G.S2, G.K, rhs!, dt, sim.t)
    damp = exp(-0.1 * sim.ch * dt / sim.box.dx)
    G.S[MPSI] .*= damp
    G.S[MRHO] .= max.(G.S[MRHO], sim.rho_floor)
    if G.mask !== nothing
        @. G.tmp = G.mask^dt
        removed = 0.0
        for q in (MMX, MMY, MMZ, MBX, MBY, MBZ, MPSI)
            A = G.S[q]
            removed += mapreduce((a, m) -> 0.5 * a * a * (1 - m * m), +,
                                 A, G.tmp)
            A .*= G.tmp
        end
        sim.absorbed += removed * cellvol(sim.box)
    end
    sim.t += dt
    sim
end

end # module
