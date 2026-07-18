# CUDA backend for the MHD solver (mhd.jl), FP64 or FP32 (to_gpu's T).
# Not included by the FractalToroid module: run scripts include this file
# only when a gpu flag is passed, so CPU-only runs never load CUDA. The
# step mirrors
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
                      MHDSim, cellvol, center, ssprk3!, _lincomb!, _VIRIDIS

export to_gpu, download!, gpu_step!

struct GPUState{A<:AbstractArray{<:AbstractFloat,3}}
    S::Vector{A}
    S1::Vector{A}                   # RK scratch; frame path reuses 1:3 for
    S2::Vector{A}                   # J, 4:6 for ω, S2[1:4] for magnitudes
    K::Vector{A}
    F::Vector{A}                    # face fluxes for the current axis
    tmp::A                          # scratch: signal speeds, mask^dt
    mask::Union{Nothing,A}
    xr::A                           # cell-center coords, (n,1,1)/(1,n,1)/
    yr::A                           # (1,1,n) shaped for broadcasting
    zr::A
    cmap::AbstractMatrix{Float64}   # viridis table on the device
end

function to_gpu(sim::MHDSim; T::Type{<:AbstractFloat} = Float64)
    CUDA.functional() || error("CUDA is not functional on this machine")
    n = sim.box.n
    z8() = [CUDA.zeros(T, n, n, n) for _ in 1:8]
    c = T[center(sim.box, i) for i in 1:n]
    GPUState([CuArray{T,3}(A) for A in sim.S], z8(), z8(), z8(), z8(),
             CUDA.zeros(T, n, n, n),
             sim.mask === nothing ? nothing : CuArray{T,3}(sim.mask),
             CuArray(reshape(c, n, 1, 1)), CuArray(reshape(c, 1, n, 1)),
             CuArray(reshape(c, 1, 1, n)), CuArray(_VIRIDIS))
end

"Copy the evolving fields back into the CPU sim (diagnostics/checkpoint).
For T=Float32 this widens; checkpoints are always written FP64, so a
resumed gpu32 run restarts from the truncated FP32 state."
function download!(sim::MHDSim, G::GPUState)
    for q in 1:8
        A = G.S[q]
        eltype(A) === Float64 ? copyto!(sim.S[q], A) :
                                copyto!(sim.S[q], Array(A))
    end
    sim
end

# ssprk3! works on the GPU state vectors once _lincomb! knows CuArrays
function FractalToroid._lincomb!(Y::Vector{<:CuArray}, a, A, b, B)
    aT = eltype(Y[1])(a)
    bT = eltype(Y[1])(b)
    for f in eachindex(Y)
        @. Y[f] = aT * A[f] + bT * B[f]
    end
end

# typed-zero minmod (the CPU _minmod returns a Float64 literal zero, which
# would promote FP32 kernels back to FP64)
@inline _mm(a, b) = a * b <= 0 ? zero(a) : (abs(a) < abs(b) ? a : b)

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
@inline function _flux8(U::NTuple{8,T}, ax, cs2::T, ch2::T, fl::T) where {T}
    ρ = max(U[1], fl)
    vd = U[1+ax] / ρ
    Bd = U[4+ax]
    b2 = U[5]^2 + U[6]^2 + U[7]^2
    pt = cs2 * ρ + T(0.5) * b2
    F = (U[1+ax],
         U[2] * vd - U[5] * Bd + (ax == 1 ? pt : zero(T)),
         U[3] * vd - U[6] * Bd + (ax == 2 ? pt : zero(T)),
         U[4] * vd - U[7] * Bd + (ax == 3 ? pt : zero(T)),
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
    h = eltype(S[1])(0.5)
    UL = ntuple(Val(8)) do q
        @inbounds u0 = S[q][c0]; @inbounds u1 = S[q][c1]
        @inbounds u2 = S[q][c2]
        u1 + h * _mm(u1 - u0, u2 - u1)
    end
    UR = ntuple(Val(8)) do q
        @inbounds u1 = S[q][c1]; @inbounds u2 = S[q][c2]
        @inbounds u3 = S[q][c3]
        u2 - h * _mm(u2 - u1, u3 - u2)
    end
    FL, sL = _flux8(UL, ax, cs2, ch2, fl)
    FR, sR = _flux8(UR, ax, cs2, ch2, fl)
    α = max(sL, sR, ch)
    for q in 1:8
        @inbounds F[q][c1] = h * (FL[q] + FR[q]) - h * α * (UR[q] - UL[q])
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
    T = eltype(G.tmp)
    idx = T(1 / sim.box.dx)
    cs2 = T(sim.cs^2)
    ch = T(sim.ch)
    fl = T(sim.rho_floor)
    for q in 1:8
        fill!(K[q], zero(T))
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
        ηi2 = T(sim.eta / sim.box.dx^2)
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
    T = eltype(G.tmp)
    cs2 = T(sim.cs^2)
    fl = T(sim.rho_floor)
    @. G.tmp = max(abs(S[MMX]), abs(S[MMY]), abs(S[MMZ])) / max(S[MRHO], fl) +
               sqrt(cs2 + (S[MBX]^2 + S[MBY]^2 + S[MBZ]^2) / max(S[MRHO], fl))
    Float64(maximum(G.tmp))
end

"GPU mirror of mhd_step!: scalars stay in `sim`, fields stay on the device."
function gpu_step!(sim::MHDSim, G::GPUState)
    smax = max(1e-10, gpu_max_speed(G, sim))
    sim.ch = smax
    dt = sim.cfl * sim.box.dx / smax
    sim.dt = dt
    rhs! = (K, S, t) -> gpu_rhs!(K, S, G, sim)
    ssprk3!(G.S, G.S1, G.S2, G.K, rhs!, dt, sim.t)
    T = eltype(G.tmp)
    damp = T(exp(-0.1 * sim.ch * dt / sim.box.dx))
    G.S[MPSI] .*= damp
    G.S[MRHO] .= max.(G.S[MRHO], T(sim.rho_floor))
    if G.mask !== nothing
        @. G.tmp = G.mask^T(dt)
        removed = 0.0
        for q in (MMX, MMY, MMZ, MBX, MBY, MBZ, MPSI)
            A = G.S[q]
            removed += mapreduce((a, m) -> (a * a * (1 - m * m)) / 2, +,
                                 A, G.tmp)
            A .*= G.tmp
        end
        sim.absorbed += removed * cellvol(sim.box)
    end
    sim.t += dt
    sim
end

# ---- frame-time products on the device -----------------------------------
# Everything the per-frame path needs (magnitudes, curls, moments, energies,
# the 3D volume render) computed on the GPU into the RK scratch arrays,
# with only small results downloaded. CPU-side mirror: the frame branch of
# scripts/v2_pathB_selfassembly.jl.

"Central-difference curl (mirror of curl_central), periodic."
function k_curl!(Cx, Cy, Cz, Fx, Fy, Fz, n, h)
    t = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    t > n * n * n && return
    i, j, k = _decode(t, n)
    ip = _wp(i, n); im = _wm(i, n)
    jp = _wp(j, n); jm = _wm(j, n)
    kp = _wp(k, n); km = _wm(k, n)
    @inbounds begin
        Cx[i, j, k] = (Fz[i, jp, k] - Fz[i, jm, k]) * h -
                      (Fy[i, j, kp] - Fy[i, j, km]) * h
        Cy[i, j, k] = (Fx[i, j, kp] - Fx[i, j, km]) * h -
                      (Fz[ip, j, k] - Fz[im, j, k]) * h
        Cz[i, j, k] = (Fy[ip, j, k] - Fy[im, j, k]) * h -
                      (Fx[i, jp, k] - Fx[i, jm, k]) * h
    end
    return
end

# Orthographic emission–absorption raycaster (mirror of volume_render;
# ray math in FP64 for parity with the CPU renderer regardless of state T)
function k_volren!(img, color, opacity, cmap, n, lo, dx, opmax,
                   res, nstep, κ, op_gamma, clo, span, logcolor, bg,
                   dirx, diry, dirz, ex1, ex2, ex3, ey1, ey2, ey3,
                   cx, diag, ds)
    t = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    t > res * res && return
    px = (t - 1) % res + 1
    py = (t - 1) ÷ res + 1
    u = (2 * (px - 0.5) / res - 1) * diag
    v = (2 * (py - 0.5) / res - 1) * diag
    x = cx + u * ex1 + v * ey1 - diag * dirx
    y = cx + u * ex2 + v * ey2 - diag * diry
    z = cx + u * ex3 + v * ey3 - diag * dirz
    r = 0.0; g = 0.0; b = 0.0; T = 1.0
    for _ in 1:nstep
        x += dirx * ds; y += diry * ds; z += dirz * ds
        fx = (x - lo) / dx - 0.5
        fy = (y - lo) / dx - 0.5
        fz = (z - lo) / dx - 0.5
        (0.0 <= fx <= n - 1.001 && 0.0 <= fy <= n - 1.001 &&
         0.0 <= fz <= n - 1.001) || continue
        i0 = unsafe_trunc(Int32, fx) + Int32(1); tx = fx - (i0 - 1)
        j0 = unsafe_trunc(Int32, fy) + Int32(1); ty = fy - (j0 - 1)
        k0 = unsafe_trunc(Int32, fz) + Int32(1); tz = fz - (k0 - 1)
        op = 0.0; cv = 0.0
        for dk in Int32(0):Int32(1), dj in Int32(0):Int32(1),
            di in Int32(0):Int32(1)
            w = (di == 1 ? tx : 1 - tx) * (dj == 1 ? ty : 1 - ty) *
                (dk == 1 ? tz : 1 - tz)
            @inbounds op += w * Float64(opacity[i0+di, j0+dj, k0+dk])
            @inbounds cv += w * Float64(color[i0+di, j0+dj, k0+dk])
        end
        op <= 0 && continue
        a = 1 - exp(-κ * (op / opmax)^op_gamma * ds)
        a < 1e-4 && continue
        vv = logcolor ? log10(max(cv, 1e-300)) : cv
        vn = clamp((vv - clo) / span, 0.0, 1.0)
        nc = size(cmap, 1)
        vc = vn * (nc - 1)
        ic = min(unsafe_trunc(Int32, vc), Int32(nc - 2))
        f = vc - ic
        @inbounds cr = cmap[ic+1, 1] * (1 - f) + cmap[ic+2, 1] * f
        @inbounds cg = cmap[ic+1, 2] * (1 - f) + cmap[ic+2, 2] * f
        @inbounds cb = cmap[ic+1, 3] * (1 - f) + cmap[ic+2, 3] * f
        w = T * a
        r += w * 255cr; g += w * 255cg; b += w * 255cb
        T *= 1 - a
        T < 0.01 && break
    end
    row = res - py + 1
    @inbounds img[1, px, row] = UInt8(round(clamp(r + T * 255bg, 0.0, 255.0)))
    @inbounds img[2, px, row] = UInt8(round(clamp(g + T * 255bg, 0.0, 255.0)))
    @inbounds img[3, px, row] = UInt8(round(clamp(b + T * 255bg, 0.0, 255.0)))
    return
end

"Device version of volume_render; returns the RGB image on the host."
function gpu_volume_render(color, opacity, G::GPUState, box;
                           azim = 0.6, elev = 0.45, res = 384, nstep = 256,
                           op_scale = 18.0, op_gamma = 0.5, clo = 0.0,
                           chi = 1.0, logcolor = false, bg = 0.07)
    n = box.n
    L = n * box.dx
    opmax = Float64(maximum(opacity))
    opmax <= 0 && (opmax = 1.0)
    κ = op_scale / L
    ca, sa = cos(azim), sin(azim)
    ce, se = cos(elev), sin(elev)
    half = L / 2
    diag = sqrt(3.0) * half
    ds = 2 * diag / nstep
    span = chi - clo > 0 ? chi - clo : 1.0
    img = CuArray{UInt8,3}(undef, 3, res, res)
    @cuda threads = _nthreads() blocks = cld(res * res, _nthreads()) k_volren!(
        img, color, opacity, G.cmap, n, box.lo, box.dx, opmax,
        res, nstep, κ, op_gamma, clo, span, logcolor, bg,
        -ca * ce, -sa * ce, -se, -sa, ca, 0.0, -ca * se, -sa * se, ce,
        box.lo + half, diag, ds)
    Array(img)
end

"""
    frame_render_products!(G, sim; chi, res, azim, elev)

Compute |B|, |ω|=|∇×(ρv)|, |J|²=|∇×B|² and the 3D volume render on the
device (J is left in S1[1:3] for gpu_frame_scalars!). Returns host-side
FP64 copies (Bmag, Wmag, J2, img3d).
"""
function frame_render_products!(G::GPUState, sim::MHDSim;
                                chi, res = 448, azim = 0.6, elev = 0.45)
    T = eltype(G.tmp)
    n = sim.box.n
    h = T(1 / (2 * sim.box.dx))
    S = G.S
    @. G.S2[1] = sqrt(S[MBX]^2 + S[MBY]^2 + S[MBZ]^2)
    @cuda threads = _nthreads() blocks = _nblocks(n) k_curl!(
        G.S1[4], G.S1[5], G.S1[6], S[MMX], S[MMY], S[MMZ], n, h)
    @. G.S2[2] = sqrt(G.S1[4]^2 + G.S1[5]^2 + G.S1[6]^2)
    @cuda threads = _nthreads() blocks = _nblocks(n) k_curl!(
        G.S1[1], G.S1[2], G.S1[3], S[MBX], S[MBY], S[MBZ], n, h)
    @. G.S2[3] = G.S1[1]^2 + G.S1[2]^2 + G.S1[3]^2
    img = gpu_volume_render(G.S2[2], G.S2[1], G, sim.box;
                            res, chi, azim, elev)
    Float64.(Array(G.S2[1])), Float64.(Array(G.S2[2])),
    Float64.(Array(G.S2[3])), img
end

"""
    gpu_frame_scalars!(G, sim) -> (m, T, E_kin, E_mag)

Dipole/anapole moments of J (from S1[1:3], filled by
frame_render_products!) and the energies, reduced on the device.
Mirrors grid_moments / mhd_kinetic_energy / mhd_magnetic_energy.
"""
function gpu_frame_scalars!(G::GPUState, sim::MHDSim)
    S = G.S
    Jx, Jy, Jz = G.S1[1], G.S1[2], G.S1[3]
    x, y, z = G.xr, G.yr, G.zr
    w = G.S2[4]
    dV = cellvol(sim.box)
    @. w = y * Jz - z * Jy; m1 = Float64(sum(w))
    @. w = z * Jx - x * Jz; m2 = Float64(sum(w))
    @. w = x * Jy - y * Jx; m3 = Float64(sum(w))
    @. w = (x * Jx + y * Jy + z * Jz) * x - 2 * (x^2 + y^2 + z^2) * Jx
    T1 = Float64(sum(w))
    @. w = (x * Jx + y * Jy + z * Jz) * y - 2 * (x^2 + y^2 + z^2) * Jy
    T2 = Float64(sum(w))
    @. w = (x * Jx + y * Jy + z * Jz) * z - 2 * (x^2 + y^2 + z^2) * Jz
    T3 = Float64(sum(w))
    fl = eltype(G.tmp)(sim.rho_floor)
    @. w = (S[MMX]^2 + S[MMY]^2 + S[MMZ]^2) / max(S[MRHO], fl)
    ekin = 0.5 * Float64(sum(w)) * dV
    emag = 0.5 * dV * Float64(sum(abs2, S[MBX]) + sum(abs2, S[MBY]) +
                              sum(abs2, S[MBZ]))
    (dV / 2) .* [m1, m2, m3], (dV / 10) .* [T1, T2, T3], ekin, emag
end

end # module
