# Maxwell on a periodic Yee-staggered cubic grid, dimensionless (c = ε0 = μ0
# = 1), integrated with SSP-RK3 (method of lines; stable for the staggered
# curl spectrum at CFL ≲ 0.5). Open boundaries are emulated by a "sponge"
# damping layer; the energy it removes is the radiated energy.
#
# Staggering (cell size dx, array index i ↔ coordinate lo + (i-1)dx):
#   Ex(i+½,j,k)  Ey(i,j+½,k)  Ez(i,j,k+½)
#   Bx(i,j+½,k+½)  By(i+½,j,k+½)  Bz(i+½,j+½,k)
# Charge/fluid quantities live on nodes (i,j,k).

# State-vector layout: a state is a Vector{Array{Float64,3}}.
const EX, EY, EZ, BX, BY, BZ = 1, 2, 3, 4, 5, 6
const FN, FPX, FPY, FPZ = 7, 8, 9, 10   # fluid: n, momentum density p = n·u

zero_state(box::Box, nfields::Int) =
    [zeros(box.n, box.n, box.n) for _ in 1:nfields]

_wrap(n) = ([2:n; 1], [n; 1:n-1])   # ip, im index maps (periodic)

"""
    em_rhs!(dS, S, box, ip, im; Jx=nothing, Jy=nothing, Jz=nothing, gJ=1.0)

dE = ∇×B − gJ·J, dB = −∇×E into dS[EX..BZ]. J arrays are optional Yee-edge
current patterns scaled by gJ.
"""
function em_rhs!(dS, S, box::Box, ip, im;
                 Jx = nothing, Jy = nothing, Jz = nothing, gJ = 1.0)
    n = box.n
    idx = 1 / box.dx
    Ex, Ey, Ez = S[EX], S[EY], S[EZ]
    Bx, By, Bz = S[BX], S[BY], S[BZ]
    dEx, dEy, dEz = dS[EX], dS[EY], dS[EZ]
    dBx, dBy, dBz = dS[BX], dS[BY], dS[BZ]
    Threads.@threads for k in 1:n
        kp = ip[k]; km = im[k]
        @inbounds for j in 1:n
            jp = ip[j]; jm = im[j]
            for i in 1:n
                i_p = ip[i]; i_m = im[i]
                dEx[i, j, k] = (Bz[i, j, k] - Bz[i, jm, k]) * idx -
                               (By[i, j, k] - By[i, j, km]) * idx
                dEy[i, j, k] = (Bx[i, j, k] - Bx[i, j, km]) * idx -
                               (Bz[i, j, k] - Bz[i_m, j, k]) * idx
                dEz[i, j, k] = (By[i, j, k] - By[i_m, j, k]) * idx -
                               (Bx[i, j, k] - Bx[i, jm, k]) * idx
                dBx[i, j, k] = -((Ez[i, jp, k] - Ez[i, j, k]) * idx -
                                 (Ey[i, j, kp] - Ey[i, j, k]) * idx)
                dBy[i, j, k] = -((Ex[i, j, kp] - Ex[i, j, k]) * idx -
                                 (Ez[i_p, j, k] - Ez[i, j, k]) * idx)
                dBz[i, j, k] = -((Ey[i_p, j, k] - Ey[i, j, k]) * idx -
                                 (Ex[i, jp, k] - Ex[i, j, k]) * idx)
            end
        end
    end
    if Jx !== nothing
        Threads.@threads for k in 1:n
            @inbounds for j in 1:n, i in 1:n
                dEx[i, j, k] -= gJ * Jx[i, j, k]
                dEy[i, j, k] -= gJ * Jy[i, j, k]
                dEz[i, j, k] -= gJ * Jz[i, j, k]
            end
        end
    end
    dS
end

"""
    make_sponge(box, width, dt; sigma0=12/(width*box.dx)) -> mask

Per-step multiplicative damping factor exp(−σ dt), with σ ramping cubically
over `width` cells from each boundary. Returns an n³ array.
"""
function make_sponge(box::Box, width::Int, dt::Float64;
                     sigma0 = 12 / (width * box.dx))
    n = box.n
    σ1 = zeros(n)
    for i in 1:n
        d = min(i - 1, n - i)
        d < width && (σ1[i] = sigma0 * ((width - d) / width)^3)
    end
    mask = Array{Float64,3}(undef, n, n, n)
    for k in 1:n, j in 1:n, i in 1:n
        mask[i, j, k] = exp(-(σ1[i] + σ1[j] + σ1[k]) * dt)
    end
    mask
end

"""
    apply_sponge!(S, mask, box; fields=EX:BZ) -> energy removed

Damp the given fields and return the field energy removed this step.
"""
function apply_sponge!(S, mask, box::Box; fields = EX:BZ)
    removed = Threads.Atomic{Float64}(0.0)
    n = box.n
    for f in fields
        A = S[f]
        Threads.@threads for k in 1:n
            acc = 0.0
            @inbounds for j in 1:n, i in 1:n
                m = mask[i, j, k]
                a = A[i, j, k]
                acc += 0.5 * a * a * (1 - m * m)
                A[i, j, k] = a * m
            end
            Threads.atomic_add!(removed, acc)
        end
    end
    removed[] * cellvol(box)
end

field_energy(S, box::Box) =
    0.5 * cellvol(box) *
    sum(sum(abs2, S[f]) for f in EX:BZ)

"""
    div_B(S, box, ip) -> n³ array of ∇·B at cell centers (forward diffs).
"""
function div_B(S, box::Box, ip)
    n = box.n
    idx = 1 / box.dx
    Bx, By, Bz = S[BX], S[BY], S[BZ]
    D = zeros(n, n, n)
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        D[i, j, k] = (Bx[ip[i], j, k] - Bx[i, j, k]) * idx +
                     (By[i, ip[j], k] - By[i, j, k]) * idx +
                     (Bz[i, j, ip[k]] - Bz[i, j, k]) * idx
    end
    D
end

"""
    div_E(S, box, im) -> n³ array of ∇·E at nodes (backward diffs).
"""
function div_E(S, box::Box, im)
    n = box.n
    idx = 1 / box.dx
    Ex, Ey, Ez = S[EX], S[EY], S[EZ]
    D = zeros(n, n, n)
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        D[i, j, k] = (Ex[i, j, k] - Ex[im[i], j, k]) * idx +
                     (Ey[i, j, k] - Ey[i, im[j], k]) * idx +
                     (Ez[i, j, k] - Ez[i, j, im[k]]) * idx
    end
    D
end

# --- generic SSP-RK3 over a vector-of-arrays state ------------------------

function _lincomb!(Y, a, A, b, B)
    for f in eachindex(Y)
        Yf, Af, Bf = Y[f], A[f], B[f]
        Threads.@threads for idx in eachindex(Yf)
            @inbounds Yf[idx] = a * Af[idx] + b * Bf[idx]
        end
    end
end

"""
    ssprk3!(S, S1, S2, K, rhs!, dt, t)

Advance S by one SSP-RK3 step. `rhs!(K, S, t)` writes tendencies into K.
S1, S2, K are scratch states of the same shape.
"""
function ssprk3!(S, S1, S2, K, rhs!, dt, t)
    rhs!(K, S, t)
    _lincomb!(S1, 1.0, S, dt, K)                 # S1 = S + dt K
    rhs!(K, S1, t + dt)
    _lincomb!(S1, 1.0, S1, dt, K)                # S1 + dt K
    _lincomb!(S2, 0.75, S, 0.25, S1)             # S2 = ¾S + ¼(S1 + dt K)
    rhs!(K, S2, t + dt / 2)
    _lincomb!(S2, 1.0, S2, dt, K)
    _lincomb!(S, 1 / 3, S, 2 / 3, S2)            # S = ⅓S + ⅔(S2 + dt K)
    S
end
