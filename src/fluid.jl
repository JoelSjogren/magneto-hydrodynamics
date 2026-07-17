# Cold (pressureless) electron fluid coupled to the Yee EM grid — the
# dimensionless Euler–Maxwell system of README §3.2–3.3:
#   ∂n/∂t + ∇·(n u)          = 0
#   ∂u/∂t + (u·∇)u           = −(E + u×B)        (electrons, e/m → 1)
#   ∂E/∂t = ∇×B − J,  J = −n u,  ρ = f·n_b − n
# Fluid variables (n and momentum density p = n·u) live on nodes; advection
# uses a Rusanov finite-volume flux (the cold system is weakly hyperbolic and
# forms density caustics — Rusanov dissipation regularizes them). Gauss's law
# is maintained by Marder divergence cleaning.

mutable struct FluidSim
    box::Box
    ip::Vector{Int}
    im::Vector{Int}
    dt::Float64
    t::Float64
    f_neut::Float64                 # neutralization fraction f
    nb::Array{Float64,3}            # frozen ion background density n_b
    S::Vector{Array{Float64,3}}     # 10 fields: E, B, n, p
    S1::Vector{Array{Float64,3}}
    S2::Vector{Array{Float64,3}}
    K::Vector{Array{Float64,3}}
    Jx::Array{Float64,3}            # scratch: edge currents
    Jy::Array{Float64,3}
    Jz::Array{Float64,3}
    mask::Union{Nothing,Array{Float64,3}}
    absorbed::Float64               # cumulative energy removed by sponge
    alpha_floor::Float64
    n_floor::Float64
    marder::Float64                 # cleaning coefficient κ (0 disables)
end

function FluidSim(box::Box; dt = 0.3 * box.dx, f_neut = 1.0,
                  sponge_width = 0, alpha_floor = 0.02, n_floor = 1e-6,
                  marder = 0.15 * box.dx^2)
    ip, im = _wrap(box.n)
    S = zero_state(box, 10)
    fill!(S[FN], 1.0)
    mask = sponge_width > 0 ? make_sponge(box, sponge_width, dt) : nothing
    FluidSim(box, ip, im, dt, 0.0, f_neut, ones(box.n, box.n, box.n),
             S, zero_state(box, 10), zero_state(box, 10), zero_state(box, 10),
             zeros(box.n, box.n, box.n), zeros(box.n, box.n, box.n),
             zeros(box.n, box.n, box.n), mask, 0.0, alpha_floor, n_floor, marder)
end

"Edge currents J = −p averaged from nodes to Yee edge positions."
function _deposit_current!(sim::FluidSim, S)
    n = sim.box.n
    ip = sim.ip
    px, py, pz = S[FPX], S[FPY], S[FPZ]
    Jx, Jy, Jz = sim.Jx, sim.Jy, sim.Jz
    Threads.@threads for k in 1:n
        kp = ip[k]
        @inbounds for j in 1:n
            jp = ip[j]
            for i in 1:n
                Jx[i, j, k] = -0.5 * (px[i, j, k] + px[ip[i], j, k])
                Jy[i, j, k] = -0.5 * (py[i, j, k] + py[i, jp, k])
                Jz[i, j, k] = -0.5 * (pz[i, j, k] + pz[i, j, kp])
            end
        end
    end
end

"""
Rusanov flux update for the pressureless system along one axis.
`ax` = 1,2,3; adds −∂F/∂x_ax to the fluid tendencies in dS.
"""
function _advect_axis!(dS, S, sim::FluidSim, ax::Int)
    n = sim.box.n
    idx = 1 / sim.box.dx
    ip = sim.ip
    nf = sim.n_floor
    αf = sim.alpha_floor
    N = S[FN]
    P = (S[FPX], S[FPY], S[FPZ])
    dN = dS[FN]
    dP = (dS[FPX], dS[FPY], dS[FPZ])
    pd = P[ax]
    # `o` maps to an axis ≠ `ax` below, so face updates never race across
    # threads
    Threads.@threads for o in 1:n
        @inbounds for m in 1:n, i in 1:n
            # (i, m, o) mapped so that `i` runs along `ax`
            if ax == 1
                a, b, c = i, m, o
                a2, b2, c2 = ip[i], m, o
            elseif ax == 2
                a, b, c = m, i, o
                a2, b2, c2 = m, ip[i], o
            else
                a, b, c = m, o, i
                a2, b2, c2 = m, o, ip[i]
            end
            nL = max(N[a, b, c], nf); nR = max(N[a2, b2, c2], nf)
            uL = pd[a, b, c] / nL;    uR = pd[a2, b2, c2] / nR
            α = max(abs(uL), abs(uR)) + αf
            # face flux F̂ = ½(F_L+F_R) − ½α(U_R−U_L), F = u_ax·U
            fN = 0.5 * (uL * N[a, b, c] + uR * N[a2, b2, c2]) -
                 0.5 * α * (N[a2, b2, c2] - N[a, b, c])
            dN[a, b, c] -= fN * idx
            dN[a2, b2, c2] += fN * idx
            for q in 1:3
                Pq = P[q]
                fP = 0.5 * (uL * Pq[a, b, c] + uR * Pq[a2, b2, c2]) -
                     0.5 * α * (Pq[a2, b2, c2] - Pq[a, b, c])
                dPq = dP[q]
                dPq[a, b, c] -= fP * idx
                dPq[a2, b2, c2] += fP * idx
            end
        end
    end
end

"Lorentz source −(n E + p × B) with E, B averaged to nodes."
function _lorentz!(dS, S, sim::FluidSim)
    n = sim.box.n
    im = sim.im
    Ex, Ey, Ez = S[EX], S[EY], S[EZ]
    Bx, By, Bz = S[BX], S[BY], S[BZ]
    N, px, py, pz = S[FN], S[FPX], S[FPY], S[FPZ]
    dpx, dpy, dpz = dS[FPX], dS[FPY], dS[FPZ]
    Threads.@threads for k in 1:n
        km = im[k]
        @inbounds for j in 1:n
            jm = im[j]
            for i in 1:n
                i_m = im[i]
                ex = 0.5 * (Ex[i, j, k] + Ex[i_m, j, k])
                ey = 0.5 * (Ey[i, j, k] + Ey[i, jm, k])
                ez = 0.5 * (Ez[i, j, k] + Ez[i, j, km])
                bx = 0.25 * (Bx[i, j, k] + Bx[i, jm, k] +
                             Bx[i, j, km] + Bx[i, jm, km])
                by = 0.25 * (By[i, j, k] + By[i_m, j, k] +
                             By[i, j, km] + By[i_m, j, km])
                bz = 0.25 * (Bz[i, j, k] + Bz[i_m, j, k] +
                             Bz[i, jm, k] + Bz[i_m, jm, k])
                nn = N[i, j, k]
                dpx[i, j, k] -= nn * ex + (py[i, j, k] * bz - pz[i, j, k] * by)
                dpy[i, j, k] -= nn * ey + (pz[i, j, k] * bx - px[i, j, k] * bz)
                dpz[i, j, k] -= nn * ez + (px[i, j, k] * by - py[i, j, k] * bx)
            end
        end
    end
end

function coupled_rhs!(K, S, sim::FluidSim)
    _deposit_current!(sim, S)
    em_rhs!(K, S, sim.box, sim.ip, sim.im; Jx = sim.Jx, Jy = sim.Jy, Jz = sim.Jz)
    for f in FN:FPZ
        fill!(K[f], 0.0)
    end
    for ax in 1:3
        _advect_axis!(K, S, sim, ax)
    end
    _lorentz!(K, S, sim)
    K
end

"Marder divergence cleaning: E += κ ∇(∇·E − ρ)."
function marder_clean!(sim::FluidSim)
    sim.marder == 0 && return
    box = sim.box
    n = box.n
    idx = 1 / box.dx
    D = div_E(sim.S, box, sim.im)
    N = sim.S[FN]
    nb = sim.nb
    f = sim.f_neut
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        D[i, j, k] -= f * nb[i, j, k] - N[i, j, k]
    end
    κ = sim.marder
    Ex, Ey, Ez = sim.S[EX], sim.S[EY], sim.S[EZ]
    ip = sim.ip
    Threads.@threads for k in 1:n
        kp = ip[k]
        @inbounds for j in 1:n
            jp = ip[j]
            for i in 1:n
                Ex[i, j, k] += κ * (D[ip[i], j, k] - D[i, j, k]) * idx
                Ey[i, j, k] += κ * (D[i, jp, k] - D[i, j, k]) * idx
                Ez[i, j, k] += κ * (D[i, j, kp] - D[i, j, k]) * idx
            end
        end
    end
end

function step!(sim::FluidSim)
    rhs! = (K, S, t) -> coupled_rhs!(K, S, sim)
    ssprk3!(sim.S, sim.S1, sim.S2, sim.K, rhs!, sim.dt, sim.t)
    if sim.mask !== nothing
        sim.absorbed += apply_sponge!(sim.S, sim.mask, sim.box)
    end
    marder_clean!(sim)
    sim.t += sim.dt
    sim
end

kinetic_energy(sim::FluidSim) = kinetic_energy(sim.S, sim.box, sim.n_floor)

function kinetic_energy(S, box::Box, n_floor)
    N, px, py, pz = S[FN], S[FPX], S[FPY], S[FPZ]
    e = 0.0
    @inbounds for idx in eachindex(N)
        nn = max(N[idx], n_floor)
        e += (px[idx]^2 + py[idx]^2 + pz[idx]^2) / nn
    end
    0.5 * e * cellvol(box)
end

"RMS Gauss-law residual ‖∇·E − ρ‖ (diagnostic)."
function gauss_residual(sim::FluidSim)
    D = div_E(sim.S, sim.box, sim.im)
    N = sim.S[FN]
    r = 0.0
    @inbounds for idx in eachindex(D)
        r += (D[idx] - (sim.f_neut * sim.nb[idx] - N[idx]))^2
    end
    sqrt(r / length(D))
end
