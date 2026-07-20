# v3 — potentials-first, coherent matter (README §8). A Gross–Pitaevskii
# condensate minimally coupled to the electromagnetic potentials; E and B
# are derived diagnostics only, never dynamical inputs:
#
#   i ∂t ψ = [ (−i∇ − qA)² / 2m + qφ + g|ψ|² ] ψ
#   (∂tt − ∇²) φ = ρ = q (|ψ|² − n̄)
#   (∂tt − ∇²) A = J = (q/m) [ Im(ψ*∇ψ) − q|ψ|² A ]
#
# Lorenz gauge: the wave equations preserve ∂tφ + ∇·A = 0 exactly when
# charge continuity holds, which the GPE guarantees analytically (and the
# central-difference discretization approximately). n̄ is a uniform
# neutralizing background (jellium) so a periodic box carries no net
# charge. Units c = μ0 = ε0 = ℏ = 1: m sets the dispersion, q²n/m the
# plasma frequency — the same normalization family as v1.
#
# State (10 real cell-centered fields):
#   GR ψ_re | GI ψ_im | GPH φ | GPP ∂tφ | GAX..GAZ A | GQX..GQZ ∂tA
#
# Explicit method of lines on the shared SSP-RK3 harness. Fixed dt from
# the Schrödinger stability bound |λ|dt ≤ √3 with |λ| ≈ 6/(m dx²)
# (dominates the wave CFL dt ≤ dx at any resolved dx).

const GR, GI, GPH, GPP, GAX, GAY, GAZ, GQX, GQY, GQZ =
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10

mutable struct GPESim
    box::Box
    ip::Vector{Int}
    im::Vector{Int}
    t::Float64
    dt::Float64
    m::Float64                      # particle mass
    q::Float64                      # particle charge
    g::Float64                      # contact interaction strength
    nbar::Float64                   # neutralizing background density
    S::Vector{Array{Float64,3}}
    S1::Vector{Array{Float64,3}}
    S2::Vector{Array{Float64,3}}
    K::Vector{Array{Float64,3}}
end

function GPESim(box::Box; m = 1.0, q = 0.0, g = 0.0,
                cfl_s = 0.25, cfl_w = 0.5)
    ip, im_ = _wrap(box.n)
    dt = min(cfl_s * m * box.dx^2, cfl_w * box.dx)
    GPESim(box, ip, im_, 0.0, dt, m, q, g, 0.0,
           zero_state(box, 10), zero_state(box, 10),
           zero_state(box, 10), zero_state(box, 10))
end

"Set the jellium background to the current mean density (call after ICs)."
function gpe_neutralize!(sim::GPESim)
    sim.nbar = (sum(abs2, sim.S[GR]) + sum(abs2, sim.S[GI])) / sim.box.n^3
    sim
end

function gpe_rhs!(K, S, sim::GPESim)
    n = sim.box.n
    ih2 = 1 / sim.box.dx^2
    h = 1 / (2 * sim.box.dx)
    i2m = 1 / (2 * sim.m)
    q = sim.q
    g = sim.g
    qm = q / sim.m
    nbar = sim.nbar
    a = S[GR]; b = S[GI]; φ = S[GPH]
    Ax = S[GAX]; Ay = S[GAY]; Az = S[GAZ]
    Threads.@threads for k in 1:n
        kp = sim.ip[k]; km = sim.im[k]
        @inbounds for j in 1:n
            jp = sim.ip[j]; jm = sim.im[j]
            for i in 1:n
                i_p = sim.ip[i]; i_m = sim.im[i]
                lap(F) = (F[i_p, j, k] + F[i_m, j, k] + F[i, jp, k] +
                          F[i, jm, k] + F[i, j, kp] + F[i, j, km] -
                          6 * F[i, j, k]) * ih2
                dx(F) = (F[i_p, j, k] - F[i_m, j, k]) * h
                dy(F) = (F[i, jp, k] - F[i, jm, k]) * h
                dz(F) = (F[i, j, kp] - F[i, j, km]) * h

                ac = a[i, j, k]; bc = b[i, j, k]
                n2 = ac * ac + bc * bc
                V = q * φ[i, j, k] + g * n2

                if q == 0
                    Hr = -i2m * lap(a) + V * ac
                    Hi = -i2m * lap(b) + V * bc
                    K[GR][i, j, k] = Hi
                    K[GI][i, j, k] = -Hr
                else
                    axc = Ax[i, j, k]; ayc = Ay[i, j, k]; azc = Az[i, j, k]
                    divA = dx(Ax) + dy(Ay) + dz(Az)
                    A2 = axc * axc + ayc * ayc + azc * azc
                    dax = dx(a); day = dy(a); daz = dz(a)
                    dbx = dx(b); dby = dy(b); dbz = dz(b)
                    Adot_a = axc * dax + ayc * day + azc * daz
                    Adot_b = axc * dbx + ayc * dby + azc * dbz
                    Hr = i2m * (-lap(a) - q * divA * bc - 2q * Adot_b +
                                q * q * A2 * ac) + V * ac
                    Hi = i2m * (-lap(b) + q * divA * ac + 2q * Adot_a +
                                q * q * A2 * bc) + V * bc
                    K[GR][i, j, k] = Hi
                    K[GI][i, j, k] = -Hr
                    # EM sector: □φ = ρ, □A = J
                    K[GPH][i, j, k] = S[GPP][i, j, k]
                    K[GPP][i, j, k] = lap(φ) + q * (n2 - nbar)
                    K[GAX][i, j, k] = S[GQX][i, j, k]
                    K[GAY][i, j, k] = S[GQY][i, j, k]
                    K[GAZ][i, j, k] = S[GQZ][i, j, k]
                    K[GQX][i, j, k] = lap(Ax) +
                                      qm * (ac * dbx - bc * dax - q * n2 * axc)
                    K[GQY][i, j, k] = lap(Ay) +
                                      qm * (ac * dby - bc * day - q * n2 * ayc)
                    K[GQZ][i, j, k] = lap(Az) +
                                      qm * (ac * dbz - bc * daz - q * n2 * azc)
                end
            end
        end
    end
    if sim.q == 0
        for f in GPH:GQZ
            fill!(K[f], 0.0)
        end
    end
    K
end

function gpe_step!(sim::GPESim)
    rhs! = (K, S, t) -> gpe_rhs!(K, S, sim)
    ssprk3!(sim.S, sim.S1, sim.S2, sim.K, rhs!, sim.dt, sim.t)
    sim.t += sim.dt
    sim
end

gpe_norm(sim::GPESim) =
    (sum(abs2, sim.S[GR]) + sum(abs2, sim.S[GI])) * cellvol(sim.box)

"""
    gpe_energy(sim) -> (; kinetic, interaction, field)

Kinetic (1/2m)∫|(−i∇−qA)ψ|², interaction (g/2)∫|ψ|⁴, and EM field energy
(1/2)∫(E² + B²) with E = −∇φ − ∂tA, B = ∇×A.
"""
function gpe_energy(sim::GPESim)
    box = sim.box
    n = box.n
    h = 1 / (2 * box.dx)
    q = sim.q
    a = sim.S[GR]; b = sim.S[GI]
    ekin = 0.0; eint = 0.0
    ip = sim.ip; im_ = sim.im
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        kp = ip[k]; km = im_[k]; jp = ip[j]; jm = im_[j]
        i_p = ip[i]; i_m = im_[i]
        for (da, db, Aq) in (((a[i_p, j, k] - a[i_m, j, k]) * h,
                              (b[i_p, j, k] - b[i_m, j, k]) * h,
                              sim.S[GAX][i, j, k]),
                             ((a[i, jp, k] - a[i, jm, k]) * h,
                              (b[i, jp, k] - b[i, jm, k]) * h,
                              sim.S[GAY][i, j, k]),
                             ((a[i, j, kp] - a[i, j, km]) * h,
                              (b[i, j, kp] - b[i, j, km]) * h,
                              sim.S[GAZ][i, j, k]))
            re = db - q * Aq * a[i, j, k]
            im2 = -da - q * Aq * b[i, j, k]
            ekin += re * re + im2 * im2
        end
        eint += (a[i, j, k]^2 + b[i, j, k]^2)^2
    end
    efld = 0.0
    if q != 0
        Bx, By, Bz = curl_central(sim.S[GAX], sim.S[GAY], sim.S[GAZ],
                                  box, ip, im_)
        @inbounds for k in 1:n, j in 1:n, i in 1:n
            kp = ip[k]; km = im_[k]; jp = ip[j]; jm = im_[j]
            i_p = ip[i]; i_m = im_[i]
            Ex = -(sim.S[GPH][i_p, j, k] - sim.S[GPH][i_m, j, k]) * h -
                 sim.S[GQX][i, j, k]
            Ey = -(sim.S[GPH][i, jp, k] - sim.S[GPH][i, jm, k]) * h -
                 sim.S[GQY][i, j, k]
            Ez = -(sim.S[GPH][i, j, kp] - sim.S[GPH][i, j, km]) * h -
                 sim.S[GQZ][i, j, k]
            efld += Ex^2 + Ey^2 + Ez^2 +
                    Bx[i, j, k]^2 + By[i, j, k]^2 + Bz[i, j, k]^2
        end
    end
    dV = cellvol(box)
    (kinetic = ekin * dV / (2 * sim.m), interaction = 0.5 * sim.g * eint * dV,
     field = 0.5 * efld * dV)
end

"Uniform condensate ψ = √n0 (and matching jellium background)."
function gpe_uniform!(sim::GPESim; n0 = 1.0)
    fill!(sim.S[GR], sqrt(n0))
    fill!(sim.S[GI], 0.0)
    gpe_neutralize!(sim)
end

"Gaussian wave packet: |ψ|² has std `sigma` per axis; carrier momentum k0."
function gpe_packet!(sim::GPESim; x0 = (0.0, 0.0, 0.0), k0 = (0.0, 0.0, 0.0),
                     sigma = 1.0, amp = 1.0)
    box = sim.box
    inv4s2 = 1 / (4 * sigma^2)
    for k in 1:box.n, j in 1:box.n, i in 1:box.n
        x = center(box, i) - x0[1]
        y = center(box, j) - x0[2]
        z = center(box, k) - x0[3]
        env = amp * exp(-(x * x + y * y + z * z) * inv4s2)
        ph = k0[1] * x + k0[2] * y + k0[3] * z
        sim.S[GR][i, j, k] = env * cos(ph)
        sim.S[GI][i, j, k] = env * sin(ph)
    end
    gpe_neutralize!(sim)
end

"""
    gpe_vortex_pair!(sim; n0=1.0, d=4.0, xi=1.0)

Straight vortex–antivortex line pair along z at x = ±d/2 (net circulation
zero, so the phase is periodic-compatible); tanh(r/ξ) core profile.
"""
function gpe_vortex_pair!(sim::GPESim; n0 = 1.0, d = 4.0, xi = 1.0)
    box = sim.box
    for k in 1:box.n, j in 1:box.n, i in 1:box.n
        x = center(box, i); y = center(box, j)
        r1 = sqrt((x - d / 2)^2 + y * y)
        r2 = sqrt((x + d / 2)^2 + y * y)
        θ = atan(y, x - d / 2) - atan(y, x + d / 2)
        f = sqrt(n0) * tanh(r1 / xi) * tanh(r2 / xi)
        sim.S[GR][i, j, k] = f * cos(θ)
        sim.S[GI][i, j, k] = f * sin(θ)
    end
    gpe_neutralize!(sim)
end

"""
    gpe_winding(sim, xc, yc; radius=1.0, nseg=64) -> Int

Phase winding number of ψ around a circle (xc, yc) in the z-midplane.
"""
function gpe_winding(sim::GPESim, xc, yc; radius = 1.0, nseg = 64)
    box = sim.box
    k0 = box.n ÷ 2
    prev = 0.0
    tot = 0.0
    for s in 0:nseg
        θ = 2π * s / nseg
        x = xc + radius * cos(θ); y = yc + radius * sin(θ)
        fi = clamp((x - box.lo) / box.dx + 0.5, 1.0, Float64(box.n))
        fj = clamp((y - box.lo) / box.dx + 0.5, 1.0, Float64(box.n))
        i = round(Int, fi); j = round(Int, fj)
        ph = atan(sim.S[GI][i, j, k0], sim.S[GR][i, j, k0])
        if s > 0
            dp = ph - prev
            dp > π && (dp -= 2π)
            dp < -π && (dp += 2π)
            tot += dp
        end
        prev = ph
    end
    round(Int, tot / 2π)
end
