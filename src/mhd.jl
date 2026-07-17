# v2 — single-fluid, quasineutral, compressible, visco-resistive MHD in
# Alfvén units (README §7.3). Isothermal EOS p = cs²ρ. Cell-centered
# finite volume, Rusanov fluxes, GLM/Dedner hyperbolic divergence cleaning,
# explicit resistivity, adaptive dt, SSP-RK3 via the shared harness.
#
# State layout (Vector{Array{Float64,3}}):
#   MRHO ρ | MMX..MMZ momentum ρv | MBX..MBZ B | MPSI ψ (GLM scalar)

const MRHO, MMX, MMY, MMZ, MBX, MBY, MBZ, MPSI = 1, 2, 3, 4, 5, 6, 7, 8

mutable struct MHDSim
    box::Box
    ip::Vector{Int}
    im::Vector{Int}
    t::Float64
    dt::Float64                     # last dt used (adaptive)
    cfl::Float64
    cs::Float64                     # isothermal sound speed (units of v_A0)
    eta::Float64                    # resistivity = 1/Lundquist
    ch::Float64                     # GLM wave speed (updated each step)
    rho_floor::Float64
    S::Vector{Array{Float64,3}}
    S1::Vector{Array{Float64,3}}
    S2::Vector{Array{Float64,3}}
    K::Vector{Array{Float64,3}}
    mask::Union{Nothing,Array{Float64,3}}
    absorbed::Float64
end

function MHDSim(box::Box; cs = 0.5, eta = 2e-3, cfl = 0.3,
                sponge_width = 0, rho_floor = 1e-3)
    ip, im = _wrap(box.n)
    S = zero_state(box, 8)
    fill!(S[MRHO], 1.0)
    # mask is rebuilt against the current dt each step (dt is adaptive), so
    # store the sponge σ-profile via a fixed reference dt of 1 and
    # exponentiate: mask_dt = mask_ref^dt
    mask = sponge_width > 0 ? make_sponge(box, sponge_width, 1.0) : nothing
    MHDSim(box, ip, im, 0.0, 0.0, cfl, cs, eta, 1.0, rho_floor,
           S, zero_state(box, 8), zero_state(box, 8), zero_state(box, 8),
           mask, 0.0)
end

"Maximum signal speed max(|v_d|) + fast speed, over the grid."
function max_speed(sim::MHDSim, S)
    ρ = S[MRHO]
    vmax = Threads.Atomic{Float64}(0.0)
    n = sim.box.n
    cs2 = sim.cs^2
    Threads.@threads for k in 1:n
        local m = 0.0
        @inbounds for j in 1:n, i in 1:n
            r = max(ρ[i, j, k], sim.rho_floor)
            vv = max(abs(S[MMX][i, j, k]), abs(S[MMY][i, j, k]),
                     abs(S[MMZ][i, j, k])) / r
            b2 = S[MBX][i, j, k]^2 + S[MBY][i, j, k]^2 + S[MBZ][i, j, k]^2
            m = max(m, vv + sqrt(cs2 + b2 / r))
        end
        Threads.atomic_max!(vmax, m)
    end
    vmax[]
end

@inline _minmod(a, b) = a * b <= 0 ? 0.0 : (abs(a) < abs(b) ? a : b)

"Flux vector of the GLM-MHD system from a state 8-vector, along axis ax."
@inline function _flux8!(F, U, ax, cs2, ch2, fl)
    ρ = max(U[MRHO], fl)
    vd = U[MMX+ax-1] / ρ
    Bd = U[MBX+ax-1]
    b2 = U[MBX]^2 + U[MBY]^2 + U[MBZ]^2
    pt = cs2 * ρ + 0.5 * b2
    F[MRHO] = U[MMX+ax-1]
    for q in 0:2
        F[MMX+q] = U[MMX+q] * vd - U[MBX+q] * Bd
        F[MBX+q] = U[MBX+q] * vd - (U[MMX+q] / ρ) * Bd
    end
    F[MMX+ax-1] += pt
    F[MBX+ax-1] = U[MPSI]            # GLM
    F[MPSI] = ch2 * Bd
    abs(vd) + sqrt(cs2 + b2 / ρ)     # returns signal speed of this state
end

# MUSCL (minmod-limited linear reconstruction) + Rusanov flux along one
# axis for the full 8-component system. Second-order in smooth regions —
# first-order Rusanov destroyed few-cell ring cores within an Alfvén time.
function _mhd_flux_axis!(K, S, sim::MHDSim, ax::Int)
    n = sim.box.n
    idx = 1 / sim.box.dx
    ip = sim.ip; im = sim.im
    cs2 = sim.cs^2
    ch = sim.ch
    ch2 = ch * ch
    fl = sim.rho_floor
    Threads.@threads for o in 1:n
        FL = zeros(8); FR = zeros(8); UL = zeros(8); UR = zeros(8)
        @inbounds for m2 in 1:n, i in 1:n
            i2 = ip[i]
            if ax == 1
                a, b, c = i, m2, o;      a2, b2, c2 = i2, m2, o
                am, bm, cm = im[i], m2, o; a3, b3, c3 = ip[i2], m2, o
            elseif ax == 2
                a, b, c = m2, i, o;      a2, b2, c2 = m2, i2, o
                am, bm, cm = m2, im[i], o; a3, b3, c3 = m2, ip[i2], o
            else
                a, b, c = m2, o, i;      a2, b2, c2 = m2, o, i2
                am, bm, cm = m2, o, im[i]; a3, b3, c3 = m2, o, ip[i2]
            end
            for q in 1:8
                u0 = S[q][am, bm, cm]; u1 = S[q][a, b, c]
                u2 = S[q][a2, b2, c2]; u3 = S[q][a3, b3, c3]
                UL[q] = u1 + 0.5 * _minmod(u1 - u0, u2 - u1)
                UR[q] = u2 - 0.5 * _minmod(u2 - u1, u3 - u2)
            end
            sL = _flux8!(FL, UL, ax, cs2, ch2, fl)
            sR = _flux8!(FR, UR, ax, cs2, ch2, fl)
            α = max(sL, sR, ch)
            for q in 1:8
                f = 0.5 * (FL[q] + FR[q]) - 0.5 * α * (UR[q] - UL[q])
                K[q][a, b, c] -= f * idx
                K[q][a2, b2, c2] += f * idx
            end
        end
    end
end

"Explicit resistive term η∇²B."
function _resistive!(K, S, sim::MHDSim)
    sim.eta == 0 && return
    n = sim.box.n
    ηi2 = sim.eta / sim.box.dx^2
    ip = sim.ip; im = sim.im
    for q in MBX:MBZ
        B = S[q]; dB = K[q]
        Threads.@threads for k in 1:n
            kp = ip[k]; km = im[k]
            @inbounds for j in 1:n
                jp = ip[j]; jm = im[j]
                for i in 1:n
                    dB[i, j, k] += ηi2 * (B[ip[i], j, k] + B[im[i], j, k] +
                                          B[i, jp, k] + B[i, jm, k] +
                                          B[i, j, kp] + B[i, j, km] -
                                          6 * B[i, j, k])
                end
            end
        end
    end
end

function mhd_rhs!(K, S, sim::MHDSim)
    for q in 1:8
        fill!(K[q], 0.0)
    end
    for ax in 1:3
        _mhd_flux_axis!(K, S, sim, ax)
    end
    _resistive!(K, S, sim)
    K
end

function mhd_step!(sim::MHDSim)
    smax = max(1e-10, max_speed(sim, sim.S))
    sim.ch = smax
    dt = sim.cfl * sim.box.dx / smax
    sim.dt = dt
    rhs! = (K, S, t) -> mhd_rhs!(K, S, sim)
    ssprk3!(sim.S, sim.S1, sim.S2, sim.K, rhs!, dt, sim.t)
    # GLM damping (Dedner mixed correction) and density floor
    damp = exp(-0.1 * sim.ch * dt / sim.box.dx)
    ψ = sim.S[MPSI]; ρ = sim.S[MRHO]
    Threads.@threads for idx in eachindex(ψ)
        @inbounds ψ[idx] *= damp
        @inbounds ρ[idx] < sim.rho_floor && (ρ[idx] = sim.rho_floor)
    end
    if sim.mask !== nothing
        # mask was built for dt=1; raise to current dt: m^dt
        n = sim.box.n
        removed = 0.0
        for q in (MMX, MMY, MMZ, MBX, MBY, MBZ, MPSI)
            A = sim.S[q]
            @inbounds for k in 1:n, j in 1:n, i in 1:n
                m = sim.mask[i, j, k]^dt
                a = A[i, j, k]
                removed += 0.5 * a * a * (1 - m * m)
                A[i, j, k] = a * m
            end
        end
        sim.absorbed += removed * cellvol(sim.box)
    end
    sim.t += dt
    sim
end

mhd_kinetic_energy(sim::MHDSim) = begin
    ρ = sim.S[MRHO]
    e = 0.0
    @inbounds for idx in eachindex(ρ)
        e += (sim.S[MMX][idx]^2 + sim.S[MMY][idx]^2 + sim.S[MMZ][idx]^2) /
             max(ρ[idx], sim.rho_floor)
    end
    0.5 * e * cellvol(sim.box)
end

mhd_magnetic_energy(sim::MHDSim) =
    0.5 * cellvol(sim.box) *
    (sum(abs2, sim.S[MBX]) + sum(abs2, sim.S[MBY]) + sum(abs2, sim.S[MBZ]))

"Central-difference curl of fields (fx,fy,fz) — used for J = ∇×B and ω = ∇×v."
function curl_central(Fx, Fy, Fz, box::Box, ip, im)
    n = box.n
    h = 1 / (2 * box.dx)
    Cx = zeros(n, n, n); Cy = zeros(n, n, n); Cz = zeros(n, n, n)
    Threads.@threads for k in 1:n
        kp = ip[k]; km = im[k]
        @inbounds for j in 1:n
            jp = ip[j]; jm = im[j]
            for i in 1:n
                i_p = ip[i]; i_m = im[i]
                Cx[i, j, k] = (Fz[i, jp, k] - Fz[i, jm, k]) * h -
                              (Fy[i, j, kp] - Fy[i, j, km]) * h
                Cy[i, j, k] = (Fx[i, j, kp] - Fx[i, j, km]) * h -
                              (Fz[i_p, j, k] - Fz[i_m, j, k]) * h
                Cz[i, j, k] = (Fy[i_p, j, k] - Fy[i_m, j, k]) * h -
                              (Fx[i, jp, k] - Fx[i, jm, k]) * h
            end
        end
    end
    Cx, Cy, Cz
end

"Dipole and anapole moments of the grid current J = ∇×B."
function grid_moments(sim::MHDSim)
    Jx, Jy, Jz = curl_central(sim.S[MBX], sim.S[MBY], sim.S[MBZ],
                              sim.box, sim.ip, sim.im)
    box = sim.box
    n = box.n
    m = zeros(3); T = zeros(3)
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        x = center(box, i); y = center(box, j); z = center(box, k)
        jx = Jx[i, j, k]; jy = Jy[i, j, k]; jz = Jz[i, j, k]
        m[1] += y * jz - z * jy
        m[2] += z * jx - x * jz
        m[3] += x * jy - y * jx
        rj = x * jx + y * jy + z * jz
        r2 = x * x + y * y + z * z
        T[1] += rj * x - 2 * r2 * jx
        T[2] += rj * y - 2 * r2 * jy
        T[3] += rj * z - 2 * r2 * jz
    end
    dV = cellvol(box)
    (dV / 2) .* m, (dV / 10) .* T
end

"""
    azimuthal_spectrum(A, box, radius; z=0.0, nmodes=12, nphi=256)

|DFT| amplitudes of the field A sampled on the circle of given radius in the
plane z — mode 0..nmodes. Detects azimuthal breakup ("ring of rings").
"""
function azimuthal_spectrum(A, box::Box, radius; z = 0.0, nmodes = 12,
                            nphi = 256)
    samp = Vector{Float64}(undef, nphi)
    for q in 1:nphi
        φ = 2π * (q - 1) / nphi
        x = radius * cos(φ); y = radius * sin(φ)
        # trilinear interpolation at (x, y, z), cell-centered grid
        fx = (x - box.lo) / box.dx - 0.5; fy = (y - box.lo) / box.dx - 0.5
        fz = (z - box.lo) / box.dx - 0.5
        i0 = clamp(floor(Int, fx) + 1, 1, box.n - 1)
        j0 = clamp(floor(Int, fy) + 1, 1, box.n - 1)
        k0 = clamp(floor(Int, fz) + 1, 1, box.n - 1)
        tx = clamp(fx - (i0 - 1), 0.0, 1.0)
        ty = clamp(fy - (j0 - 1), 0.0, 1.0)
        tz = clamp(fz - (k0 - 1), 0.0, 1.0)
        s = 0.0
        for dk in 0:1, dj in 0:1, di in 0:1
            w = (di == 1 ? tx : 1 - tx) * (dj == 1 ? ty : 1 - ty) *
                (dk == 1 ? tz : 1 - tz)
            s += w * A[i0+di, j0+dj, k0+dk]
        end
        samp[q] = s
    end
    amps = Vector{Float64}(undef, nmodes + 1)
    for m in 0:nmodes
        cr = 0.0; ci = 0.0
        for q in 1:nphi
            φ = 2π * (q - 1) / nphi
            cr += samp[q] * cos(m * φ)
            ci += samp[q] * sin(m * φ)
        end
        amps[m+1] = sqrt(cr^2 + ci^2) / nphi
    end
    amps
end

# --- initial-condition builders (cell-centered, curl of a potential so the
# --- seeded fields are discretely divergence-free under central differences)

"Gaussian ring profile exp(−s²/2a²), s = distance to the circle (R, z0)."
@inline function _ring_s2(x, y, z, R, z0)
    rc = sqrt(x * x + y * y)
    (rc - R)^2 + (z - z0)^2
end

"""
    add_flux_ring!(sim; R=1, a=0.25, z0=0, A0=1.0, Bt0=0.0)

Add a magnetic flux ring: poloidal field from the curl of a toroidal vector
potential A_φ = A0·exp(−s²/2a²) (ring current sense = sign(A0)), plus an
optional direct toroidal component B_φ = Bt0·exp(−s²/2a²) (twist/helicity
sign = sign(Bt0) relative to A0).
"""
function add_flux_ring!(sim::MHDSim; R = 1.0, a = 0.25, z0 = 0.0,
                        A0 = 1.0, Bt0 = 0.0)
    box = sim.box
    n = box.n
    Ax = zeros(n, n, n); Ay = zeros(n, n, n); Az = zeros(n, n, n)
    inv2a2 = 1 / (2a^2)
    for k in 1:n, j in 1:n, i in 1:n
        x = center(box, i); y = center(box, j); z = center(box, k)
        rc = max(sqrt(x * x + y * y), 1e-12)
        f = A0 * exp(-_ring_s2(x, y, z, R, z0) * inv2a2)
        Ax[i, j, k] = -f * y / rc
        Ay[i, j, k] = f * x / rc
    end
    if Bt0 != 0
        # toroidal component from a potential too, so the discrete central
        # divergence stays at machine zero: B_φ = −∂A_z/∂r_c with
        # A_z = −Bt0·e^{−(z−z0)²/2a²}·∫₀^{r_c} e^{−(r−R)²/2a²} dr
        pre = Bt0 * a * sqrt(π / 2)
        s2 = a * sqrt(2.0)
        for k in 1:n, j in 1:n, i in 1:n
            x = center(box, i); y = center(box, j); z = center(box, k)
            rc = sqrt(x * x + y * y)
            Az[i, j, k] -= pre * exp(-(z - z0)^2 * inv2a2) *
                           (_erf((rc - R) / s2) + _erf(R / s2))
        end
    end
    Bx, By, Bz = curl_central(Ax, Ay, Az, box, sim.ip, sim.im)
    sim.S[MBX] .+= Bx
    sim.S[MBY] .+= By
    sim.S[MBZ] .+= Bz
    sim
end

"Abramowitz–Stegun 7.1.26 rational approximation of erf (|err| < 1.5e-7)."
function _erf(x)
    s = x < 0 ? -1.0 : 1.0
    x = abs(x)
    t = 1 / (1 + 0.3275911x)
    y = 1 - (((((1.061405429t - 1.453152027)t) + 1.421413741)t -
              0.284496736)t + 0.254829592)t * exp(-x * x)
    s * y
end

"""
    add_vortex_ring!(sim; R=1, a=0.25, z0=0, P0=1.0)

Add a hydrodynamic vortex ring: v = ∇×(ψ_φ e_φ), ψ_φ = P0·exp(−s²/2a²).
P0 > 0 propagates toward +z. Momentum is ρ·v with the current ρ.
"""
function add_vortex_ring!(sim::MHDSim; R = 1.0, a = 0.25, z0 = 0.0, P0 = 1.0)
    box = sim.box
    n = box.n
    Px = zeros(n, n, n); Py = zeros(n, n, n); Pz = zeros(n, n, n)
    inv2a2 = 1 / (2a^2)
    for k in 1:n, j in 1:n, i in 1:n
        x = center(box, i); y = center(box, j); z = center(box, k)
        rc = max(sqrt(x * x + y * y), 1e-12)
        f = P0 * exp(-_ring_s2(x, y, z, R, z0) * inv2a2)
        Px[i, j, k] = -f * y / rc
        Py[i, j, k] = f * x / rc
    end
    vx, vy, vz = curl_central(Px, Py, Pz, box, sim.ip, sim.im)
    ρ = sim.S[MRHO]
    @inbounds for idx in eachindex(ρ)
        sim.S[MMX][idx] += ρ[idx] * vx[idx]
        sim.S[MMY][idx] += ρ[idx] * vy[idx]
        sim.S[MMZ][idx] += ρ[idx] * vz[idx]
    end
    sim
end
