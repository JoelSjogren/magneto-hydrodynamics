# v2 MHD validation (README §7.6). The scheme is first-order Rusanov, so
# k≈grid-scale waves damp numerically; tests check propagation physics and
# the *resistivity-dependent increment* of damping, not absolute rates.

function _cp_alfven_sim(; eta = 0.0, n = 32)
    box = Box(n, π)                  # z ∈ [−π, π), k = 1, v_A = 1
    sim = MHDSim(box; cs = 0.5, eta = eta)
    ε = 0.01
    for k in 1:box.n, j in 1:box.n, i in 1:box.n
        z = center(box, k)
        sim.S[MBZ][i, j, k] = 1.0
        sim.S[MMX][i, j, k] = ε * cos(z)
        sim.S[MMY][i, j, k] = ε * sin(z)
        sim.S[MBX][i, j, k] = -ε * cos(z)   # δb = −δv → +z propagation
        sim.S[MBY][i, j, k] = -ε * sin(z)
    end
    sim, ε
end

@testset "MHD: circularly polarized Alfvén wave" begin
    sim, ε = _cp_alfven_sim()
    box = sim.box
    e0 = mhd_kinetic_energy(sim)
    while sim.t < π                  # half period: v → −A·cos z, A > 0
        mhd_step!(sim)
    end
    # project v_x(z) onto cos z and sin z: forward propagation at v_A = 1
    # means the cos-projection is negative and dominates the sin one
    pc = 0.0; ps = 0.0
    for k in 1:box.n
        z = center(box, k)
        pc += sim.S[MMX][1, 1, k] * cos(z)
        ps += sim.S[MMX][1, 1, k] * sin(z)
    end
    pc /= box.n / 2 * ε; ps /= box.n / 2 * ε   # ≈ −A/ε and phase error
    @test pc < -0.3                  # arrived with the right sign (speed ok)
    @test abs(ps) < 0.5 * abs(pc)    # phase error < ~27°
    # diffusive but propagating: energy within the 1st-order-scheme window
    @test 0.1 < mhd_kinetic_energy(sim) / e0 < 0.95
end

@testset "MHD: resistive damping increment" begin
    # identical waves with η = 0 and η = 0.05: the damping-rate difference
    # must be ≈ ηk² = 0.05 (numerical diffusion cancels in the ratio).
    amp(sim) = maximum(abs, sim.S[MMX])
    rates = Float64[]
    for η in (0.0, 0.05)
        sim, ε = _cp_alfven_sim(; eta = η)
        a0 = amp(sim)
        while sim.t < 8.0
            mhd_step!(sim)
        end
        push!(rates, -log(amp(sim) / a0) / sim.t)
    end
    dγ = rates[2] - rates[1]
    @test 0.02 < dγ < 0.125          # ηk² = 0.05 within a factor ~2.5
end

@testset "MHD: vortex ring propagates" begin
    box = Box(32, 2.0)
    sim = MHDSim(box; cs = 1.0, eta = 0.0, sponge_width = 4)
    add_vortex_ring!(sim; R = 0.8, a = 0.25, z0 = -0.8, P0 = 0.3)
    function zcent()
        ωx, ωy, ωz = curl_central(sim.S[MMX], sim.S[MMY], sim.S[MMZ],
                                  box, sim.ip, sim.im)
        num = 0.0; den = 0.0
        for k in 1:box.n, j in 1:box.n, i in 1:box.n
            w = ωx[i, j, k]^2 + ωy[i, j, k]^2 + ωz[i, j, k]^2
            num += w * center(box, k)
            den += w
        end
        num / den
    end
    z0 = zcent()
    while sim.t < 4.0
        mhd_step!(sim)
    end
    @test zcent() > z0 + 0.1
end

@testset "MHD: flux ring is divergence-free and GLM keeps it so" begin
    box = Box(32, 2.0)
    sim = MHDSim(box; cs = 0.5, eta = 1e-3, sponge_width = 4)
    add_flux_ring!(sim; R = 0.8, a = 0.25, A0 = 0.5, Bt0 = 0.3)
    ip, im = sim.ip, sim.im
    function divb(S)
        n = box.n; h = 1 / (2box.dx); mx = 0.0
        for k in 1:n, j in 1:n, i in 1:n
            d = (S[MBX][ip[i], j, k] - S[MBX][im[i], j, k]) * h +
                (S[MBY][i, ip[j], k] - S[MBY][i, im[j], k]) * h +
                (S[MBZ][i, j, ip[k]] - S[MBZ][i, j, im[k]]) * h
            mx = max(mx, abs(d))
        end
        mx
    end
    @test divb(sim.S) < 1e-10        # curl-of-potential construction
    bref = maximum(abs, sim.S[MBX]) / box.dx
    while sim.t < 3.0
        mhd_step!(sim)
    end
    @test divb(sim.S) < 0.05 * bref  # GLM keeps residual small vs |B|/dx
    @test mhd_magnetic_energy(sim) > 0
end
