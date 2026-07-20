# Validation of the v3 GPE + potentials module (src/gpe.jl).

@testset "GPE free packet dispersion" begin
    box = Box(32, 8.0)                    # dx = 0.5
    sim = GPESim(box)                     # q = 0, g = 0: free Schrödinger
    gpe_packet!(sim; sigma = 1.5)
    n0 = gpe_norm(sim)
    while sim.t < 2.0
        gpe_step!(sim)
    end
    @test abs(gpe_norm(sim) - n0) / n0 < 1e-6
    # |ψ|² width: σ(t)² = σ0²(1 + (t/2mσ0²)²), per axis
    ρtot = 0.0; s2 = 0.0
    for k in 1:box.n, j in 1:box.n, i in 1:box.n
        ρ = sim.S[GR][i, j, k]^2 + sim.S[GI][i, j, k]^2
        ρtot += ρ
        s2 += ρ * center(box, i)^2        # packet is centered at 0
    end
    σmeas = sqrt(s2 / ρtot)
    σth = 1.5 * sqrt(1 + (sim.t / (2 * 1.5^2))^2)
    @test abs(σmeas - σth) / σth < 0.03
end

@testset "GPE uniform chemical potential" begin
    # ψ = √n0 e^{−iμt} with μ = g n0
    box = Box(16, 4.0)
    sim = GPESim(box; g = 0.5)
    gpe_uniform!(sim; n0 = 1.0)
    while sim.t < 1.0
        gpe_step!(sim)
    end
    ph = atan(sim.S[GI][8, 8, 8], sim.S[GR][8, 8, 8])
    @test abs(ph - (-0.5 * sim.t)) < 1e-3
    # RK3's |R(−iμdt)| = 1 − (μdt)⁴/24 damps the norm ~1e-6 over 16 steps
    @test abs(gpe_norm(sim) - (2 * 4.0)^3) / (2 * 4.0)^3 < 1e-5
end

@testset "GPE plasma oscillation" begin
    # charged uniform condensate + uniform A kick: ∂ttA = −(q²n/m)A,
    # ω_p = 1 — the coherent-matter analogue of the v1 Langmuir test
    box = Box(16, 4.0)
    sim = GPESim(box; q = 1.0, m = 1.0)
    gpe_uniform!(sim; n0 = 1.0)
    A0 = 1e-3
    sim.S[GAX] .= A0
    tz = 0.0                              # first zero crossing of ⟨Ax⟩
    prev = A0
    while sim.t < 3.5
        gpe_step!(sim)
        cur = sum(sim.S[GAX]) / box.n^3
        if tz == 0.0 && prev > 0 && cur <= 0
            # interpolate the crossing inside the step
            tz = sim.t - sim.dt * cur / (cur - prev)
        end
        prev = cur
    end
    @test abs(tz - π / 2) / (π / 2) < 0.03
    @test abs(prev) < A0                  # bounded oscillation, no growth
end

@testset "GPE vortex pair winding" begin
    box = Box(32, 8.0)
    sim = GPESim(box; g = 0.25)           # ξ = 1/√(2gn0) = √2 ≈ 2.8 dx
    gpe_vortex_pair!(sim; n0 = 1.0, d = 8.0, xi = sqrt(2.0))
    @test gpe_winding(sim, 4.0, 0.0; radius = 2.0) == 1
    @test gpe_winding(sim, -4.0, 0.0; radius = 2.0) == -1
    n0 = gpe_norm(sim)
    while sim.t < 2.0
        gpe_step!(sim)
    end
    # quantized circulation is topological: it must survive evolution
    @test gpe_winding(sim, 4.0, 0.0; radius = 2.0) == 1
    @test gpe_winding(sim, -4.0, 0.0; radius = 2.0) == -1
    # cores are ~3 cells wide: RK3 grid-scale dissipation costs ~1e-3 of
    # the norm over t = 2 (the winding above is exact — it's topological)
    @test abs(gpe_norm(sim) - n0) / n0 < 5e-3
    e = gpe_energy(sim)
    @test isfinite(e.kinetic) && isfinite(e.interaction)
end
