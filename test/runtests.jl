using Test
using FractalToroid

@testset "geometry" begin
    # every level closes: last point connects near first
    for k in 0:2
        c = fractal_coil(k; windings = 6, ppt = 32)
        P = c.P
        gap = sqrt(sum(abs2, P[:, 1] .- P[:, end]))
        seg = sqrt(sum(abs2, P[:, 2] .- P[:, 1]))
        @test gap < 3 * seg
    end
    # level-1 coil stays within tube radius a1 of the base circle
    c = fractal_coil(1; R0 = 1.0, ratio = 0.25, windings = 8, ppt = 64)
    for i in 1:npoints(c)
        ρ = sqrt(c.P[1, i]^2 + c.P[2, i]^2)
        d = sqrt((ρ - 1.0)^2 + c.P[3, i]^2)   # distance to base circle
        @test d < 0.25 + 1e-6
    end
    # winding number: z changes sign 2*w times for w poloidal turns
    # (phase offset so no sample lands exactly on z = 0)
    c = fractal_coil(1; R0 = 1.0, ratio = 0.25, windings = 8, ppt = 64,
                     phase = 0.3)
    zs = c.P[3, :]
    crossings = count(i -> zs[i] * zs[i%length(zs)+1] < 0, 1:length(zs))
    @test crossings == 2 * 8
end

@testset "biot-savart loop" begin
    # circular loop R=1, I=1: B(0) = 1/2 (μ0=1); on-axis B(z) = R²/2(R²+z²)^{3/2}
    c = fractal_coil(0; R0 = 1.0, minpts = 2048)
    mid, dl = segments(c)
    box = Box(33, 1.6)   # odd n so a cell center sits at the origin
    Bx, By, Bz = biot_savart(mid, dl, box; eps = 1e-9)
    i0 = 17
    @test abs(Bz[i0, i0, i0] - 0.5) < 2e-3
    @test abs(Bx[i0, i0, i0]) < 1e-6
    z = center(box, 25)
    @test abs(Bz[i0, i0, 25] - 0.5 / (1 + z^2)^1.5) < 2e-3
end

@testset "moments" begin
    c = fractal_coil(0; R0 = 1.0, minpts = 2048)
    mid, dl = segments(c)
    m, T = current_moments(mid, dl)
    @test abs(m[3] - π) < 1e-3           # m = I π R²
    @test abs(m[1]) + abs(m[2]) < 1e-9
    @test maximum(abs, T) < 1e-9          # planar centered loop: T = 0
    # level-1 toroidal coil: dipole stays ≈ that of the net one-turn loop,
    # and a finite anapole moment along z appears
    c1 = fractal_coil(1; R0 = 1.0, ratio = 0.25, windings = 12, ppt = 64)
    mid1, dl1 = segments(c1)
    m1, T1 = current_moments(mid1, dl1)
    @test abs(m1[3] - π) / π < 0.15
    @test abs(T1[3]) > 10 * max(abs(T1[1]), abs(T1[2]))
    @test abs(T1[3]) > 0.01
end

@testset "vacuum FDTD" begin
    box = Box(32, 1.0)
    ip, im = _wrap(box.n)
    S = zero_state(box, 6)
    # zero-mean Ez pulse (the k=0 mode of E does not propagate and would
    # never reach the sponge). Energy must be conserved (periodic, no
    # sponge) and div B must stay at machine zero.
    for k in 1:box.n, j in 1:box.n, i in 1:box.n
        x, y = center(box, i), center(box, j)
        S[EZ][i, j, k] = x * exp(-(x^2 + y^2) / 0.1)
    end
    e0 = field_energy(S, box)
    S1 = zero_state(box, 6); S2 = zero_state(box, 6); K = zero_state(box, 6)
    dt = 0.3 * box.dx
    rhs! = (Kv, Sv, t) -> em_rhs!(Kv, Sv, box, ip, im)
    for s in 1:200
        ssprk3!(S, S1, S2, K, rhs!, dt, (s - 1) * dt)
    end
    e1 = field_energy(S, box)
    # SSP-RK3 is weakly dissipative at grid-scale modes (unlike leapfrog);
    # well-resolved fields lose ≪1% over hundreds of steps.
    @test abs(e1 - e0) / e0 < 1e-2
    @test maximum(abs, div_B(S, box, ip)) < 1e-10
    # with a sponge, the pulse energy must be absorbed almost entirely
    mask = make_sponge(box, 6, dt)
    absorbed = 0.0
    for s in 1:400
        ssprk3!(S, S1, S2, K, rhs!, dt, 0.0)
        absorbed += apply_sponge!(S, mask, box)
    end
    @test field_energy(S, box) / e0 < 0.02
    @test absorbed / e0 > 0.9
end

@testset "Langmuir oscillation" begin
    # uniform quasineutral plasma, sinusoidal velocity kick along x:
    # cold plasma oscillates at ω_p = 1 → kinetic energy at ω = 2.
    box = Box(24, π)     # domain [−π, π)
    sim = FluidSim(box; dt = 0.25 * box.dx)
    u0 = 1e-3
    for k in 1:box.n, j in 1:box.n, i in 1:box.n
        sim.S[FPX][i, j, k] = u0 * sin(node(box, i))
    end
    ke0 = kinetic_energy(sim)
    kes = Float64[]
    ts = Float64[]
    while sim.t < 2.2π
        step!(sim)
        push!(kes, kinetic_energy(sim))
        push!(ts, sim.t)
    end
    # first minimum of KE should sit at t = π/2 (quarter period of ω_p)
    imin = argmin(abs.(ts .- π / 2))
    @test kes[imin] / ke0 < 0.05
    # KE back near max at t = π
    imax = argmin(abs.(ts .- π))
    @test kes[imax] / ke0 > 0.8
    # Gauss law maintained
    @test gauss_residual(sim) < 1e-6
end

@testset "png writer" begin
    path = joinpath(mktempdir(), "t.png")
    heatmap_png(path, rand(20, 15))
    bytes = read(path)
    @test bytes[1:8] == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
    @test length(bytes) > 100
end
include("mhd_tests.jl")
include("gpe_tests.jl")
