# v2 path A — the v1 fractal coil evolved under fluid MHD (README §7.4).
# B from Biot–Savart of the coil, density blob following the tube, v = 0.
# Decay channels: reconnection + resistivity. Key output: τ(k) under MHD.
#
#   julia -t auto --project=. scripts/v2_pathA_fractal_mhd.jl [K] [ngrid] [t_end]
#
# Outputs: out/v2/pathA_N<ngrid>/

using FractalToroid
using Printf

const K = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
const NGRID = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 64
const T_END = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 20.0
const HALF = 2.0
const R0 = 1.0
const RATIO = 0.25
const WINDINGS = [16, 8, 6]
const B0 = 1.0              # peak |B| after rescaling
const RHO_AMB = 0.05
const ETA = 2e-3
const CS = 0.55
const SPONGE = 6
const OUT = joinpath(@__DIR__, "..", "out", "v2", "pathA_N$(NGRID)")

mkpath(OUT)
box = Box(NGRID, HALF)

println("v2 path A: fractal coil under MHD, $(NGRID)³, t_end=$T_END, ",
        "S=$(1/ETA), $(Threads.nthreads()) threads")
@printf("%2s %12s %12s %12s %12s\n", "k", "E_mag(0)", "E_mag(end)", "tau",
        "|T|(end)/|T|(0)")

summary = ["k,Emag0,EmagEnd,tau,Tratio"]
for k in 0:K
    coil = fractal_coil(k; R0, ratio = RATIO, windings = WINDINGS[1:k],
                        ppt = 24)
    mid, dl = segments(coil)
    sim = MHDSim(box; cs = CS, eta = ETA, sponge_width = SPONGE,
                 rho_floor = 1e-3)

    Bx, By, Bz = biot_savart(mid, dl, box; eps = 2 * box.dx)
    bmax = maximum(sqrt.(Bx .^ 2 .+ By .^ 2 .+ Bz .^ 2))
    s = B0 / bmax
    sim.S[MBX] .= s .* Bx
    sim.S[MBY] .= s .* By
    sim.S[MBZ] .= s .* Bz
    # density blob follows |B| (matter where the field is), ambient floor
    for kk in 1:box.n, j in 1:box.n, i in 1:box.n
        bb = s * sqrt(Bx[i, j, kk]^2 + By[i, j, kk]^2 + Bz[i, j, kk]^2)
        sim.S[MRHO][i, j, kk] = RHO_AMB + bb / B0
    end

    m0, T0 = grid_moments(sim)
    e0 = mhd_magnetic_energy(sim)
    ts = Float64[]; es = Float64[]
    rows = ["t,E_kin,E_mag,Tnorm"]
    tlog = 0.0
    while sim.t < T_END
        mhd_step!(sim)
        push!(ts, sim.t); push!(es, mhd_magnetic_energy(sim))
        if sim.t >= tlog
            m, T = grid_moments(sim)
            push!(rows, join((round(sim.t, digits = 3),
                              mhd_kinetic_energy(sim), es[end],
                              sqrt(sum(abs2, T))), ","))
            tlog += 0.25
        end
    end
    m1, T1 = grid_moments(sim)
    τ = efold_time(ts, es)
    Tr = sqrt(sum(abs2, T1)) / max(sqrt(sum(abs2, T0)), 1e-30)
    write(joinpath(OUT, "timeseries_k$(k).csv"), join(rows, "\n") * "\n")
    @printf("%2d %12.4g %12.4g %12.4g %12.4g\n", k, e0, es[end], τ, Tr)
    push!(summary, join((k, e0, es[end], τ, Tr), ","))
end
write(joinpath(OUT, "summary.csv"), join(summary, "\n") * "\n")
println("Wrote $(abspath(OUT))")
