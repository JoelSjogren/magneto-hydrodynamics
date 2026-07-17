# Phase 2 — vacuum FDTD with a prescribed AC current through the fractal
# coil, then switch-off. Two measurements per nesting depth k:
#   1. steady-state radiated power during the AC drive (does deeper nesting
#      radiate less, as the anapole "self-containment" claim implies?)
#   2. field energy left in the core after the current is switched off
#      (vacuum Maxwell: it must vanish on the light-crossing time).
#
#   julia -t auto --project=. scripts/phase2_fdtd.jl [K] [ngrid]
#
# Outputs land in out/phase2/.

using FractalToroid
using Printf

const K = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
const NGRID = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 64
const R0 = 1.0
const RATIO = 0.25
const WINDINGS = [16, 8, 6]
const HALF = 2.0
const SPONGE = 8            # sponge width, cells
const OMEGA = 2π / 1.0      # drive frequency (vacuum wavelength = 1)
const T_RAMP = 2.0
const T_OFF = 10.0          # drive switch-off time
const T_END = 16.0
const OUT = joinpath(@__DIR__, "..", "out", "phase2")

mkpath(OUT)

"Smooth ramp 0→1 over [0, w]."
smoothstep(x) = x <= 0 ? 0.0 : x >= 1 ? 1.0 : x * x * (3 - 2x)

"Drive waveform: ramped sinusoid, smoothly off after T_OFF."
g(t) = sin(OMEGA * t) * smoothstep(t / T_RAMP) *
       (1 - smoothstep((t - T_OFF) / T_RAMP))

"Field energy outside the sponge (the physical core region)."
function core_energy(S, box, margin)
    e = 0.0
    n = box.n
    for f in EX:BZ
        A = S[f]
        @inbounds for k in margin+1:n-margin, j in margin+1:n-margin,
                      i in margin+1:n-margin
            e += A[i, j, k]^2
        end
    end
    0.5 * e * cellvol(box)
end

box = Box(NGRID, HALF)
ip, im = _wrap(box.n)
dt = 0.3 * box.dx
mask = make_sponge(box, SPONGE, dt)
margin = SPONGE + 2
nsteps = ceil(Int, T_END / dt)

println("Phase 2: driven coil in vacuum, $(NGRID)³ grid, dt=$(round(dt, sigdigits=3)), ",
        "$nsteps steps, k = 0..$K, $(Threads.nthreads()) threads")
println()
@printf("%2s %12s %14s %14s\n",
        "k", "P_rad (AC)", "E_core(t_off)", "E_core(t_end)")

summary = ["k,P_rad_steady,E_core_at_off,E_core_end"]
for k in 0:K
    coil = fractal_coil(k; R0, ratio = RATIO, windings = WINDINGS[1:k],
                        ppt = 24)
    mid, dl = segments(coil)
    Jx = zeros(box.n, box.n, box.n)
    Jy = zeros(box.n, box.n, box.n)
    Jz = zeros(box.n, box.n, box.n)
    splat_current!(Jx, Jy, Jz, mid, dl, box; sigma = 1.5 * box.dx)

    S = zero_state(box, 6)
    S1 = zero_state(box, 6); S2 = zero_state(box, 6); Kt = zero_state(box, 6)
    rhs! = (Kv, Sv, t) -> em_rhs!(Kv, Sv, box, ip, im;
                                  Jx, Jy, Jz, gJ = g(t))

    absorbed = 0.0
    prev_absorbed = 0.0
    prad_window = Float64[]
    e_at_off = 0.0
    rows = ["t,E_core,absorbed_cum"]
    t = 0.0
    for s in 1:nsteps
        ssprk3!(S, S1, S2, Kt, rhs!, dt, t)
        absorbed += apply_sponge!(S, mask, box)
        t += dt
        # steady-state radiated power: averaged over the last 2 drive
        # periods before switch-off
        if T_OFF - 4π / OMEGA < t <= T_OFF
            push!(prad_window, (absorbed - prev_absorbed) / dt)
        end
        prev_absorbed = absorbed
        if abs(t - T_OFF) < dt / 2
            e_at_off = core_energy(S, box, margin)
        end
        if s % 10 == 0
            push!(rows, join((round(t, digits = 3),
                              core_energy(S, box, margin), absorbed), ","))
        end
    end
    e_end = core_energy(S, box, margin)
    prad = sum(prad_window) / max(1, length(prad_window))

    j0 = box.n ÷ 2
    slice = [sqrt(S[BX][i, j0, kk]^2 + S[BY][i, j0, kk]^2 + S[BZ][i, j0, kk]^2)
             for i in 1:box.n, kk in 1:box.n]
    heatmap_png(joinpath(OUT, "Bmag_end_k$(k).png"), slice;
                logscale = true, upscale = max(1, 512 ÷ box.n))
    write(joinpath(OUT, "timeseries_k$(k).csv"), join(rows, "\n") * "\n")

    @printf("%2d %12.4g %14.4g %14.4g\n", k, prad, e_at_off, e_end)
    push!(summary, join((k, prad, e_at_off, e_end), ","))
end

write(joinpath(OUT, "summary.csv"), join(summary, "\n") * "\n")
println("\nWrote time series and images to $(abspath(OUT))")
