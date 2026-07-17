# Phase 3 — self-consistent cold Euler–Maxwell: seed the electron fluid with
# the fractal-coil velocity pattern (quasineutral, f = 1), initialize B from
# Biot–Savart of that current, then release everything and watch the
# relaxation. Key output: decay time τ of the core energy vs nesting depth k.
#
#   julia -t auto --project=. scripts/phase3_euler_maxwell.jl [K] [ngrid] [t_end] [mode]
#
# mode = uniform : coil embedded in a uniform quasineutral plasma (n = 1
#                  everywhere). Radiation below ω_p is trapped by the medium.
# mode = ball    : plasma density follows the coil tube (n = n_b = floor +
#                  normalized tube profile) — an isolated object in
#                  near-vacuum, the configuration the persistence claim is
#                  actually about.
#
# Outputs land in out/phase3/ (files suffixed with the mode).

using FractalToroid
using Printf

const K = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
const NGRID = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 48
const T_END = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 60.0
const MODE = length(ARGS) >= 4 ? ARGS[4] : "uniform"
const N_OUT = 1e-3          # density floor outside the tube (ball mode)
const R0 = 1.0
const RATIO = 0.25
const WINDINGS = [16, 8, 6]
const HALF = 2.0
const SPONGE = 6
const U0 = 0.05             # peak fluid speed, units of c
const FRAME_DT = 0.25       # render a video frame every Δt (0 = off)
const UP = 4                # frame upscale factor (pixels per cell)
const GAP = 8               # gap between the |B| and n panels, px
const NSCALE = 2.5          # density colormap saturates at this value
const OUT = joinpath(@__DIR__, "..", "out", "phase3")

mkpath(OUT)

function core_energy(sim::FluidSim, margin)
    box = sim.box
    n = box.n
    e = 0.0
    for f in EX:BZ
        A = sim.S[f]
        @inbounds for k in margin+1:n-margin, j in margin+1:n-margin,
                      i in margin+1:n-margin
            e += A[i, j, k]^2
        end
    end
    N, px, py, pz = sim.S[FN], sim.S[FPX], sim.S[FPY], sim.S[FPZ]
    @inbounds for k in margin+1:n-margin, j in margin+1:n-margin,
                  i in margin+1:n-margin
        nn = max(N[i, j, k], sim.n_floor)
        e += (px[i, j, k]^2 + py[i, j, k]^2 + pz[i, j, k]^2) / nn
    end
    0.5 * e * cellvol(box)
end

box = Box(NGRID, HALF)
margin = SPONGE + 2
println("Phase 3: Euler–Maxwell relaxation, $(NGRID)³ grid, t_end=$T_END, ",
        "u0=$U0 c, k = 0..$K, $(Threads.nthreads()) threads")
println()
@printf("%2s %12s %12s %12s %10s\n",
        "k", "E_core(0)", "E_core(end)", "tau", "gauss_res")

summary = ["k,E0,Eend,tau,gauss_residual"]
for k in 0:K
    coil = fractal_coil(k; R0, ratio = RATIO, windings = WINDINGS[1:k],
                        ppt = 24)
    mid, dl = segments(coil)

    sim = FluidSim(box; sponge_width = SPONGE)

    # velocity seed: Gaussian-tube current pattern, scaled to peak speed U0
    jx = zeros(box.n, box.n, box.n)
    jy = zeros(box.n, box.n, box.n)
    jz = zeros(box.n, box.n, box.n)
    sigma = 1.5 * box.dx
    splat_current!(jx, jy, jz, mid, dl, box; sigma, stagger = :node)
    jmag = sqrt.(jx .^ 2 .+ jy .^ 2 .+ jz .^ 2)
    jmax = maximum(jmag)
    s = U0 / jmax
    @. sim.S[FPX] = -s * jx           # p = n·u = −s·j  ⇒ J = −p = s·j
    @. sim.S[FPY] = -s * jy
    @. sim.S[FPZ] = -s * jz
    if MODE == "ball"
        # electron density (and matching ion background) confined to the
        # tube: quasineutral everywhere, near-vacuum outside. p is kept
        # proportional to j, so J and the Biot–Savart B stay consistent.
        @. sim.S[FN] = N_OUT + jmag / jmax
        sim.nb .= sim.S[FN]
    end

    # B from Biot–Savart of the same current (∇×B = J_fluid), E = 0:
    # a consistent magnetostatic start for the quasineutral f = 1 state.
    Bx, By, Bz = biot_savart_yee(mid, dl, box; I = s, eps = sigma)
    sim.S[BX] .= Bx
    sim.S[BY] .= By
    sim.S[BZ] .= Bz

    e0 = core_energy(sim, margin)
    rows = ["t,E_core,E_kin,E_field,absorbed_cum,gauss_res"]
    ts = Float64[]
    es = Float64[]
    nsnap = 0

    # video frames rendered inline (no intermediate data files): PNG with
    # two panels, left = log10 |B| (fixed window set by the initial slice
    # maximum), right = electron density n (fixed scale 0..NSCALE).
    # Slice: xz-plane through y ≈ 0.
    framedir = joinpath(OUT, "frames", "k$(k)_$(MODE)")
    frame_every = FRAME_DT > 0 ? max(1, round(Int, FRAME_DT / sim.dt)) : 0
    nframe = 0
    j0 = box.n ÷ 2
    logBhi = log10(max(1e-300,
        maximum(sqrt(sim.S[BX][i, j0, kk]^2 + sim.S[BY][i, j0, kk]^2 +
                     sim.S[BZ][i, j0, kk]^2) for i in 1:box.n, kk in 1:box.n)))
    npx = box.n * UP
    function render_frame()
        frame_every > 0 || return
        rgb = fill(0x12, 3, 2npx + GAP, npx)
        for kk in 1:box.n, ii in 1:box.n
            b = sqrt(sim.S[BX][ii, j0, kk]^2 + sim.S[BY][ii, j0, kk]^2 +
                     sim.S[BZ][ii, j0, kk]^2)
            vB = clamp((log10(max(b, 1e-300)) - (logBhi - 3.5)) / 3.5, 0.0, 1.0)
            vN = clamp(sim.S[FN][ii, j0, kk] / NSCALE, 0.0, 1.0)
            cB = FractalToroid._colormap(vB)
            cN = FractalToroid._colormap(vN)
            for dj in 1:UP, di in 1:UP
                row = (box.n - kk) * UP + dj
                col = (ii - 1) * UP + di
                for ch in 1:3
                    rgb[ch, col, row] = cB[ch]
                    rgb[ch, npx + GAP + col, row] = cN[ch]
                end
            end
        end
        save_png(joinpath(framedir, "frame_$(lpad(nframe, 5, '0')).png"), rgb)
        nframe += 1
    end
    if frame_every > 0
        rm(framedir; force = true, recursive = true)
        mkpath(framedir)
        render_frame()
    end

    nstep = 0
    while sim.t < T_END
        step!(sim)
        nstep += 1
        frame_every > 0 && nstep % frame_every == 0 && render_frame()
        if length(ts) % 5 == 0
            push!(rows, join((round(sim.t, digits = 3),
                              core_energy(sim, margin),
                              kinetic_energy(sim),
                              field_energy(sim.S, sim.box),
                              sim.absorbed,
                              round(gauss_residual(sim), sigdigits = 4)), ","))
        end
        push!(ts, sim.t)
        push!(es, core_energy(sim, margin))
        if sim.t >= nsnap * T_END / 4
            j0 = box.n ÷ 2
            slice = [sqrt(sim.S[BX][i, j0, kk]^2 + sim.S[BY][i, j0, kk]^2 +
                          sim.S[BZ][i, j0, kk]^2)
                     for i in 1:box.n, kk in 1:box.n]
            heatmap_png(joinpath(OUT, "Bmag_k$(k)_$(MODE)_t$(round(Int, sim.t)).png"),
                        slice; logscale = true,
                        upscale = max(1, 512 ÷ box.n))
            nsnap += 1
        end
    end

    τ = efold_time(ts, es)
    gres = gauss_residual(sim)
    write(joinpath(OUT, "timeseries_k$(k)_$(MODE).csv"), join(rows, "\n") * "\n")
    @printf("%2d %12.4g %12.4g %12.4g %10.3g\n", k, e0, es[end], τ, gres)
    push!(summary, join((k, e0, es[end], τ, gres), ","))
end

write(joinpath(OUT, "summary_$(MODE).csv"), join(summary, "\n") * "\n")
println("\nWrote time series and images to $(abspath(OUT))")
