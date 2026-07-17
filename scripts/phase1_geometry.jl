# Phase 1 — fractal coil geometry + Biot–Savart magnetostatics.
# For nesting depths k = 0..K: build the coil, render it, compute the static
# B field, the magnetic dipole and toroidal (anapole) moments, magnetic
# helicity, and the far-field radial profile.
#
#   julia -t auto --project=. scripts/phase1_geometry.jl [K] [ngrid]
#
# Outputs land in out/phase1/.

using FractalToroid
using Printf

const K = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
const NGRID = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 64
const R0 = 1.0
const RATIO = 0.25          # tube radius shrink per level (Clarage: 1/4)
const WINDINGS = [16, 8, 6] # turns per parent turn at each level
const HALF = 2.0            # grid spans [-HALF, HALF]³
const OUT = joinpath(@__DIR__, "..", "out", "phase1")

mkpath(OUT)

"Project curve points into a 2D density image."
function curve_image(P::Matrix{Float64}, axes::Tuple{Int,Int}; res = 700,
                     lim = 1.5)
    img = zeros(res, res)
    for s in 1:size(P, 2)
        i = round(Int, (P[axes[1], s] + lim) / 2lim * (res - 1)) + 1
        j = round(Int, (P[axes[2], s] + lim) / 2lim * (res - 1)) + 1
        (1 <= i <= res && 1 <= j <= res) && (img[i, j] += 1)
    end
    log10.(img .+ 1)
end

box = Box(NGRID, HALF)
println("Phase 1: fractal coil magnetostatics on a $(NGRID)³ grid, ",
        "depths k = 0..$K, $(Threads.nthreads()) threads")
println()
@printf("%2s %9s %11s %11s %11s %11s %9s\n",
        "k", "segments", "m_z", "|T|", "|T|/|m|", "helicity", "far-slope")

summary = ["k,segments,mz,Tnorm,helicity,far_slope"]
for k in 0:K
    coil = fractal_coil(k; R0, ratio = RATIO, windings = WINDINGS[1:k],
                        ppt = 24)
    mid, dl = segments(coil)

    heatmap_png(joinpath(OUT, "coil_k$(k)_xy.png"),
                curve_image(coil.P, (1, 2)))
    heatmap_png(joinpath(OUT, "coil_k$(k)_xz.png"),
                curve_image(coil.P, (1, 3)))

    m, T = current_moments(mid, dl)
    Bx, By, Bz = biot_savart(mid, dl, box)
    Ax, Ay, Az = vector_potential(mid, dl, box)
    H = helicity(Ax, Ay, Az, Bx, By, Bz, box)

    # |B| slice in the xz-plane through y ≈ 0
    j0 = box.n ÷ 2
    slice = [sqrt(Bx[i, j0, kk]^2 + By[i, j0, kk]^2 + Bz[i, j0, kk]^2)
             for i in 1:box.n, kk in 1:box.n]
    heatmap_png(joinpath(OUT, "Bmag_k$(k)_xz.png"), slice;
                logscale = true, upscale = max(1, 512 ÷ box.n))

    r, prof = shell_profile(Bx, By, Bz, box; nbins = 40)
    slope = powerlaw_slope(r, prof; rmin = 1.4, rmax = 0.95 * HALF)
    open(joinpath(OUT, "profile_k$(k).csv"), "w") do io
        println(io, "r,mean_absB")
        for b in eachindex(r)
            println(io, r[b], ",", prof[b])
        end
    end

    mn = sqrt(sum(abs2, m))
    Tn = sqrt(sum(abs2, T))
    @printf("%2d %9d %11.4g %11.4g %11.4g %11.4g %9.2f\n",
            k, size(mid, 2), m[3], Tn, Tn / mn, H, slope)
    push!(summary, join((k, size(mid, 2), m[3], Tn, H, slope), ","))
end

write(joinpath(OUT, "summary.csv"), join(summary, "\n") * "\n")
println("\nWrote images and CSVs to $(abspath(OUT))")
