# v3 groundwork, experiment (a) of README §8: the anapole null, quantified.
# How completely do the fractal coil's external *fields* cancel while its
# external *potentials* survive?
#
# Observables, all evaluated pointwise from the line current (no grid):
#   1. mean |B| over direction-averaged shells of radius r  (synthetic
#      magnetometer — the field-side far signature)
#   2. mean |A| over the same shells (Coulomb gauge; gauge-dependent, shown
#      for orientation only)
#   3. Φ_link(ρ) = ∮ A·dl around a circle of radius ρ linking the torus
#      tube (centered on the tube at (R0,0,0), in the xz-plane). This is
#      gauge-invariant (= linked flux) and is exactly what an
#      Aharonov–Bohm-sensitive matter-wave probe reads. For an ideal
#      toroidal solenoid it stays at the tube flux for every ρ while B on
#      the loop → 0: potentials measurable where fields vanish.
#
# Variants: k0 (plain ring), k1 (16-turn coil), k2 (16×8 coil of coils),
# k1r/k2r (same + counter-oriented plain ring superposed, canceling the
# net loop current → the far dipole, i.e. the proper "anapole coil").
#
#   julia -t auto --project=. scripts/v3_anapole_null.jl
#
# Writes out/v3/anapole_null_{shells,link,moments}.csv.

using FractalToroid
using Printf

const OUT = joinpath(@__DIR__, "..", "out", "v3")
mkpath(OUT)

const EPS2 = 1e-6            # wire-core regularization (off-wire sampling)

"B(x) of the closed polyline current, direct segment sum (μ0 = I = 1)."
function bs_at(mid, dl, x, y, z)
    bx = by = bz = 0.0
    @inbounds for s in 1:size(mid, 2)
        rx = x - mid[1, s]; ry = y - mid[2, s]; rz = z - mid[3, s]
        r2 = rx * rx + ry * ry + rz * rz + EPS2
        iv = 1 / (r2 * sqrt(r2))
        bx += (dl[2, s] * rz - dl[3, s] * ry) * iv
        by += (dl[3, s] * rx - dl[1, s] * rz) * iv
        bz += (dl[1, s] * ry - dl[2, s] * rx) * iv
    end
    (bx / 4π, by / 4π, bz / 4π)
end

"A(x) of the closed polyline current, Coulomb gauge."
function A_at(mid, dl, x, y, z)
    ax = ay = az = 0.0
    @inbounds for s in 1:size(mid, 2)
        rx = x - mid[1, s]; ry = y - mid[2, s]; rz = z - mid[3, s]
        iv = 1 / sqrt(rx * rx + ry * ry + rz * rz + EPS2)
        ax += dl[1, s] * iv
        ay += dl[2, s] * iv
        az += dl[3, s] * iv
    end
    (ax / 4π, ay / 4π, az / 4π)
end

"~uniform unit directions (Fibonacci sphere)."
function fib_dirs(n)
    ga = π * (3 - sqrt(5.0))
    [begin
         zc = 1 - 2 * (i - 0.5) / n
         rc = sqrt(max(0.0, 1 - zc^2))
         (rc * cos(ga * i), rc * sin(ga * i), zc)
     end for i in 1:n]
end

"Reverse a curve's orientation: current sense flips, geometry unchanged."
revsegs(mid, dl) = (mid, -dl)

function variants()
    c0 = segments(fractal_coil(0))
    c1 = segments(fractal_coil(1; windings = 16))
    c2 = segments(fractal_coil(2; windings = [16, 8]))
    r0 = revsegs(c0...)
    cat2(a, b) = (hcat(a[1], b[1]), hcat(a[2], b[2]))
    ["k0" => c0, "k1" => c1, "k2" => c2,
     "k1r" => cat2(c1, r0), "k2r" => cat2(c2, r0)]
end

const DIRS = fib_dirs(192)
const NLOOP = 512
const RADII = [1.5 * 1.12^i for i in 0:24]          # shells 1.5 → ~17 R0
# AB loops: circles in the xz-plane through the hole center (0,0,0),
# centered at (s,0,0) with radius s — they thread the hole and link the
# tube for every s > 0.625, growing arbitrarily large. (A loop can't link
# and be entirely far away — that's the topology of the observable — but
# its far half recedes to distance 2s while Φ_link must persist.)
const SVALS = [0.7 * 1.14^i for i in 0:24]          # loop radius 0.7 → ~16

shellrows = ["variant,r,meanB,meanA"]
linkrows = ["variant,kind,s,phi_link,meanB_loop,meanB_farhalf,loop_len"]
momrows = ["variant,mz,Tz,absm,absT"]

for (name, (mid, dl)) in variants()
    m, T = current_moments(mid, dl)
    push!(momrows, @sprintf("%s,%.6g,%.6g,%.6g,%.6g", name, m[3], T[3],
                            sqrt(sum(abs2, m)), sqrt(sum(abs2, T))))
    for r in RADII
        sB = 0.0; sA = 0.0
        for d in DIRS
            b = bs_at(mid, dl, r * d[1], r * d[2], r * d[3])
            a = A_at(mid, dl, r * d[1], r * d[2], r * d[3])
            sB += sqrt(sum(abs2, b))
            sA += sqrt(sum(abs2, a))
        end
        push!(shellrows, @sprintf("%s,%.4f,%.6g,%.6g", name, r,
                                  sB / length(DIRS), sA / length(DIRS)))
    end
    for s in SVALS, ctrl in (false, true)
        # linking loop: center (s,0,0) in the xz-plane, radius s (through
        # the hole). Control: same circle shifted to the y = 3 plane —
        # same size and comparable distances, linking number 0.
        y = ctrl ? 3.0 : 0.0
        Φ = 0.0; sB = 0.0; sBfar = 0.0; nfar = 0
        for q in 1:NLOOP
            θ = 2π * (q - 0.5) / NLOOP
            x = s + s * cos(θ); z = s * sin(θ)
            tx = -sin(θ); tz = cos(θ)           # unit tangent
            dlq = 2π * s / NLOOP
            a = A_at(mid, dl, x, y, z)
            Φ += (a[1] * tx + a[3] * tz) * dlq
            b = bs_at(mid, dl, x, y, z)
            bb = sqrt(sum(abs2, b))
            sB += bb
            x > s && (sBfar += bb; nfar += 1)   # far half of the loop
        end
        push!(linkrows, @sprintf("%s,%s,%.4f,%.6g,%.6g,%.6g,%.6g",
                                 name, ctrl ? "ctrl" : "link", s, Φ,
                                 sB / NLOOP, sBfar / max(nfar, 1), 2π * s))
    end
    println("$name done ($(size(mid, 2)) segments)")
end

write(joinpath(OUT, "anapole_null_shells.csv"), join(shellrows, "\n") * "\n")
write(joinpath(OUT, "anapole_null_link.csv"), join(linkrows, "\n") * "\n")
write(joinpath(OUT, "anapole_null_moments.csv"), join(momrows, "\n") * "\n")
println("wrote out/v3/anapole_null_{shells,link,moments}.csv")
