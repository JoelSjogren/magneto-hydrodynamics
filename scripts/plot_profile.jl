# Component profile of a 256³ GPU frame cycle (scripts/bench_gpu_profile.jl):
# is the bottleneck physics (GPU stepping) or rendering/diagnostics/I/O?
# Dependency-free SVG, same convention as plot_anapole.jl. Categorical
# colors from the project's validated palette (dataviz skill).
#
#   julia --project=. scripts/plot_profile.jl
#
# Reads  out/v2/profile_N256.csv
# Writes out/plots/profile_pie.svg, out/plots/profile_bar.svg

const OUT = joinpath(@__DIR__, "..", "out", "plots")
mkpath(OUT)

const BLUE, GREEN, MAGENTA, YELLOW, AQUA, ORANGE, VIOLET, RED =
    "#2a78d6", "#008300", "#e87ba4", "#eda100",
    "#1baf7a", "#eb6834", "#4a3aa7", "#e34948"
const INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
const SURFACE = "#fcfcfb"

function read_kv(path)
    d = Dict{String,Float64}()
    for l in readlines(path)[2:end]
        k, v = split(l, ",")
        d[k] = parse(Float64, v)
    end
    d
end

d = read_kv(joinpath(@__DIR__, "..", "out", "v2", "profile_N256.csv"))

physics = d["step_stretch"]
diag = d["frame_products"] + d["frame_scalars"] + d["spectrum_cpu"] +
       d["render2d_cpu"]
io_ = d["png_panel"] + d["png_vol3d"] + d["voldump"]
total = physics + diag + io_

# ---------------------------------------------------------------- pie ----
function svg_pie(path; title, subtitle, slices, W = 760, H = 460)
    cx, cy, R = 210, H ÷ 2 + 6, 150
    tot = sum(v for (_, v, _) in slices)
    io = IOBuffer()
    print(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$H" font-family="sans-serif">
<rect width="$W" height="$H" fill="$SURFACE"/>
<text x="$(W ÷ 2)" y="26" text-anchor="middle" font-size="16" fill="$INK">$title</text>
<text x="$(W ÷ 2)" y="46" text-anchor="middle" font-size="12" fill="$MUTED">$subtitle</text>
""")
    a0 = -π / 2
    for (i, (label, v, color)) in enumerate(slices)
        frac = v / tot
        a1 = a0 + 2π * frac
        large = frac > 0.5 ? 1 : 0
        x0, y0 = cx + R * cos(a0), cy + R * sin(a0)
        x1, y1 = cx + R * cos(a1), cy + R * sin(a1)
        print(io, """<path d="M $cx,$cy L $x0,$y0 A $R,$R 0 $large,1 $x1,$y1 Z" fill="$color" stroke="$SURFACE" stroke-width="2"/>
""")
        # direct label for slices ≥ 3%; small slices get a leader line
        amid = (a0 + a1) / 2
        pct = round(100frac, digits = frac < 0.01 ? 2 : 1)
        if frac >= 0.03
            lx, ly = cx + 0.62R * cos(amid), cy + 0.62R * sin(amid)
            print(io, """<text x="$lx" y="$ly" text-anchor="middle" font-size="13" fill="white" font-weight="600">$(pct)%</text>
""")
        end
        a0 = a1
    end
    # legend
    ly0 = 90
    for (i, (label, v, color)) in enumerate(slices)
        frac = 100v / tot
        y = ly0 + 26i
        pct = frac < 0.1 ? round(frac, digits = 2) : round(frac, digits = 1)
        print(io, """<rect x="430" y="$(y - 12)" width="14" height="14" fill="$color"/>
<text x="452" y="$y" font-size="13" fill="$INK">$label</text>
<text x="$(W - 20)" y="$y" text-anchor="end" font-size="13" fill="$INK2">$(pct)%</text>
""")
    end
    print(io, "</svg>\n")
    write(path, take!(io))
end

svg_pie(joinpath(OUT, "profile_pie.svg");
    title = "Where does the frame-cycle time go? (256³, GTX 1080, FP64)",
    subtitle = "$(round(Int, total)) ms per frame (t += 0.25, ~$(round(Int, d["steps_per_frame"])) physics steps)",
    slices = [("GPU physics stepping", physics, BLUE),
              ("Frame diagnostics (render+moments)", diag, ORANGE),
              ("File I/O (PNG + volume)", io_, GREEN)])

# ---------------------------------------------------------------- bar ----
function svg_bar_log(path; title, subtitle, bars, W = 760, H = 400)
    ml, mr, mt, mb = 230, 100, 56, 40
    pw, ph = W - ml - mr, H - mt - mb
    vals = [v for (_, v, _) in bars]
    lo, hi = 0.05, 10^ceil(log10(maximum(vals)))
    X(v) = ml + (log10(max(v, lo)) - log10(lo)) / (log10(hi) - log10(lo)) * pw
    io = IOBuffer()
    print(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$H" font-family="sans-serif">
<rect width="$W" height="$H" fill="$SURFACE"/>
<text x="$(W ÷ 2)" y="24" text-anchor="middle" font-size="16" fill="$INK">$title</text>
<text x="$(W ÷ 2)" y="42" text-anchor="middle" font-size="12" fill="$MUTED">$subtitle</text>
""")
    e0, e1 = ceil(Int, log10(lo)), floor(Int, log10(hi))
    for e in e0:e1
        x = X(10.0^e)
        print(io, """<line x1="$x" y1="$mt" x2="$x" y2="$(mt + ph)" stroke="$GRID"/>
<text x="$x" y="$(mt + ph + 20)" text-anchor="middle" font-size="11" fill="$MUTED">10^$e</text>
""")
    end
    bh = ph / length(bars) * 0.6
    for (i, (label, v, color)) in enumerate(bars)
        y = mt + (i - 0.8) * ph / length(bars)
        x0 = X(lo)
        x1 = X(v)
        print(io, """<rect x="$x0" y="$y" width="$(x1 - x0)" height="$bh" rx="3" fill="$color"/>
<text x="$(ml - 10)" y="$(y + bh / 2 + 4)" text-anchor="end" font-size="12" fill="$INK">$label</text>
<text x="$(x1 + 8)" y="$(y + bh / 2 + 4)" font-size="12" fill="$INK2">$(round(v, sigdigits = 3)) ms</text>
""")
    end
    print(io, """<line x1="$ml" y1="$(mt + ph)" x2="$(ml + pw)" y2="$(mt + ph)" stroke="$MUTED"/>
<text x="$(ml + pw ÷ 2)" y="$(H - 8)" text-anchor="middle" font-size="12" fill="$INK2">milliseconds per frame (log scale)</text>
</svg>
""")
    write(path, take!(io))
end

svg_bar_log(joinpath(OUT, "profile_bar.svg");
    title = "Per-frame component times, physics vs diagnostics vs I/O",
    subtitle = "256³ limnickels, GTX 1080 FP64 — each component, log scale",
    bars = [("GPU physics stepping (~$(round(Int, d["steps_per_frame"])) steps)", physics, BLUE),
            ("frame_products (curls, |B|/|ω|, volren)", d["frame_products"], ORANGE),
            ("frame_scalars (moments, energies)", d["frame_scalars"], ORANGE),
            ("azimuthal spectrum (CPU)", d["spectrum_cpu"], ORANGE),
            ("2D panel render (CPU)", d["render2d_cpu"], ORANGE),
            ("PNG encode (panel + 3D)", d["png_panel"] + d["png_vol3d"], GREEN),
            ("volume dump (uint8, raycaster)", d["voldump"], GREEN)])

println("wrote out/plots/profile_{pie,bar}.svg")
println("physics=$(round(physics, digits=1))ms diag=$(round(diag, digits=1))ms io=$(round(io_, digits=1))ms  total=$(round(total, digits=1))ms")
println("physics fraction: $(round(100*physics/total, digits=2))%")
