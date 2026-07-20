# Plot the v3 anapole-null experiment (scripts/v3_anapole_null.jl, README §8
# experiment (a)): how completely the fractal coil's external *fields*
# cancel while its external *potentials* survive. Dependency-free SVG,
# same convention as plot_anapole.jl; palette from the dataviz skill.
#
#   julia --project=. scripts/plot_anapole_null.jl
#
# Reads  out/v3/anapole_null_{shells,link}.csv
# Writes out/plots/anapole_null_falloff.svg, out/plots/anapole_null_link.svg

const OUT = joinpath(@__DIR__, "..", "out", "plots")
mkpath(OUT)
const BLUE, RED, GREEN = "#2a78d6", "#e34948", "#008300"
const INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
const SURFACE = "#fcfcfb"

function read_csv(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    rows = [split(l, ",") for l in lines[2:end]]
    Dict(h => i == 1 || isletter(first(rows[1][i])) ?
             [r[i] for r in rows] : [parse(Float64, r[i]) for r in rows]
         for (i, h) in enumerate(header))
end

"Log-log (or log-linear) multi-series line plot, direct-labeled."
function svg_loglog(path, series; xlabel, ylabel, title, subtitle = "",
                    logx = true, W = 720, H = 460)
    ml, mr, mt, mb = 74, 24, 60, 52
    pw, ph = W - ml - mr, H - mt - mb
    xs = [s[1] for s in series]; ys = [s[2] for s in series]
    xmin = minimum(minimum, xs); xmax = maximum(maximum, xs)
    yall = reduce(vcat, ys)
    yall = [y for y in yall if y > 0]
    lx0, lx1 = logx ? (log10(xmin), log10(xmax)) : (xmin, xmax)
    ly0, ly1 = log10(minimum(yall)), log10(maximum(yall))
    pad = 0.08 * (ly1 - ly0); ly0 -= pad; ly1 += pad
    X(x) = ml + ((logx ? log10(x) : x) - lx0) / (lx1 - lx0) * pw
    Y(y) = mt + ph - (log10(max(y, 1e-300)) - ly0) / (ly1 - ly0) * ph
    io = IOBuffer()
    print(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$H" font-family="sans-serif">
<rect width="$W" height="$H" fill="$SURFACE"/>
<text x="$(W ÷ 2)" y="24" text-anchor="middle" font-size="16" fill="$INK">$title</text>
<text x="$(W ÷ 2)" y="42" text-anchor="middle" font-size="12" fill="$MUTED">$subtitle</text>
""")
    for e in ceil(Int, ly0):floor(Int, ly1)
        y = Y(10.0^e)
        print(io, """<line x1="$ml" y1="$y" x2="$(W - mr)" y2="$y" stroke="$GRID"/>
<text x="$(ml - 8)" y="$(y + 4)" text-anchor="end" font-size="11" fill="$MUTED">1e$e</text>
""")
    end
    if logx
        for e in ceil(Int, lx0):floor(Int, lx1)
            x = X(10.0^e)
            print(io, """<line x1="$x" y1="$mt" x2="$x" y2="$(mt + ph)" stroke="$GRID"/>
<text x="$x" y="$(mt + ph + 18)" text-anchor="middle" font-size="11" fill="$MUTED">1e$e</text>
""")
        end
    else
        for q in 0:5
            v = lx0 + q / 5 * (lx1 - lx0)
            x = X(v)
            print(io, """<line x1="$x" y1="$mt" x2="$x" y2="$(mt + ph)" stroke="$GRID"/>
<text x="$x" y="$(mt + ph + 18)" text-anchor="middle" font-size="11" fill="$MUTED">$(round(v, digits = 1))</text>
""")
        end
    end
    print(io, """<rect x="$ml" y="$mt" width="$pw" height="$ph" fill="none" stroke="$MUTED"/>
<text x="$(ml + pw ÷ 2)" y="$(H - 10)" text-anchor="middle" font-size="13" fill="$INK2">$xlabel</text>
<text x="18" y="$(mt + ph ÷ 2)" text-anchor="middle" font-size="13" fill="$INK2" transform="rotate(-90 18 $(mt + ph ÷ 2))">$ylabel</text>
""")
    for (i, (x, y, label, color, dash)) in enumerate(series)
        pts = join(["$(round(X(x[j]), digits = 1)),$(round(Y(y[j]), digits = 1))"
                    for j in eachindex(x) if y[j] > 0], " ")
        dd = dash ? " stroke-dasharray=\"6,4\"" : ""
        print(io, """<polyline points="$pts" fill="none" stroke="$color" stroke-width="2.5"$dd/>
<line x1="$(ml + 12)" y1="$(mt + 18i)" x2="$(ml + 40)" y2="$(mt + 18i)" stroke="$color" stroke-width="3"$dd/>
<text x="$(ml + 46)" y="$(mt + 18i + 4)" font-size="12" fill="$INK">$label</text>
""")
    end
    print(io, "</svg>\n")
    write(path, take!(io))
end

sh = read_csv(joinpath(@__DIR__, "..", "out", "v3", "anapole_null_shells.csv"))
lk = read_csv(joinpath(@__DIR__, "..", "out", "v3", "anapole_null_link.csv"))

pick(d, col, var) = [d[col][i] for i in eachindex(d["variant"]) if d["variant"][i] == var]
pick(d, col, var, kind) = [d[col][i] for i in eachindex(d["variant"])
                          if d["variant"][i] == var && d["kind"][i] == kind]

svg_loglog(joinpath(OUT, "anapole_null_falloff.svg"),
    [(pick(sh, "r", "k1r"), pick(sh, "meanB", "k1r"), "shell-mean |B|", BLUE, false),
     (pick(sh, "r", "k1r"), pick(sh, "meanA", "k1r"), "shell-mean |A| (Coulomb gauge)", RED, false)];
    xlabel = "distance from coil center r / R0  (log)",
    ylabel = "field / potential magnitude  (log)",
    title = "The anapole coil's far field vs far potential (k1r: 16-turn coil + counter-loop, |m| ≈ 0)",
    subtitle = "|B| falls off faster than |A| — the field cancels more completely than the potential")

svg_loglog(joinpath(OUT, "anapole_null_link.svg"),
    [(pick(lk, "s", "k1r", "link"), abs.(pick(lk, "phi_link", "k1r", "link")),
      "linked loop: |Φ_link| (through the hole)", BLUE, false),
     (pick(lk, "s", "k1r", "link"), pick(lk, "meanB_loop", "k1r", "link"),
      "linked loop: mean |B| on the loop", RED, false),
     (pick(lk, "s", "k1r", "ctrl"), abs.(pick(lk, "phi_link", "k1r", "ctrl")),
      "control loop (same size, y=3, unlinked): |Φ|", BLUE, true),
     (pick(lk, "s", "k1r", "ctrl"), pick(lk, "meanB_loop", "k1r", "ctrl"),
      "control loop: mean |B|", RED, true)];
    logx = true,
    xlabel = "loop radius s  (log; loop recedes to distance ~2s)",
    ylabel = "flux or field magnitude  (log)",
    title = "What an Aharonov–Bohm probe sees vs a magnetometer (k1r anapole coil)",
    subtitle = "linked flux Φ stays ≈ constant as the loop grows; B on the same loop decays like the far field")

# console summary for the docs
i15 = findlast(<=(16.0), pick(sh, "r", "k1r"))
println("k1r far shell (r≈", pick(sh, "r", "k1r")[i15], " R0): |B| = ",
        pick(sh, "meanB", "k1r")[i15], "  |A| = ", pick(sh, "meanA", "k1r")[i15])
svec = pick(lk, "s", "k1r", "link")
phivec = pick(lk, "phi_link", "k1r", "link")
bvec = pick(lk, "meanB_loop", "k1r", "link")
println("k1r linked loop, s: ", svec[1], " -> ", svec[end],
        "  Φ: ", phivec[1], " -> ", phivec[end], " (ratio ", phivec[end] / phivec[1], ")",
        "  meanB at s=", svec[end], ": ", bvec[end])
println("wrote out/plots/anapole_null_{falloff,link}.svg")
