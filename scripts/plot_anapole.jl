# Plot the anapole self-assembly observations from the v2 96³ runs as SVG
# (dependency-free; text labels and crisp lines come free with SVG).
#
#   julia --project=. scripts/plot_anapole.jl
#
# Reads  out/v2/{limnickels,counterhel}_N96/timeseries.csv
# Writes out/plots/anapole_T.svg, out/plots/dipole_mz.svg

const OUT = joinpath(@__DIR__, "..", "out", "plots")
mkpath(OUT)

function read_csv(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    cols = Dict(h => Float64[] for h in header)
    for l in lines[2:end]
        for (h, v) in zip(header, split(l, ","))
            push!(cols[h], parse(Float64, v))
        end
    end
    cols
end

"Minimal SVG multi-series line plot with log-y option."
function svg_lineplot(path, series; xlabel = "t", ylabel = "", title = "",
                      logy = false, W = 760, H = 440)
    ml, mr, mt, mb = 78, 20, 44, 52          # margins
    pw, ph = W - ml - mr, H - mt - mb
    xs = [s[1] for s in series]; ys = [s[2] for s in series]
    xmin = minimum(minimum, xs); xmax = maximum(maximum, xs)
    yall = reduce(vcat, ys)
    logy && (yall = [y for y in yall if y > 0])
    tval = logy ? log10.(yall) : yall
    ymin, ymax = minimum(tval), maximum(tval)
    ymax == ymin && (ymax += 1)
    pad = 0.06 * (ymax - ymin); ymin -= pad; ymax += pad
    X(x) = ml + (x - xmin) / (xmax - xmin) * pw
    Y(y) = mt + ph - ((logy ? log10(max(y, 1e-300)) : y) - ymin) /
                    (ymax - ymin) * ph
    io = IOBuffer()
    print(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$H" font-family="sans-serif">
<rect width="$W" height="$H" fill="white"/>
<text x="$(W ÷ 2)" y="24" text-anchor="middle" font-size="16">$title</text>
""")
    # y ticks: decades if log, else 5 linear ticks
    if logy
        for e in ceil(Int, ymin):floor(Int, ymax)
            y = Y(10.0^e)
            print(io, """<line x1="$ml" y1="$y" x2="$(W - mr)" y2="$y" stroke="#ddd"/>
<text x="$(ml - 8)" y="$(y + 4)" text-anchor="end" font-size="12">1e$e</text>
""")
        end
    else
        for q in 0:4
            v = ymin + q / 4 * (ymax - ymin)
            y = Y(v)
            lbl = round(v, sigdigits = 2)
            print(io, """<line x1="$ml" y1="$y" x2="$(W - mr)" y2="$y" stroke="#ddd"/>
<text x="$(ml - 8)" y="$(y + 4)" text-anchor="end" font-size="12">$lbl</text>
""")
        end
    end
    for q in 0:5
        v = xmin + q / 5 * (xmax - xmin)
        x = X(v)
        print(io, """<line x1="$x" y1="$mt" x2="$x" y2="$(mt + ph)" stroke="#eee"/>
<text x="$x" y="$(H - mb + 18)" text-anchor="middle" font-size="12">$(round(v, digits = 1))</text>
""")
    end
    print(io, """<rect x="$ml" y="$mt" width="$pw" height="$ph" fill="none" stroke="#888"/>
<text x="$(ml + pw ÷ 2)" y="$(H - 12)" text-anchor="middle" font-size="13">$xlabel</text>
<text x="20" y="$(mt + ph ÷ 2)" text-anchor="middle" font-size="13" transform="rotate(-90 20 $(mt + ph ÷ 2))">$ylabel</text>
""")
    for (i, (x, y, label, color)) in enumerate(series)
        pts = join(["$(round(X(x[j]), digits = 1)),$(round(Y(y[j]), digits = 1))"
                    for j in eachindex(x) if !logy || y[j] > 0], " ")
        print(io, """<polyline points="$pts" fill="none" stroke="$color" stroke-width="2"/>
<line x1="$(ml + 12)" y1="$(mt + 16i)" x2="$(ml + 40)" y2="$(mt + 16i)" stroke="$color" stroke-width="3"/>
<text x="$(ml + 46)" y="$(mt + 16i + 4)" font-size="12">$label</text>
""")
    end
    print(io, "</svg>\n")
    write(path, take!(io))
    path
end

lim = read_csv(joinpath(@__DIR__, "..", "out", "v2", "limnickels_N96", "timeseries.csv"))
cnh = read_csv(joinpath(@__DIR__, "..", "out", "v2", "counterhel_N96", "timeseries.csv"))

svg_lineplot(joinpath(OUT, "anapole_T.svg"),
    [(lim["t"], lim["Tnorm"], "limnickels (ring collision)", "#1f77b4"),
     (cnh["t"], cnh["Tnorm"], "counterhel (merging)", "#d62728")];
    ylabel = "|T|  (anapole moment)", title =
    "Spontaneous anapole growth from anapole-free initial conditions (96³, 2% noise)",
    logy = true)

svg_lineplot(joinpath(OUT, "dipole_mz.svg"),
    [(lim["t"], abs.(lim["mz"]), "limnickels (ring collision)", "#1f77b4"),
     (cnh["t"], abs.(cnh["mz"]), "counterhel (merging)", "#d62728")];
    ylabel = "|m_z|  (magnetic dipole)", title =
    "Seed field wound into a net dipole (96³, 2% noise)", logy = true)

println("Wrote ", joinpath(OUT, "anapole_T.svg"), " and dipole_mz.svg")
