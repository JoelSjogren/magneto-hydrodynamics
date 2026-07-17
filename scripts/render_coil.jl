# Render one still of the fractal "coil of coils" as a PNG: orthographic 3D
# projection, far-to-near painter's algorithm, viridis depth cueing.
#
#   julia --project=. scripts/render_coil.jl [k] [pixels] [azim_deg] [elev_deg]
#
# Output: out/render/coil_k<k>.png

using FractalToroid

const K = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2
const RES = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1200
const AZIM = deg2rad(length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 35.0)
const ELEV = deg2rad(length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 28.0)
const WINDINGS = [16, 8, 6]
const OUT = joinpath(@__DIR__, "..", "out", "render")

mkpath(OUT)

coil = fractal_coil(K; R0 = 1.0, ratio = 0.25, windings = WINDINGS[1:K],
                    ppt = 64, minpts = 4096)
P = coil.P
N = size(P, 2)
println("Rendering k=$K coil: $N points at $(RES)×$(RES) px")

# view rotation: azimuth about z, then elevation about x
ca, sa = cos(AZIM), sin(AZIM)
ce, se = cos(ELEV), sin(ELEV)
p_scr = Matrix{Float64}(undef, 3, N)      # (screen x, screen y, depth)
for i in 1:N
    x = ca * P[1, i] - sa * P[2, i]
    y = sa * P[1, i] + ca * P[2, i]
    z = P[3, i]
    p_scr[1, i] = x
    p_scr[2, i] = ce * y + se * z         # screen up
    p_scr[3, i] = -se * y + ce * z        # depth toward viewer
end

lim = 1.06 * maximum(abs, @view p_scr[1:2, :])
tosc(v) = (v + lim) / 2lim * (RES - 1) + 1
zlo, zhi = extrema(@view p_scr[3, :])

# far-to-near segment order
order = sortperm([0.5 * (p_scr[3, i] + p_scr[3, i % N + 1]) for i in 1:N])

img = fill(0.055f0, 3, RES, RES)          # dark background
"Stamp a filled disc, overwriting (painter's occlusion)."
function stamp!(img, x, y, r, rgb)
    for dj in -ceil(Int, r):ceil(Int, r), di in -ceil(Int, r):ceil(Int, r)
        di * di + dj * dj <= r * r || continue
        i = round(Int, x) + di
        j = round(Int, y) + dj
        (1 <= i <= RES && 1 <= j <= RES) || continue
        img[1, i, RES-j+1] = rgb[1]
        img[2, i, RES-j+1] = rgb[2]
        img[3, i, RES-j+1] = rgb[3]
    end
end

for s in order
    s2 = s % N + 1
    x1, y1, z1 = p_scr[1, s], p_scr[2, s], p_scr[3, s]
    x2, y2, z2 = p_scr[1, s2], p_scr[2, s2], p_scr[3, s2]
    t = clamp(((z1 + z2) / 2 - zlo) / (zhi - zlo), 0.0, 1.0)
    r8, g8, b8 = FractalToroid._colormap(0.15 + 0.8 * t)
    shade = Float32(0.45 + 0.55 * t)      # nearer = brighter
    rgb = (r8 / 255 * shade, g8 / 255 * shade, b8 / 255 * shade)
    rad = 0.0012 * RES * (0.7 + 0.9 * t)  # nearer = slightly thicker
    steps = max(2, ceil(Int, hypot(tosc(x2) - tosc(x1), tosc(y2) - tosc(y1))))
    for q in 0:steps
        f = q / steps
        stamp!(img, tosc(x1 + f * (x2 - x1)), tosc(y1 + f * (y2 - y1)), rad, rgb)
    end
end

out = joinpath(OUT, "coil_k$(K).png")
save_png(out, UInt8.(round.(clamp.(img, 0, 1) .* 255)))
println("Wrote $(abspath(out))")
