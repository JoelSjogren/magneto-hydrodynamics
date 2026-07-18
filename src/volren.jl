# Orthographic emission–absorption volume rendering, dependency-free.
# One scalar field controls opacity (where the structure is), another
# controls color (what property it has there), composited front-to-back.

"""
    volume_render(color, opacity, box; azim=0.6, elev=0.45, res=384,
                  nstep=256, op_scale=8.0, clo=0.0, chi=1.0,
                  logcolor=false) -> Array{UInt8,3}

Raycast the cube: `opacity` (≥0) sets extinction per unit length
(scaled by `op_scale` relative to its maximum), `color` is normalized to
[clo, chi] (after log10 if `logcolor`) through the viridis map. Returns a
3×res×res RGB image (pass to `save_png`).
"""
function volume_render(color::Array{Float64,3}, opacity::Array{Float64,3},
                       box::Box; azim = 0.6, elev = 0.45, res = 384,
                       nstep = 256, op_scale = 18.0, op_gamma = 0.5,
                       clo = 0.0, chi = 1.0, logcolor = false, bg = 0.07)
    n = box.n
    L = n * box.dx
    opmax = maximum(opacity)
    opmax <= 0 && (opmax = 1.0)
    # opacity transfer: (op/opmax)^γ — γ < 1 lifts the diffuse mid-range so
    # weak extended fields stay visible (raw values left the limnickels
    # video almost fully transparent)
    κ = op_scale / L
    ca, sa = cos(azim), sin(azim)
    ce, se = cos(elev), sin(elev)
    # view basis: dir = into screen, (ex, ey) = screen axes (world coords)
    dir = (-ca * ce, -sa * ce, -se)
    ex = (-sa, ca, 0.0)
    ey = (-ca * se, -sa * se, ce)
    half = L / 2
    cx = box.lo + half
    diag = sqrt(3.0) * half
    ds = 2 * diag / nstep
    span = chi - clo > 0 ? chi - clo : 1.0

    img = Array{UInt8,3}(undef, 3, res, res)
    bg8 = UInt8(round(255 * bg))
    Threads.@threads for py in 1:res
        @inbounds for px in 1:res
            u = (2 * (px - 0.5) / res - 1) * diag
            v = (2 * (py - 0.5) / res - 1) * diag
            # start point: center + u·ex + v·ey − diag·dir
            x = cx + u * ex[1] + v * ey[1] - diag * dir[1]
            y = cx + u * ex[2] + v * ey[2] - diag * dir[2]
            z = cx + u * ex[3] + v * ey[3] - diag * dir[3]
            r = 0.0; g = 0.0; b = 0.0; T = 1.0
            for _ in 1:nstep
                x += dir[1] * ds; y += dir[2] * ds; z += dir[3] * ds
                fx = (x - box.lo) / box.dx - 0.5
                fy = (y - box.lo) / box.dx - 0.5
                fz = (z - box.lo) / box.dx - 0.5
                (0.0 <= fx <= n - 1.001 && 0.0 <= fy <= n - 1.001 &&
                 0.0 <= fz <= n - 1.001) || continue
                i0 = floor(Int, fx) + 1; tx = fx - (i0 - 1)
                j0 = floor(Int, fy) + 1; ty = fy - (j0 - 1)
                k0 = floor(Int, fz) + 1; tz = fz - (k0 - 1)
                # trilinear opacity
                op = 0.0; cv = 0.0
                for dk in 0:1, dj in 0:1, di in 0:1
                    w = (di == 1 ? tx : 1 - tx) * (dj == 1 ? ty : 1 - ty) *
                        (dk == 1 ? tz : 1 - tz)
                    op += w * opacity[i0+di, j0+dj, k0+dk]
                    cv += w * color[i0+di, j0+dj, k0+dk]
                end
                op <= 0 && continue
                a = 1 - exp(-κ * (op / opmax)^op_gamma * ds)
                a < 1e-4 && continue
                vv = logcolor ? log10(max(cv, 1e-300)) : cv
                vn = clamp((vv - clo) / span, 0.0, 1.0)
                cr, cg, cb = _colormap(vn)
                w = T * a
                r += w * cr; g += w * cg; b += w * cb
                T *= 1 - a
                T < 0.01 && break
            end
            row = res - py + 1
            img[1, px, row] = UInt8(round(clamp(r + T * 255bg, 0, 255)))
            img[2, px, row] = UInt8(round(clamp(g + T * 255bg, 0, 255)))
            img[3, px, row] = UInt8(round(clamp(b + T * 255bg, 0, 255)))
        end
    end
    img
end
