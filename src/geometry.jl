# Recursive "coil of coils" space curve, built with a parallel-transport
# (rotation-minimizing) frame. Winding numbers are integers so every level
# closes; the frame's holonomy around the closed loop is measured and
# distributed as a corrective twist so the offset curve closes too.

"""
A closed polyline: 3×N matrix of points; segment N connects back to point 1.
"""
struct Curve
    P::Matrix{Float64}
end

npoints(c::Curve) = size(c.P, 2)

@inline _norm3(x, y, z) = sqrt(x * x + y * y + z * z)

"""
    frames(P) -> (T, N, B)

Unit tangent plus a closed parallel-transport normal/binormal frame along the
closed polyline `P` (3×N). Holonomy is removed by distributing a compensating
twist, so column N+1 would coincide with column 1.
"""
function frames(P::Matrix{Float64})
    N = size(P, 2)
    T = zeros(3, N)
    @inbounds for i in 1:N
        ip = i == N ? 1 : i + 1
        im = i == 1 ? N : i - 1
        tx = P[1, ip] - P[1, im]
        ty = P[2, ip] - P[2, im]
        tz = P[3, ip] - P[3, im]
        s = 1 / _norm3(tx, ty, tz)
        T[1, i] = tx * s; T[2, i] = ty * s; T[3, i] = tz * s
    end

    Nrm = zeros(3, N)
    # seed normal: any unit vector ⊥ T₁
    let tx = T[1, 1], ty = T[2, 1], tz = T[3, 1]
        ax, ay, az = abs(tx) < 0.9 ? (1.0, 0.0, 0.0) : (0.0, 1.0, 0.0)
        d = ax * tx + ay * ty + az * tz
        nx, ny, nz = ax - d * tx, ay - d * ty, az - d * tz
        s = 1 / _norm3(nx, ny, nz)
        Nrm[1, 1] = nx * s; Nrm[2, 1] = ny * s; Nrm[3, 1] = nz * s
    end
    # transport by projecting the previous normal off the new tangent
    @inbounds for i in 2:N
        nx, ny, nz = Nrm[1, i-1], Nrm[2, i-1], Nrm[3, i-1]
        tx, ty, tz = T[1, i], T[2, i], T[3, i]
        d = nx * tx + ny * ty + nz * tz
        nx -= d * tx; ny -= d * ty; nz -= d * tz
        s = 1 / _norm3(nx, ny, nz)
        Nrm[1, i] = nx * s; Nrm[2, i] = ny * s; Nrm[3, i] = nz * s
    end
    # holonomy: transport once more from N back to 1 and compare
    hx, hy, hz = Nrm[1, N], Nrm[2, N], Nrm[3, N]
    let tx = T[1, 1], ty = T[2, 1], tz = T[3, 1]
        d = hx * tx + hy * ty + hz * tz
        hx -= d * tx; hy -= d * ty; hz -= d * tz
        s = 1 / _norm3(hx, hy, hz)
        hx *= s; hy *= s; hz *= s
    end
    n1x, n1y, n1z = Nrm[1, 1], Nrm[2, 1], Nrm[3, 1]
    b1x = T[2, 1] * n1z - T[3, 1] * n1y
    b1y = T[3, 1] * n1x - T[1, 1] * n1z
    b1z = T[1, 1] * n1y - T[2, 1] * n1x
    α = atan(hx * b1x + hy * b1y + hz * b1z, hx * n1x + hy * n1y + hz * n1z)

    B = zeros(3, N)
    @inbounds for i in 1:N
        # binormal, then untwist both by −α·(i−1)/N about T
        tx, ty, tz = T[1, i], T[2, i], T[3, i]
        nx, ny, nz = Nrm[1, i], Nrm[2, i], Nrm[3, i]
        bx = ty * nz - tz * ny
        by = tz * nx - tx * nz
        bz = tx * ny - ty * nx
        φ = -α * (i - 1) / N
        c, s = cos(φ), sin(φ)
        Nrm[1, i] = c * nx + s * bx; Nrm[2, i] = c * ny + s * by; Nrm[3, i] = c * nz + s * bz
        B[1, i] = -s * nx + c * bx; B[2, i] = -s * ny + c * by; B[3, i] = -s * nz + c * bz
    end
    T, Nrm, B
end

"""
    fractal_coil(levels; R0=1.0, ratio=0.25, windings=8, phase=0.0,
                 ppt=24, minpts=512) -> Curve

Level 0 is a circle of radius `R0` in the xy-plane. Each further level winds
`windings[l]` integer turns (per parent turn) around the parent tube at radius
`ratio[l] * (parent tube radius)`. `ppt` = sample points per innermost turn.
"""
function fractal_coil(levels::Int; R0 = 1.0, ratio = 0.25, windings = 8,
                      phase = 0.0, ppt = 24, minpts = 512)
    w = windings isa Integer ? fill(Int(windings), levels) : collect(Int, windings)
    r = ratio isa Real ? fill(Float64(ratio), levels) : collect(Float64, ratio)
    length(w) == levels || throw(ArgumentError("need $levels winding numbers"))
    length(r) == levels || throw(ArgumentError("need $levels ratios"))

    turns = 1
    for l in 1:levels
        turns *= w[l]
    end
    N = max(minpts, ppt * turns)
    θ = [2π * (i - 1) / N for i in 1:N]

    P = zeros(3, N)
    @inbounds for i in 1:N
        P[1, i] = R0 * cos(θ[i])
        P[2, i] = R0 * sin(θ[i])
    end

    a = R0
    rate = 1
    for l in 1:levels
        a *= r[l]
        rate *= w[l]
        _, Nrm, B = frames(P)
        Q = similar(P)
        @inbounds for i in 1:N
            c = cos(rate * θ[i] + phase)
            s = sin(rate * θ[i] + phase)
            for ax in 1:3
                Q[ax, i] = P[ax, i] + a * (Nrm[ax, i] * c + B[ax, i] * s)
            end
        end
        P = Q
    end
    Curve(P)
end

"""
    segments(c::Curve) -> (mid, dl)

Midpoints and vectors of the N closed segments, both 3×N.
"""
function segments(c::Curve)
    P = c.P
    N = size(P, 2)
    mid = similar(P)
    dl = similar(P)
    @inbounds for i in 1:N
        j = i == N ? 1 : i + 1
        for a in 1:3
            dl[a, i] = P[a, j] - P[a, i]
            mid[a, i] = 0.5 * (P[a, j] + P[a, i])
        end
    end
    mid, dl
end
