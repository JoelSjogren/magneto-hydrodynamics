# Magnetostatics of a line current (Biot–Savart), multipole moments, and
# deposition of the line current onto a grid as a Gaussian tube.
# Dimensionless units: μ0 = ε0 = c = 1.

"""
Uniform cubic grid: `n`³ cells over [lo, lo+n*dx)³. Node i sits at
lo + (i-1)*dx; cell centers at lo + (i-1/2)*dx.
"""
struct Box
    n::Int
    dx::Float64
    lo::Float64
end

Box(n::Int, halfwidth::Real) = Box(n, 2halfwidth / n, -Float64(halfwidth))

node(b::Box, i::Int) = b.lo + (i - 1) * b.dx
center(b::Box, i::Int) = b.lo + (i - 0.5) * b.dx
cellvol(b::Box) = b.dx^3

"""
    biot_savart(mid, dl, box; I=1.0, eps) -> (Bx, By, Bz)

B(x) = (I/4π) Σ dl × (x−r) / (|x−r|² + eps²)^{3/2}, sampled at cell centers.
`eps` regularizes the on-wire singularity (default: one cell size).
"""
function biot_savart(mid::Matrix{Float64}, dl::Matrix{Float64}, box::Box;
                     I = 1.0, eps = box.dx)
    n = box.n
    Bx = zeros(n, n, n); By = zeros(n, n, n); Bz = zeros(n, n, n)
    ns = size(mid, 2)
    pre = I / (4π)
    e2 = eps * eps
    Threads.@threads for k in 1:n
        zk = center(box, k)
        @inbounds for j in 1:n
            yj = center(box, j)
            for i in 1:n
                xi = center(box, i)
                bx = 0.0; by = 0.0; bz = 0.0
                for s in 1:ns
                    rx = xi - mid[1, s]
                    ry = yj - mid[2, s]
                    rz = zk - mid[3, s]
                    r2 = rx * rx + ry * ry + rz * rz + e2
                    invr3 = 1 / (r2 * sqrt(r2))
                    bx += (dl[2, s] * rz - dl[3, s] * ry) * invr3
                    by += (dl[3, s] * rx - dl[1, s] * rz) * invr3
                    bz += (dl[1, s] * ry - dl[2, s] * rx) * invr3
                end
                Bx[i, j, k] = pre * bx
                By[i, j, k] = pre * by
                Bz[i, j, k] = pre * bz
            end
        end
    end
    Bx, By, Bz
end

"""
    biot_savart_component(mid, dl, box, comp; off=(0,0,0), I=1.0, eps)

One component of B sampled on the grid whose point (i,j,k) sits at
lo + (i-1+off)*dx per axis (for staggered/Yee sampling).
"""
function biot_savart_component(mid::Matrix{Float64}, dl::Matrix{Float64},
                               box::Box, comp::Int;
                               off = (0.0, 0.0, 0.0), I = 1.0, eps = box.dx)
    n = box.n
    B = zeros(n, n, n)
    ns = size(mid, 2)
    pre = I / (4π)
    e2 = eps * eps
    c1 = comp % 3 + 1        # cross product: B_c = dl_{c1} r_{c2} − dl_{c2} r_{c1}
    c2 = (comp + 1) % 3 + 1
    Threads.@threads for k in 1:n
        zk = box.lo + (k - 1 + off[3]) * box.dx
        @inbounds for j in 1:n
            yj = box.lo + (j - 1 + off[2]) * box.dx
            for i in 1:n
                xi = box.lo + (i - 1 + off[1]) * box.dx
                acc = 0.0
                for s in 1:ns
                    r1 = (xi, yj, zk)[c1] - mid[c1, s]
                    r2c = (xi, yj, zk)[c2] - mid[c2, s]
                    rx = xi - mid[1, s]
                    ry = yj - mid[2, s]
                    rz = zk - mid[3, s]
                    rr = rx * rx + ry * ry + rz * rz + e2
                    acc += (dl[c1, s] * r2c - dl[c2, s] * r1) / (rr * sqrt(rr))
                end
                B[i, j, k] = pre * acc
            end
        end
    end
    B
end

"""
    biot_savart_yee(mid, dl, box; I=1.0, eps) -> (Bx, By, Bz)

B on the Yee face positions: Bx(i,j+½,k+½), By(i+½,j,k+½), Bz(i+½,j+½,k).
"""
biot_savart_yee(mid, dl, box::Box; I = 1.0, eps = box.dx) = (
    biot_savart_component(mid, dl, box, 1; off = (0.0, 0.5, 0.5), I, eps),
    biot_savart_component(mid, dl, box, 2; off = (0.5, 0.0, 0.5), I, eps),
    biot_savart_component(mid, dl, box, 3; off = (0.5, 0.5, 0.0), I, eps),
)

"""
    vector_potential(mid, dl, box; I=1.0, eps) -> (Ax, Ay, Az)

A(x) = (I/4π) Σ dl / sqrt(|x−r|² + eps²), sampled at cell centers.
"""
function vector_potential(mid::Matrix{Float64}, dl::Matrix{Float64}, box::Box;
                          I = 1.0, eps = box.dx)
    n = box.n
    Ax = zeros(n, n, n); Ay = zeros(n, n, n); Az = zeros(n, n, n)
    ns = size(mid, 2)
    pre = I / (4π)
    e2 = eps * eps
    Threads.@threads for k in 1:n
        zk = center(box, k)
        @inbounds for j in 1:n
            yj = center(box, j)
            for i in 1:n
                xi = center(box, i)
                ax = 0.0; ay = 0.0; az = 0.0
                for s in 1:ns
                    rx = xi - mid[1, s]
                    ry = yj - mid[2, s]
                    rz = zk - mid[3, s]
                    invr = 1 / sqrt(rx * rx + ry * ry + rz * rz + e2)
                    ax += dl[1, s] * invr
                    ay += dl[2, s] * invr
                    az += dl[3, s] * invr
                end
                Ax[i, j, k] = pre * ax
                Ay[i, j, k] = pre * ay
                Az[i, j, k] = pre * az
            end
        end
    end
    Ax, Ay, Az
end

"""
    current_moments(mid, dl; I=1.0) -> (m, T)

Magnetic dipole m = (I/2) Σ r × dl and toroidal (anapole) moment
T = (I/10) Σ [(r·dl) r − 2 r² dl] of the closed line current.
"""
function current_moments(mid::Matrix{Float64}, dl::Matrix{Float64}; I = 1.0)
    m = zeros(3)
    T = zeros(3)
    @inbounds for s in 1:size(mid, 2)
        rx, ry, rz = mid[1, s], mid[2, s], mid[3, s]
        lx, ly, lz = dl[1, s], dl[2, s], dl[3, s]
        m[1] += ry * lz - rz * ly
        m[2] += rz * lx - rx * lz
        m[3] += rx * ly - ry * lx
        rdl = rx * lx + ry * ly + rz * lz
        r2 = rx * rx + ry * ry + rz * rz
        T[1] += rdl * rx - 2 * r2 * lx
        T[2] += rdl * ry - 2 * r2 * ly
        T[3] += rdl * rz - 2 * r2 * lz
    end
    (I / 2) .* m, (I / 10) .* T
end

"""
    helicity(Ax, Ay, Az, Bx, By, Bz, box) -> H = ∫ A·B dV
"""
function helicity(Ax, Ay, Az, Bx, By, Bz, box::Box)
    h = 0.0
    @inbounds for idx in eachindex(Bx)
        h += Ax[idx] * Bx[idx] + Ay[idx] * By[idx] + Az[idx] * Bz[idx]
    end
    h * cellvol(box)
end

"""
    shell_profile(Bx, By, Bz, box; nbins=32) -> (r, mean|B|)

Mean |B| over spherical shells around the origin — the "synthetic
magnetometer" radial profile.
"""
function shell_profile(Bx, By, Bz, box::Box; nbins = 32)
    n = box.n
    rmax = -box.lo
    acc = zeros(nbins)
    cnt = zeros(Int, nbins)
    @inbounds for k in 1:n, j in 1:n, i in 1:n
        r = sqrt(center(box, i)^2 + center(box, j)^2 + center(box, k)^2)
        b = floor(Int, r / rmax * nbins) + 1
        1 <= b <= nbins || continue
        acc[b] += sqrt(Bx[i, j, k]^2 + By[i, j, k]^2 + Bz[i, j, k]^2)
        cnt[b] += 1
    end
    r = [(b - 0.5) * rmax / nbins for b in 1:nbins]
    prof = [cnt[b] > 0 ? acc[b] / cnt[b] : NaN for b in 1:nbins]
    r, prof
end

"""
    splat_current!(Jx, Jy, Jz, mid, dl, box; I=1.0, sigma)

Deposit the line current as a Gaussian tube of radius `sigma` onto the three
staggered (Yee edge) component grids: Jx at (i+½,j,k), Jy at (i,j+½,k),
Jz at (i,j,k+½). ∫J dV = I ∮ dl by construction.
"""
function splat_current!(Jx, Jy, Jz, mid::Matrix{Float64}, dl::Matrix{Float64},
                        box::Box; I = 1.0, sigma = 2 * box.dx, stagger = :yee)
    n = box.n
    w = ceil(Int, 3.5 * sigma / box.dx)
    norm = I / ((2π)^1.5 * sigma^3)
    inv2s2 = 1 / (2 * sigma^2)
    h = stagger === :yee ? 0.5 : 0.0
    for s in 1:size(mid, 2)
        rx, ry, rz = mid[1, s], mid[2, s], mid[3, s]
        # component offsets: Yee edges (½,0,0),(0,½,0),(0,0,½) or nodes (0,0,0)
        for (J, comp, ox, oy, oz) in ((Jx, 1, h, 0.0, 0.0),
                                      (Jy, 2, 0.0, h, 0.0),
                                      (Jz, 3, 0.0, 0.0, h))
            amp = norm * dl[comp, s]
            amp == 0 && continue
            ic = round(Int, (rx - box.lo) / box.dx - ox) + 1
            jc = round(Int, (ry - box.lo) / box.dx - oy) + 1
            kc = round(Int, (rz - box.lo) / box.dx - oz) + 1
            @inbounds for k in kc-w:kc+w
                (1 <= k <= n) || continue
                dz = rz - (box.lo + (k - 1 + oz) * box.dx)
                for j in jc-w:jc+w
                    (1 <= j <= n) || continue
                    dy = ry - (box.lo + (j - 1 + oy) * box.dx)
                    for i in ic-w:ic+w
                        (1 <= i <= n) || continue
                        dx_ = rx - (box.lo + (i - 1 + ox) * box.dx)
                        J[i, j, k] += amp * exp(-(dx_^2 + dy^2 + dz^2) * inv2s2)
                    end
                end
            end
        end
    end
    Jx, Jy, Jz
end
