# Small fitting helpers for the run diagnostics.

"""
    efold_time(t, e; tail=0.5) -> τ

E-folding decay time of a positive time series, from a least-squares fit of
log(e) over the last `tail` fraction of samples. Returns Inf if not decaying.
"""
function efold_time(t::AbstractVector, e::AbstractVector; tail = 0.5)
    i0 = max(1, round(Int, length(t) * (1 - tail)))
    ts = Float64[]
    ls = Float64[]
    for i in i0:length(t)
        e[i] > 0 || continue
        push!(ts, t[i])
        push!(ls, log(e[i]))
    end
    length(ts) < 2 && return Inf
    tm = sum(ts) / length(ts)
    lm = sum(ls) / length(ls)
    num = sum((ts .- tm) .* (ls .- lm))
    den = sum(abs2, ts .- tm)
    slope = num / den
    slope < 0 ? -1 / slope : Inf
end

"""
    powerlaw_slope(r, y; rmin, rmax) -> d log y / d log r

Log-log slope over [rmin, rmax] (e.g. −3 for a dipole far field).
"""
function powerlaw_slope(r::AbstractVector, y::AbstractVector;
                        rmin = -Inf, rmax = Inf)
    xs = Float64[]
    ys = Float64[]
    for i in eachindex(r)
        (rmin <= r[i] <= rmax && isfinite(y[i]) && y[i] > 0) || continue
        push!(xs, log(r[i]))
        push!(ys, log(y[i]))
    end
    length(xs) < 2 && return NaN
    xm = sum(xs) / length(xs)
    ym = sum(ys) / length(ys)
    sum((xs .- xm) .* (ys .- ym)) / sum(abs2, xs .- xm)
end
