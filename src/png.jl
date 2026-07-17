# Minimal dependency-free PNG writer (truecolor, zlib "stored" blocks) plus a
# viridis-style colormap, so field slices can be saved without any plotting
# package.

const _CRC_TABLE = let t = Vector{UInt32}(undef, 256)
    for n in 0:255
        c = UInt32(n)
        for _ in 1:8
            c = (c & 0x1) == 0x1 ? (0xedb88320 ⊻ (c >> 1)) : (c >> 1)
        end
        t[n+1] = c
    end
    t
end

function _crc32(data::Vector{UInt8})
    c = 0xffffffff
    @inbounds for b in data
        c = _CRC_TABLE[((c ⊻ b) & 0xff)+1] ⊻ (c >> 8)
    end
    c ⊻ 0xffffffff
end

function _adler32(data::Vector{UInt8})
    a = UInt32(1); b = UInt32(0)
    @inbounds for x in data
        a = (a + x) % 65521
        b = (b + a) % 65521
    end
    (b << 16) | a
end

_be32(x::Integer) = UInt8[(x >> 24) & 0xff, (x >> 16) & 0xff, (x >> 8) & 0xff, x & 0xff]

function _chunk(io, tag::String, data::Vector{UInt8})
    write(io, _be32(length(data)))
    body = vcat(Vector{UInt8}(tag), data)
    write(io, body)
    write(io, _be32(_crc32(body)))
end

function _zlib_stored(raw::Vector{UInt8})
    out = UInt8[0x78, 0x01]
    i = 1
    while i <= length(raw) || isempty(raw)
        len = min(65535, length(raw) - i + 1)
        final = (i + len > length(raw)) ? 0x01 : 0x00
        push!(out, final)
        append!(out, UInt8[len & 0xff, (len >> 8) & 0xff,
                           ~UInt8(len & 0xff), ~UInt8((len >> 8) & 0xff)])
        append!(out, @view raw[i:i+len-1])
        i += len
        isempty(raw) && break
    end
    append!(out, _be32(_adler32(raw)))
    out
end

"""
    save_png(path, rgb::Array{UInt8,3})

Write an 8-bit truecolor PNG. `rgb` is 3×W×H (channel, column, row).
"""
function save_png(path::AbstractString, rgb::Array{UInt8,3})
    _, w, h = size(rgb)
    raw = Vector{UInt8}(undef, h * (1 + 3w))
    p = 1
    @inbounds for row in 1:h
        raw[p] = 0x00  # filter: none
        p += 1
        for col in 1:w
            raw[p] = rgb[1, col, row]
            raw[p+1] = rgb[2, col, row]
            raw[p+2] = rgb[3, col, row]
            p += 3
        end
    end
    open(path, "w") do io
        write(io, UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a])
        _chunk(io, "IHDR", vcat(_be32(w), _be32(h), UInt8[8, 2, 0, 0, 0]))
        _chunk(io, "IDAT", _zlib_stored(raw))
        _chunk(io, "IEND", UInt8[])
    end
    path
end

const _VIRIDIS = [
    0.267 0.005 0.329; 0.283 0.141 0.458; 0.254 0.265 0.530;
    0.207 0.372 0.553; 0.164 0.471 0.558; 0.128 0.567 0.551;
    0.135 0.659 0.518; 0.267 0.749 0.441; 0.478 0.821 0.318;
    0.741 0.873 0.150; 0.993 0.906 0.144
]

function _colormap(v::Float64)
    v = clamp(v, 0.0, 1.0) * (size(_VIRIDIS, 1) - 1)
    i = min(floor(Int, v), size(_VIRIDIS, 1) - 2)
    f = v - i
    r = _VIRIDIS[i+1, 1] * (1 - f) + _VIRIDIS[i+2, 1] * f
    g = _VIRIDIS[i+1, 2] * (1 - f) + _VIRIDIS[i+2, 2] * f
    b = _VIRIDIS[i+1, 3] * (1 - f) + _VIRIDIS[i+2, 3] * f
    UInt8(round(255r)), UInt8(round(255g)), UInt8(round(255b))
end

"""
    heatmap_png(path, A; logscale=false, upscale=1)

Save matrix `A` as a viridis heatmap. `A[i,j]` maps to pixel column i, row j
with row 1 at the bottom (math orientation). `upscale` repeats pixels.
"""
function heatmap_png(path::AbstractString, A::AbstractMatrix;
                     logscale = false, upscale = 1)
    V = float.(A)
    if logscale
        floor_ = maximum(V) * 1e-6
        V = log10.(max.(V, floor_))
    end
    lo, hi = extrema(V)
    span = hi > lo ? hi - lo : 1.0
    w, h = size(V)
    u = max(1, upscale)
    rgb = Array{UInt8,3}(undef, 3, w * u, h * u)
    for j in 1:h, i in 1:w
        r, g, b = _colormap((V[i, j] - lo) / span)
        for dj in 1:u, di in 1:u
            col = (i - 1) * u + di
            row = (h - j) * u + dj   # flip: row 1 of image = top
            rgb[1, col, row] = r
            rgb[2, col, row] = g
            rgb[3, col, row] = b
        end
    end
    save_png(path, rgb)
end
