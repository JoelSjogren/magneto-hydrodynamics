# Component profile of a production-size GPU run: where does the wall time
# of one frame cycle go — physics stepping, frame products (curls +
# magnitudes + volume raycast), scalar reductions, the CPU-side azimuthal
# spectrum, or the (async in production) writer work: 2D panel render, PNG
# encodes, volume dumps? Checkpoint write is timed once, separately.
#
#   julia --project=. scripts/bench_gpu_profile.jl [ngrid] [ncycles]
#
# Setup mirrors the limnickels scenario of v2_pathB_selfassembly.jl at
# production parameters. Writes out/v2/profile_N<ngrid>.csv with per-frame
# mean milliseconds per component, and prints the same as a table.

using FractalToroid
using Printf
using CUDA

const NGRID = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 256
const NCYC = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3

CUDA.functional() || error("CUDA is not functional on this machine")
include(joinpath(dirname(pathof(FractalToroid)), "mhd_cuda.jl"))

const RING_R = 0.8
const FRAME_DT = 0.25
const UP = max(1, 512 ÷ NGRID)
const GAP = 8
const OUTDIR = joinpath(@__DIR__, "..", "out", "v2")
const SCRATCH = joinpath(OUTDIR, "profile_scratch")
mkpath(SCRATCH)

box = Box(NGRID, 2.0)
sim = MHDSim(box; cs = 0.55, eta = 2e-3, sponge_width = 6)
add_vortex_ring!(sim; R = RING_R, a = 0.2, z0 = -0.9, P0 = 0.40)
add_vortex_ring!(sim; R = RING_R, a = 0.2, z0 = +0.9, P0 = -0.40)
add_flux_ring!(sim; R = RING_R, a = 0.2, z0 = -0.9, A0 = 0.10, Bt0 = 0.0)
add_flux_ring!(sim; R = RING_R, a = 0.2, z0 = +0.9, A0 = 0.10, Bt0 = 0.0)
G = MHDCuda.to_gpu(sim)

ωx, ωy, ωz = curl_central(sim.S[MMX], sim.S[MMY], sim.S[MMZ],
                          box, sim.ip, sim.im)
const OM_HI = max(maximum(sqrt.(ωx .^ 2 .+ ωy .^ 2 .+ ωz .^ 2)), 1e-12) * 1.5
const LOGB_HI = log10(max(maximum(sqrt.(sim.S[MBX] .^ 2 .+
                                        sim.S[MBY] .^ 2 .+
                                        sim.S[MBZ] .^ 2)), 1e-12))

"Wall-clock ms of f() with device sync on both sides."
function tms(f)
    CUDA.synchronize()
    t0 = time_ns()
    r = f()
    CUDA.synchronize()
    r, (time_ns() - t0) / 1e6
end

"Copy of the v2 script's 2-row panel render (the writer's CPU raster work)."
function render2d(B, W)
    j0 = box.n ÷ 2
    npx = box.n * UP
    rgb = fill(0x12, 3, 2npx + GAP, 2npx + GAP)
    for (F, islog, sl, px, py) in ((B, true, :xz, 0, 0), (B, true, :xy, 1, 0),
                                   (W, false, :xz, 0, 1), (W, false, :xy, 1, 1))
        for b in 1:box.n, a in 1:box.n
            v = sl === :xz ? F[a, j0, b] : F[a, b, j0]
            vv = islog ?
                 clamp((log10(max(v, 1e-300)) - (LOGB_HI - 3)) / 3, 0, 1) :
                 clamp(v / OM_HI, 0, 1)
            c = FractalToroid._colormap(vv)
            for dj in 1:UP, di in 1:UP
                col = px * (npx + GAP) + (a - 1) * UP + di
                row = py * (npx + GAP) + (box.n - b) * UP + dj
                for ch in 1:3
                    rgb[ch, col, row] = c[ch]
                end
            end
        end
    end
    rgb
end

"Copy of the v2 script's uint8 volume quantize + write."
function dumpvol(B, W)
    qB = Vector{UInt8}(undef, length(B))
    qW = Vector{UInt8}(undef, length(W))
    @inbounds for i in eachindex(B)
        vb = (log10(max(B[i], 1e-300)) - (LOGB_HI - 3.5)) / 3.5
        qB[i] = round(UInt8, 255 * clamp(vb, 0.0, 1.0))
        qW[i] = round(UInt8, 255 * clamp(W[i] / OM_HI, 0.0, 1.0))
    end
    write(joinpath(SCRATCH, "B.bin"), qB)
    write(joinpath(SCRATCH, "W.bin"), qW)
end

# warm-up: compile every code path once
MHDCuda.gpu_step!(sim, G)
Bf, Wf, J2f, img3d = MHDCuda.frame_render_products!(G, sim; chi = OM_HI)
MHDCuda.gpu_frame_scalars!(G, sim)
MHDCuda.gpu_volume_render(G.S2[2], G.S2[1], G, box; res = 448, chi = OM_HI)
azimuthal_spectrum(J2f, box, RING_R)
render2d(Bf, Wf)
save_png(joinpath(SCRATCH, "panel.png"), render2d(Bf, Wf))
save_png(joinpath(SCRATCH, "vol3d.png"), img3d)
dumpvol(Bf, Wf)

acc = Dict{String,Float64}()
nsteps_tot = 0
for cyc in 1:NCYC
    global Bf, Wf, J2f, img3d, nsteps_tot
    target = sim.t + FRAME_DT
    nst = 0
    _, ms = tms() do
        while sim.t < target
            MHDCuda.gpu_step!(sim, G)
            nst += 1
        end
    end
    acc["step_stretch"] = get(acc, "step_stretch", 0.0) + ms
    nsteps_tot += nst

    (Bf, Wf, J2f, img3d), ms = tms() do
        MHDCuda.frame_render_products!(G, sim; chi = OM_HI)
    end
    acc["frame_products"] = get(acc, "frame_products", 0.0) + ms

    _, ms = tms() do
        MHDCuda.gpu_volume_render(G.S2[2], G.S2[1], G, box;
                                  res = 448, chi = OM_HI)
    end
    acc["volren_only"] = get(acc, "volren_only", 0.0) + ms

    _, ms = tms() do
        Float64.(Array(G.S2[1]))
    end
    acc["download_one"] = get(acc, "download_one", 0.0) + ms

    _, ms = tms() do
        MHDCuda.gpu_frame_scalars!(G, sim)
    end
    acc["frame_scalars"] = get(acc, "frame_scalars", 0.0) + ms

    _, ms = tms() do
        azimuthal_spectrum(J2f, box, RING_R)
    end
    acc["spectrum_cpu"] = get(acc, "spectrum_cpu", 0.0) + ms

    rgb, ms = tms() do
        render2d(Bf, Wf)
    end
    acc["render2d_cpu"] = get(acc, "render2d_cpu", 0.0) + ms

    _, ms = tms() do
        save_png(joinpath(SCRATCH, "panel.png"), rgb)
    end
    acc["png_panel"] = get(acc, "png_panel", 0.0) + ms

    _, ms = tms() do
        save_png(joinpath(SCRATCH, "vol3d.png"), img3d)
    end
    acc["png_vol3d"] = get(acc, "png_vol3d", 0.0) + ms

    _, ms = tms() do
        dumpvol(Bf, Wf)
    end
    acc["voldump"] = get(acc, "voldump", 0.0) + ms
    @printf("cycle %d/%d done (t=%.3f, %d steps)\n", cyc, NCYC, sim.t, nst)
end

_, ck_ms = tms() do
    MHDCuda.download!(sim, G)
    checkpoint_save(sim, SCRATCH)
end

rows = ["component,ms_per_frame"]
order = ["step_stretch", "frame_products", "volren_only", "download_one",
         "frame_scalars", "spectrum_cpu", "render2d_cpu", "png_panel",
         "png_vol3d", "voldump"]
println("\n== per-frame means over $NCYC cycles, $(NGRID)³ ==")
for k in order
    v = acc[k] / NCYC
    @printf("%-16s %10.1f ms\n", k, v)
    push!(rows, "$k,$(round(v, digits = 2))")
end
@printf("%-16s %10.1f ms  (once per run)\n", "checkpoint", ck_ms)
push!(rows, "checkpoint_once,$(round(ck_ms, digits = 2))")
push!(rows, "steps_per_frame,$(round(nsteps_tot / NCYC, digits = 2))")
write(joinpath(OUTDIR, "profile_N$(NGRID).csv"), join(rows, "\n") * "\n")
println("wrote ", joinpath(OUTDIR, "profile_N$(NGRID).csv"))
