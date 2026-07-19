# v2 path B — self-assembly from two tori (README §7.4). Evolve a pair of
# rings and watch for spontaneous structure: anapole growth |T|(t) and
# azimuthal breakup mode counts ("ring of rings").
#
#   julia -t auto --project=. scripts/v2_pathB_selfassembly.jl [scenario] [ngrid] [t_end] [resume] [gpu]
#
# The flag "gpu" (any position) runs the solver on the GPU via CUDA, FP64;
# "gpu32" instead keeps the device state in FP32 (fast on consumer cards,
# trajectory identical only statistically). Both error out if no functional
# CUDA device is available. Without a flag the run is always CPU.
# Diagnostics/rendering stay on the CPU (FP64) either way.
#
# scenarios: counterhel — co-current, counter-helicity magnetic ring pair
#            opposed    — anti-parallel ring currents forced together
#            limnickels — counter-rotating vortex-ring collision with a
#                          weak frozen-in seed field
#
# Checkpointing: the full state is saved on normal completion AND on
# Ctrl+C / SIGINT, to <outdir>/checkpoint.{bin,json}. Pass "resume" as the
# 4th argument (with a possibly larger t_end) to continue a run.
# Outputs: out/v2/<scenario>_N<ngrid>/

using FractalToroid
using Printf
using Random

Base.exit_on_sigint(false)   # deliver SIGINT as InterruptException

const GPU32 = "gpu32" in ARGS
const GPU = GPU32 || "gpu" in ARGS
# seed=N sets the symmetry-breaking RNG seed (default 1234). Non-default
# seeds write to a "<scenario>_N<n>_seedN" dir so an ensemble doesn't
# collide with the canonical run.
const SEED = let i = findfirst(a -> startswith(a, "seed="), ARGS)
    i === nothing ? 1234 : parse(Int, ARGS[i][6:end])
end
const POSARGS =
    filter(a -> !(a in ("gpu", "gpu32")) && !startswith(a, "seed="), ARGS)
const SCEN = length(POSARGS) >= 1 ? POSARGS[1] : "counterhel"
const NGRID = length(POSARGS) >= 2 ? parse(Int, POSARGS[2]) : 64
const T_END = length(POSARGS) >= 3 ? parse(Float64, POSARGS[3]) : 25.0
const RESUME = length(POSARGS) >= 4 && POSARGS[4] == "resume"

if GPU
    @eval using CUDA
    CUDA.functional() ||
        error("gpu flag given but CUDA is not functional on this machine")
    include(joinpath(dirname(pathof(FractalToroid)), "mhd_cuda.jl"))
end
const HALF = 2.0
const R = 0.8
const A = 0.2
const D = 0.7               # ring half-separation
const ETA = 2e-3            # Lundquist S = 1/η = 500 (unit length)
const CS = 0.55
const SPONGE = 6
const FRAME_DT = 0.25
const UP = max(1, 512 ÷ NGRID)
const GAP = 8
const OUT = joinpath(@__DIR__, "..", "out", "v2",
                     "$(SCEN)_N$(NGRID)" * (SEED == 1234 ? "" : "_seed$(SEED)"))

box = Box(NGRID, HALF)
sim = MHDSim(box; cs = CS, eta = ETA, sponge_width = SPONGE)

ckmeta = Dict{String,Float64}()
if RESUME
    ckmeta = checkpoint_load!(sim, OUT)
    println("resumed $(SCEN)_N$(NGRID) at t = $(round(sim.t, digits = 3))")
else
    mkpath(OUT)
    for d in ("frames", "frames3d", "vol")
        rm(joinpath(OUT, d); force = true, recursive = true)
        mkpath(joinpath(OUT, d))
    end
    if SCEN == "counterhel"
        add_flux_ring!(sim; R, a = A, z0 = -D, A0 = 0.6, Bt0 = 0.4)
        add_flux_ring!(sim; R, a = A, z0 = +D, A0 = 0.6, Bt0 = -0.4)
        add_vortex_ring!(sim; R, a = A, z0 = -D, P0 = 0.10)
        add_vortex_ring!(sim; R, a = A, z0 = +D, P0 = -0.10)
    elseif SCEN == "opposed"
        add_flux_ring!(sim; R, a = A, z0 = -D, A0 = 0.6, Bt0 = 0.0)
        add_flux_ring!(sim; R, a = A, z0 = +D, A0 = -0.6, Bt0 = 0.0)
        add_vortex_ring!(sim; R, a = A, z0 = -D, P0 = 0.30)
        add_vortex_ring!(sim; R, a = A, z0 = +D, P0 = -0.30)
    elseif SCEN == "limnickels"
        add_vortex_ring!(sim; R, a = A, z0 = -0.9, P0 = 0.40)
        add_vortex_ring!(sim; R, a = A, z0 = +0.9, P0 = -0.40)
        add_flux_ring!(sim; R, a = A, z0 = -0.9, A0 = 0.10, Bt0 = 0.0)
        add_flux_ring!(sim; R, a = A, z0 = +0.9, A0 = 0.10, Bt0 = 0.0)
    else
        error("unknown scenario $SCEN")
    end
    # symmetry-breaking noise (axisymmetric ICs on a clean grid leave
    # azimuthal instabilities nothing physical to grow from); applied only
    # at t = 0, so resumed runs continue the same trajectory
    let rng = Xoshiro(SEED)
        amp = 0.02 * maximum(abs, sim.S[MMZ])
        for q in (MMX, MMY, MMZ)
            Aq = sim.S[q]
            for idx in eachindex(Aq)
                Aq[idx] += amp * (2rand(rng) - 1)
            end
        end
    end
end

const GSIM = GPU ? MHDCuda.to_gpu(sim; T = GPU32 ? Float32 : Float64) : nothing

_bmag() = sqrt.(sim.S[MBX] .^ 2 .+ sim.S[MBY] .^ 2 .+ sim.S[MBZ] .^ 2)
function _omag()
    ωx, ωy, ωz = curl_central(sim.S[MMX], sim.S[MMY], sim.S[MMZ],
                              box, sim.ip, sim.im)
    sqrt.(ωx .^ 2 .+ ωy .^ 2 .+ ωz .^ 2)
end

# fixed render/quantization scales: from the initial state on a fresh run,
# from the checkpoint on resume (so appended frames stay consistent)
const LOGB_HI = RESUME ? ckmeta["logb_hi"] :
                log10(max(maximum(_bmag()), 1e-12))
const OM_HI = RESUME ? ckmeta["om_hi"] : max(maximum(_omag()), 1e-12) * 1.5
nframe = RESUME ? Int(ckmeta["nframe"]) : 0
tnext = RESUME ? ckmeta["tnext"] : FRAME_DT
rows = RESUME ? readlines(joinpath(OUT, "timeseries.csv")) :
       ["t,E_kin,E_mag,mz,Tnorm,Tz,mode_max,mode_amp"]

"4-panel frame: rows = log|B|, |ω|; columns = xz slice, xy (z=0) slice."
function render_frame(nf, B, W)
    j0 = box.n ÷ 2
    npx = box.n * UP
    side = 2npx + GAP
    rgb = fill(0x12, 3, side, side)
    panels = ((B, true, :xz, 0, 0), (B, true, :xy, 1, 0),
              (W, false, :xz, 0, 1), (W, false, :xy, 1, 1))
    for (F, islog, sl, px, py) in panels
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
    save_png(joinpath(OUT, "frames", "frame_$(lpad(nf, 5, '0')).png"), rgb)
end

"3D volume render: opacity = |B| (structure), color = |ω| (flow state)."
render_frame3d(nf, img) =
    save_png(joinpath(OUT, "frames3d", "frame_$(lpad(nf, 5, '0')).png"), img)

"uint8 volume dumps for the interactive raycaster (gitignored)."
function dump_volumes(nf, B, W)
    tag = lpad(nf, 5, '0')
    qB = Vector{UInt8}(undef, length(B))
    qW = Vector{UInt8}(undef, length(W))
    @inbounds for i in eachindex(B)
        vb = (log10(max(B[i], 1e-300)) - (LOGB_HI - 3.5)) / 3.5
        qB[i] = round(UInt8, 255 * clamp(vb, 0.0, 1.0))
        qW[i] = round(UInt8, 255 * clamp(W[i] / OM_HI, 0.0, 1.0))
    end
    write(joinpath(OUT, "vol", "B_$tag.bin"), qB)
    write(joinpath(OUT, "vol", "W_$tag.bin"), qW)
end

# On GPU runs the file work (panel render, PNG encode, volume dumps) runs
# on a spawned CPU task and overlaps the next stepping stretch; one frame
# in flight at a time. img3d === nothing → CPU path, fully synchronous.
const WRITER = Ref{Any}(nothing)
flush_writer() = (WRITER[] !== nothing && wait(WRITER[]); WRITER[] = nothing)

function emit_frame(B, W, img3d)
    global nframe
    nf = nframe
    if img3d === nothing
        render_frame(nf, B, W)
        render_frame3d(nf, volume_render(W, B, box; res = 448, chi = OM_HI,
                                         azim = 0.6, elev = 0.45))
        dump_volumes(nf, B, W)
    else
        flush_writer()
        WRITER[] = Threads.@spawn begin
            render_frame(nf, B, W)
            render_frame3d(nf, img3d)
            dump_volumes(nf, B, W)
        end
    end
    nframe += 1
end

"Write timeseries + checkpoint (normal completion and SIGINT both land here)."
function finish(status)
    flush_writer()
    GPU && MHDCuda.download!(sim, GSIM)
    write(joinpath(OUT, "timeseries.csv"), join(rows, "\n") * "\n")
    checkpoint_save(sim, OUT;
                    extra = Dict("logb_hi" => LOGB_HI, "om_hi" => OM_HI,
                                 "nframe" => nframe, "tnext" => tnext))
    m1, T1 = grid_moments(sim)
    @printf("[%s] t=%.3f: E_kin=%.4g E_mag=%.4g |T|=%.4g — checkpoint saved\n",
            status, sim.t, mhd_kinetic_energy(sim), mhd_magnetic_energy(sim),
            sqrt(sum(abs2, T1)))
end

if !RESUME
    write(joinpath(OUT, "vol", "meta.json"),
          """{"n": $(box.n), "fields": ["B", "W"],
              "B": {"scale": "log10", "lo": $(LOGB_HI - 3.5), "hi": $LOGB_HI},
              "W": {"scale": "linear", "lo": 0.0, "hi": $OM_HI}}""")
    emit_frame(_bmag(), _omag(), nothing)
end
write(joinpath(@__DIR__, "..", "out", "v2", "CURRENT"), "$(SCEN)_N$(NGRID)")

println("v2 path B [$SCEN] $(NGRID)³, t: $(round(sim.t, digits=2)) → $T_END, ",
        "S=$(1/ETA), ",
        GPU ? "CUDA $(GPU32 ? "FP32" : "FP64") ($(CUDA.name(CUDA.device())))" :
              "$(Threads.nthreads()) threads")
try
    while sim.t < T_END
        global tnext
        GPU ? MHDCuda.gpu_step!(sim, GSIM) : mhd_step!(sim)
        if sim.t >= tnext
            if GPU
                Bf, Wf, J2f, img3d =
                    MHDCuda.frame_render_products!(GSIM, sim; chi = OM_HI)
                m, T, ekin, emag = MHDCuda.gpu_frame_scalars!(GSIM, sim)
            else
                Bf = _bmag(); Wf = _omag()
                Jx, Jy, Jz = curl_central(sim.S[MBX], sim.S[MBY], sim.S[MBZ],
                                          box, sim.ip, sim.im)
                J2f = Jx .^ 2 .+ Jy .^ 2 .+ Jz .^ 2
                m, T = grid_moments(sim)
                ekin = mhd_kinetic_energy(sim)
                emag = mhd_magnetic_energy(sim)
                img3d = nothing
            end
            spec = azimuthal_spectrum(J2f, box, R)
            rel = spec[2:end] ./ max(spec[1], 1e-30)
            mmax = argmax(rel)
            push!(rows, join((round(sim.t, digits = 3), ekin, emag,
                              m[3], sqrt(sum(abs2, T)), T[3],
                              mmax, rel[mmax]), ","))
            emit_frame(Bf, Wf, img3d)
            tnext += FRAME_DT
        end
    end
    finish("done")
catch e
    if e isa InterruptException
        finish("interrupted")
    else
        rethrow()
    end
end
println("Outputs in $(abspath(OUT)); resume with: ",
        "julia -t auto --project=. scripts/v2_pathB_selfassembly.jl ",
        "$SCEN $NGRID <t_end> resume")
