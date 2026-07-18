# v2 path B — self-assembly from two tori (README §7.4). Evolve a pair of
# rings and watch for spontaneous structure: anapole growth |T|(t) and
# azimuthal breakup mode counts ("ring of rings").
#
#   julia -t auto --project=. scripts/v2_pathB_selfassembly.jl [scenario] [ngrid] [t_end] [resume]
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

const SCEN = length(ARGS) >= 1 ? ARGS[1] : "counterhel"
const NGRID = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 64
const T_END = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 25.0
const RESUME = length(ARGS) >= 4 && ARGS[4] == "resume"
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
const OUT = joinpath(@__DIR__, "..", "out", "v2", "$(SCEN)_N$(NGRID)")

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
    let rng = Xoshiro(1234)
        amp = 0.02 * maximum(abs, sim.S[MMZ])
        for q in (MMX, MMY, MMZ)
            Aq = sim.S[q]
            for idx in eachindex(Aq)
                Aq[idx] += amp * (2rand(rng) - 1)
            end
        end
    end
end

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
function render_frame(nf)
    B = _bmag(); W = _omag()
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
function render_frame3d(nf)
    img = volume_render(_omag(), _bmag(), box; res = 448, chi = OM_HI,
                        azim = 0.6, elev = 0.45)
    save_png(joinpath(OUT, "frames3d", "frame_$(lpad(nf, 5, '0')).png"), img)
end

"uint8 volume dumps for the interactive raycaster (gitignored)."
function dump_volumes(nf)
    tag = lpad(nf, 5, '0')
    B = _bmag(); W = _omag()
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

function emit_frame()
    global nframe
    render_frame(nframe)
    render_frame3d(nframe)
    dump_volumes(nframe)
    nframe += 1
end

"Write timeseries + checkpoint (normal completion and SIGINT both land here)."
function finish(status)
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
    emit_frame()
end
write(joinpath(@__DIR__, "..", "out", "v2", "CURRENT"), "$(SCEN)_N$(NGRID)")

println("v2 path B [$SCEN] $(NGRID)³, t: $(round(sim.t, digits=2)) → $T_END, ",
        "S=$(1/ETA), $(Threads.nthreads()) threads")
try
    while sim.t < T_END
        global tnext
        mhd_step!(sim)
        if sim.t >= tnext
            Jx, Jy, Jz = curl_central(sim.S[MBX], sim.S[MBY], sim.S[MBZ],
                                      box, sim.ip, sim.im)
            J2 = Jx .^ 2 .+ Jy .^ 2 .+ Jz .^ 2
            spec = azimuthal_spectrum(J2, box, R)
            rel = spec[2:end] ./ max(spec[1], 1e-30)
            mmax = argmax(rel)
            m, T = grid_moments(sim)
            push!(rows, join((round(sim.t, digits = 3),
                              mhd_kinetic_energy(sim), mhd_magnetic_energy(sim),
                              m[3], sqrt(sum(abs2, T)), T[3],
                              mmax, rel[mmax]), ","))
            emit_frame()
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
