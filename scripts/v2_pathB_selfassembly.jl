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
#            random     — band-limited random flow + weak random seed field
#            bubble     — ULTR: imploding low-density cavity with a seeded
#                          re-entrant jet (cavitation collapse → vortex ring)
#            bubble2    — ULTR: two cavities imploding side by side, jets
#                          from mutual shielding only (no seeded jet)
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
# half=X sets the domain half-width (box spans [-X,X]³; default 2.0). Larger X
# at the same ngrid = bigger domain, coarser resolution; to hold resolution
# fixed, scale ngrid with X. Non-default writes to a "_halfX" dir.
const HALF = let i = findfirst(a -> startswith(a, "half="), ARGS)
    i === nothing ? 2.0 : parse(Float64, ARGS[i][6:end])
end
const POSARGS = filter(a -> !(a in ("gpu", "gpu32")) &&
                       !startswith(a, "seed=") && !startswith(a, "half="), ARGS)
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
                     "$(SCEN)_N$(NGRID)" *
                     (HALF == 2.0 ? "" : "_half$(HALF)") *
                     (SEED == 1234 ? "" : "_seed$(SEED)"))

box = Box(NGRID, HALF)
sim = MHDSim(box; cs = CS, eta = ETA, sponge_width = SPONGE)

"Scalar field = Σ of `nmodes` random cosine modes with integer wavevectors
up to `kmax` — a band-limited random field (energy at large scales, not grid
noise). The curl of three of these is a divergence-free random field."
function _rand_lowk(box::Box, rng, nmodes::Int, kmax::Int)
    n = box.n
    k0 = 2π / (n * box.dx)
    S = zeros(n, n, n)
    for _ in 1:nmodes
        kx = k0 * rand(rng, -kmax:kmax)
        ky = k0 * rand(rng, -kmax:kmax)
        kz = k0 * rand(rng, -kmax:kmax)
        φ = 2π * rand(rng)
        a = 2rand(rng) - 1
        @inbounds for k in 1:n
            z = center(box, k)
            for j in 1:n
                y = center(box, j)
                for i in 1:n
                    S[i, j, k] += a * cos(kx * center(box, i) + ky * y + kz * z + φ)
                end
            end
        end
    end
    S
end

"Carve a low-density spherical cavity into ρ: ×(ρc + (1−ρc)·step(r−Rb)).
Under the isothermal EOS p = cs²ρ the ambient pressure implodes it — the
MHD stand-in for a cavitation bubble at collapse onset."
function carve_cavity!(sim; Rb = 0.8, z0 = 0.0, rho_cav = 0.05, w = 0.08)
    ρ = sim.S[MRHO]
    s2 = w * sqrt(2.0)
    for k in 1:box.n, j in 1:box.n, i in 1:box.n
        r = sqrt(center(box, i)^2 + center(box, j)^2 +
                 (center(box, k) - z0)^2)
        f = 0.5 * (1 + FractalToroid._erf((r - Rb) / s2))
        ρ[i, j, k] *= rho_cav + (1 - rho_cav) * f
    end
end

"Gaussian momentum blob ρv_z += ρ·vz·exp(−|r−(0,0,z0)|²/2a²): the seeded
re-entrant jet. (A dense 'foil' slab can't provide the asymmetry here —
with p = cs²ρ it is high-pressure and explodes — so the jet that a nearby
boundary would produce is imposed directly.)"
function add_jet!(sim; z0 = 0.9, vz = -0.8, a = 0.35)
    ρ = sim.S[MRHO]
    inv2a2 = 1 / (2a^2)
    for k in 1:box.n, j in 1:box.n, i in 1:box.n
        d2 = center(box, i)^2 + center(box, j)^2 + (center(box, k) - z0)^2
        sim.S[MMZ][i, j, k] += ρ[i, j, k] * vz * exp(-d2 * inv2a2)
    end
end

"Rescale the (freshly seeded, previously zero) B to total energy `emag`.
A box-filling random seed would start |T| ≈ 4×10⁻² — 100× the self-assembly
signal — so bubble scenarios seed axisymmetric rims-threading flux rings
(|T|(0) = machine zero, like every other scenario) and scale them here."
function scale_seed_field!(sim; emag = 0.02)
    g = sqrt(emag / max(mhd_magnetic_energy(sim), 1e-30))
    sim.S[MBX] .*= g; sim.S[MBY] .*= g; sim.S[MBZ] .*= g
end

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
    elseif SCEN == "random"
        # Haphazard but energetic ICs: a divergence-free random velocity
        # field and a weak divergence-free random seed field, each the curl
        # of a band-limited (k ≤ 3) random vector potential — energy at large
        # scales with something to evolve, not instant grid-scale
        # dissipation. Same strong-flow/weak-field regime as limnickels, for
        # comparison; the test is whether coherent rings / anapole structure
        # crystallize out of disorder. Reproducible via seed=N.
        let rng = Xoshiro(SEED), km = 3, nm = 12
            vx, vy, vz = curl_central(_rand_lowk(box, rng, nm, km),
                                      _rand_lowk(box, rng, nm, km),
                                      _rand_lowk(box, rng, nm, km),
                                      box, sim.ip, sim.im)
            ρ0 = sim.S[MRHO]
            sim.S[MMX] .= ρ0 .* vx
            sim.S[MMY] .= ρ0 .* vy
            sim.S[MMZ] .= ρ0 .* vz
            f = sqrt(2.5 / max(mhd_kinetic_energy(sim), 1e-30))   # energetic
            sim.S[MMX] .*= f; sim.S[MMY] .*= f; sim.S[MMZ] .*= f
            bx, by, bz = curl_central(_rand_lowk(box, rng, nm, km),
                                      _rand_lowk(box, rng, nm, km),
                                      _rand_lowk(box, rng, nm, km),
                                      box, sim.ip, sim.im)
            sim.S[MBX] .= bx; sim.S[MBY] .= by; sim.S[MBZ] .= bz
            g = sqrt(0.05 / max(mhd_magnetic_energy(sim), 1e-30))  # weak seed
            sim.S[MBX] .*= g; sim.S[MBY] .*= g; sim.S[MBZ] .*= g
        end
    elseif SCEN == "bubble"
        # ULTR ultrasound-cleaner sim, simplified-MHD form: cavitation
        # collapse → re-entrant jet → toroidal vortex ring. Run in a larger
        # domain (half=3): the collapse emits shocks and the ring expands.
        carve_cavity!(sim; Rb = 0.8, rho_cav = 0.05, w = 0.08)
        add_jet!(sim; z0 = 0.9, vz = -0.8, a = 0.35)
        add_flux_ring!(sim; R = 0.8, a = 0.2, z0 = 0.0, A0 = 0.1, Bt0 = 0.0)
        scale_seed_field!(sim)
    elseif SCEN == "bubble2"
        # two cavities imploding along z: collapse asymmetry from mutual
        # shielding alone — no hand-seeded jet
        carve_cavity!(sim; Rb = 0.7, z0 = -0.9, rho_cav = 0.05, w = 0.08)
        carve_cavity!(sim; Rb = 0.7, z0 = +0.9, rho_cav = 0.05, w = 0.08)
        add_flux_ring!(sim; R = 0.7, a = 0.2, z0 = -0.9, A0 = 0.1, Bt0 = 0.0)
        add_flux_ring!(sim; R = 0.7, a = 0.2, z0 = +0.9, A0 = 0.1, Bt0 = 0.0)
        scale_seed_field!(sim)
    else
        error("unknown scenario $SCEN")
    end
    # symmetry-breaking noise (axisymmetric ICs on a clean grid leave
    # azimuthal instabilities nothing physical to grow from); applied only
    # at t = 0, so resumed runs continue the same trajectory
    let rng = Xoshiro(SEED)
        amp = 0.02 * maximum(abs, sim.S[MMZ])
        amp == 0 && (amp = 0.02 * CS)   # kick-free ICs (bubble2)
        for q in (MMX, MMY, MMZ)
            Aq = sim.S[q]
            for idx in eachindex(Aq)
                Aq[idx] += amp * (2rand(rng) - 1)
            end
        end
    end
end

const GSIM = GPU ? MHDCuda.to_gpu(sim; T = GPU32 ? Float32 : Float64) : nothing

# bubble scenarios get a third frame row showing ρ (the cavity/collapse is
# a density story; |B| and |ω| don't show it). Fixed linear scale: ambient
# ρ = 1 mid-scale, collapse compressions saturate.
const RHO_PANELS = startswith(SCEN, "bubble")
const RHO_HI = 2.5

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
       ["t,E_kin,E_mag,mz,Tnorm,Tz,mode_max,mode_amp,Jsq,Jmax"]

"Panel frame: rows = log|B|, |ω| (+ ρ if given); columns = xz, xy slices."
function render_frame(nf, B, W, R)
    j0 = box.n ÷ 2
    npx = box.n * UP
    nrow = R === nothing ? 2 : 3
    rgb = fill(0x12, 3, 2npx + GAP, nrow * npx + (nrow - 1) * GAP)
    panels = Any[(B, :logb, :xz, 0, 0), (B, :logb, :xy, 1, 0),
                 (W, :om, :xz, 0, 1), (W, :om, :xy, 1, 1)]
    R === nothing ||
        append!(panels, ((R, :rho, :xz, 0, 2), (R, :rho, :xy, 1, 2)))
    for (F, scale, sl, px, py) in panels
        for b in 1:box.n, a in 1:box.n
            v = sl === :xz ? F[a, j0, b] : F[a, b, j0]
            vv = scale === :logb ?
                 clamp((log10(max(v, 1e-300)) - (LOGB_HI - 3)) / 3, 0, 1) :
                 scale === :om ? clamp(v / OM_HI, 0, 1) :
                 clamp(v / RHO_HI, 0, 1)
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

function emit_frame(B, W, R, img3d)
    global nframe
    nf = nframe
    if img3d === nothing
        render_frame(nf, B, W, R)
        render_frame3d(nf, volume_render(W, B, box; res = 448, chi = OM_HI,
                                         azim = 0.6, elev = 0.45))
        dump_volumes(nf, B, W)
    else
        flush_writer()
        WRITER[] = Threads.@spawn begin
            render_frame(nf, B, W, R)
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
    emit_frame(_bmag(), _omag(), RHO_PANELS ? sim.S[MRHO] : nothing, nothing)
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
                Rf = RHO_PANELS ? Float64.(Array(GSIM.S[MRHO])) : nothing
            else
                Bf = _bmag(); Wf = _omag()
                Jx, Jy, Jz = curl_central(sim.S[MBX], sim.S[MBY], sim.S[MBZ],
                                          box, sim.ip, sim.im)
                J2f = Jx .^ 2 .+ Jy .^ 2 .+ Jz .^ 2
                m, T = grid_moments(sim)
                ekin = mhd_kinetic_energy(sim)
                emag = mhd_magnetic_energy(sim)
                img3d = nothing
                Rf = RHO_PANELS ? sim.S[MRHO] : nothing
            end
            spec = azimuthal_spectrum(J2f, box, R)
            rel = spec[2:end] ./ max(spec[1], 1e-30)
            mmax = argmax(rel)
            # current J = ∇×B: total ∫|J|²dV (activity/dissipation) and peak
            # |J| (current-sheet strength) — separates a uniformly decaying
            # current from one that reorganizes toward the anapole
            jsq = sum(J2f) * cellvol(box)
            jmax = sqrt(maximum(J2f))
            push!(rows, join((round(sim.t, digits = 3), ekin, emag,
                              m[3], sqrt(sum(abs2, T)), T[3],
                              mmax, rel[mmax], jsq, jmax), ","))
            emit_frame(Bf, Wf, Rf, img3d)
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
