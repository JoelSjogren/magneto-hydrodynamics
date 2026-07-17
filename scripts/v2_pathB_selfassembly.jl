# v2 path B — self-assembly from two tori (README §7.4). Evolve a pair of
# rings and watch for spontaneous structure: anapole growth |T|(t) and
# azimuthal breakup mode counts ("ring of rings").
#
#   julia -t auto --project=. scripts/v2_pathB_selfassembly.jl [scenario] [ngrid] [t_end]
#
# scenarios: counterhel — co-current, counter-helicity magnetic ring pair
#                          (lab analogue: spheromak merging → FRC)
#            opposed    — anti-parallel ring currents forced together
#            limnickels — counter-rotating vortex-ring collision with a
#                          weak frozen-in seed field
# Outputs: out/v2/<scenario>_N<ngrid>/

using FractalToroid
using Printf

const SCEN = length(ARGS) >= 1 ? ARGS[1] : "counterhel"
const NGRID = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 64
const T_END = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 25.0
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

mkpath(OUT)
rm(joinpath(OUT, "frames"); force = true, recursive = true)
mkpath(joinpath(OUT, "frames"))

box = Box(NGRID, HALF)
sim = MHDSim(box; cs = CS, eta = ETA, sponge_width = SPONGE)

if SCEN == "counterhel"
    add_flux_ring!(sim; R, a = A, z0 = -D, A0 = 0.6, Bt0 = 0.4)
    add_flux_ring!(sim; R, a = A, z0 = +D, A0 = 0.6, Bt0 = -0.4)
    add_vortex_ring!(sim; R, a = A, z0 = -D, P0 = 0.10)   # gentle approach
    add_vortex_ring!(sim; R, a = A, z0 = +D, P0 = -0.10)
elseif SCEN == "opposed"
    add_flux_ring!(sim; R, a = A, z0 = -D, A0 = 0.6, Bt0 = 0.0)
    add_flux_ring!(sim; R, a = A, z0 = +D, A0 = -0.6, Bt0 = 0.0)
    add_vortex_ring!(sim; R, a = A, z0 = -D, P0 = 0.30)   # forced together
    add_vortex_ring!(sim; R, a = A, z0 = +D, P0 = -0.30)
elseif SCEN == "limnickels"
    add_vortex_ring!(sim; R, a = A, z0 = -0.9, P0 = 0.40)
    add_vortex_ring!(sim; R, a = A, z0 = +0.9, P0 = -0.40)
    add_flux_ring!(sim; R, a = A, z0 = -0.9, A0 = 0.10, Bt0 = 0.0)
    add_flux_ring!(sim; R, a = A, z0 = +0.9, A0 = 0.10, Bt0 = 0.0)
else
    error("unknown scenario $SCEN")
end

# Symmetry-breaking noise: the initial conditions are axisymmetric, and on a
# clean grid the azimuthal instabilities behind ring breakup (Lim–Nickels)
# have nothing physical to grow from — without this, the mode spectrum only
# reads the Cartesian grid's m=4 anisotropy. Deterministic seed.
using Random
let rng = Xoshiro(1234)
    amp = 0.02 * maximum(abs, sim.S[MMZ])
    for q in (MMX, MMY, MMZ)
        A = sim.S[q]
        for idx in eachindex(A)
            A[idx] += amp * (2rand(rng) - 1)
        end
    end
end

# fixed render scales from the initial state
Bmag() = sqrt.(sim.S[MBX] .^ 2 .+ sim.S[MBY] .^ 2 .+ sim.S[MBZ] .^ 2)
function omag()
    ωx, ωy, ωz = curl_central(sim.S[MMX], sim.S[MMY], sim.S[MMZ],
                              box, sim.ip, sim.im)
    sqrt.(ωx .^ 2 .+ ωy .^ 2 .+ ωz .^ 2)
end
const LOGB_HI = log10(max(maximum(Bmag()), 1e-12))
const OM_HI = max(maximum(omag()), 1e-12) * 1.5

"4-panel frame: rows = log|B|, |ω|; columns = xz slice, xy (z=0) slice."
function render_frame(nframe)
    B = Bmag(); W = omag()
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
    save_png(joinpath(OUT, "frames", "frame_$(lpad(nframe, 5, '0')).png"), rgb)
end

"3D volume render: opacity = |B| (structure), color = |ω| (flow state)."
function render_frame3d(nframe)
    img = volume_render(omag(), Bmag(), box; res = 448, chi = OM_HI,
                        azim = 0.6, elev = 0.45)
    save_png(joinpath(OUT, "frames3d", "frame_$(lpad(nframe, 5, '0')).png"),
             img)
end
rm(joinpath(OUT, "frames3d"); force = true, recursive = true)
mkpath(joinpath(OUT, "frames3d"))

println("v2 path B [$SCEN] $(NGRID)³, t_end=$T_END, S=$(1/ETA), ",
        "$(Threads.nthreads()) threads")
rows = ["t,E_kin,E_mag,mz,Tnorm,Tz,mode_max,mode_amp"]
nframe = 0
render_frame(nframe); render_frame3d(nframe); nframe += 1
tnext = FRAME_DT
m0, T0 = grid_moments(sim)
@printf("t=0: E_kin=%.4g E_mag=%.4g |T|=%.4g\n",
        mhd_kinetic_energy(sim), mhd_magnetic_energy(sim),
        sqrt(sum(abs2, T0)))
while sim.t < T_END
    global nframe, tnext
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
        render_frame(nframe); render_frame3d(nframe); nframe += 1
        tnext += FRAME_DT
    end
end
write(joinpath(OUT, "timeseries.csv"), join(rows, "\n") * "\n")
m1, T1 = grid_moments(sim)
@printf("t=%.1f: E_kin=%.4g E_mag=%.4g |T|=%.4g (was %.4g)\n",
        sim.t, mhd_kinetic_energy(sim), mhd_magnetic_energy(sim),
        sqrt(sum(abs2, T1)), sqrt(sum(abs2, T0)))
println("Wrote $(nframe) frames and timeseries to $(abspath(OUT))")
