# Fractal Toroidal Moment — Continuum Field Simulation

> **Note:** this project was vibe coded — researched, designed, implemented,
> and run by an AI agent (Claude) under human direction. Read the physics and
> the code with the corresponding skepticism.

An attempt to reproduce, by numerical continuum field simulation, the "fractal
toroidal moment" (FTM) advocated by Bob Greenyer of the Martin Fleischmann
Memorial Project (MFMP).

## 1. Project brief

Goal: find out whether classical continuum field theory — Maxwell's equations
coupled to some definite specification of moving matter — can produce and
sustain the fractal-toroid, monopole-like magnetic structure that Greenyer
describes, and if so under what conditions.

Points from the initial brief:

- Greenyer sometimes attributes the phenomenon to the movement of **"relic
  neutrinos"**. Neutrinos are electrically (almost) neutral, so they cannot
  source the current term in Maxwell's equations. For a simulatable model the
  working assumption is that the current is carried by **massive charged
  matter** (electrons and/or ions), accelerated self-consistently by the
  electromagnetic field it generates.
- The "fractal toroid" is a **coil of coils** — "a wheel within a wheel within
  a wheel", as Greenyer puts it: a current path wound helically around a
  torus, that helix itself wound from a finer helix, recursively.
- The effect is claimed to **survive for ~2 days** after the driving apparatus
  is removed. If true, the persistent object cannot depend on the details of
  its environment, so a **single phase of matter in vacuum** may be a
  sufficient model — and persistence time becomes the primary measurable
  output of the simulation.
- Implementation language: C++, Rust, Julia, or any other high-performance
  language. GPU execution (e.g. CUDA) is especially attractive. Target
  production hardware: an **NVIDIA GTX 1080 (8 GB)** on a separate machine.

## 2. Research notes: what Greenyer actually claims

Bob Greenyer is a volunteer researcher with the MFMP (LENR / cold-fusion
replication project, active since 2012). His main written outlet is the
**Remote View** Substack; the most on-point items found:

- **"The Fractal Toroidal Moment"** (Remote View; also presented at Cosmic
  Summit 2024, slides downloadable from the post):
  <https://remoteview.substack.com/p/the-fractal-toroidal-moment-7b5>
- **"Practical Applications of the Fractal Toroidal Moment"** — ICCF-25
  conference contribution (Szczecin, August 2023). Abstract/post:
  <https://remoteview.substack.com/p/practical-applications-of-the-fractal>,
  announcement: <https://e-catworld.com/2023/04/17/abstract-practical-applications-of-the-fractal-toroidal-moment-bob-greenyer/>,
  talk video: <https://www.youtube.com/watch?v=lnOeQdUaslE>
- **"Decoding EVOs: A Deep Dive into Exotic Vacuum Objects"** — Alternative
  Propulsion Engineering Conference writeup of his framework:
  <https://www.altpropulsion.com/decoding-evos-a-deep-dive-into-exotic-vacuum-objects/>
- **"O-Day"** posts (ball-lightning "O" structures in water/plasma):
  <https://remoteview.substack.com/p/o-day>
- Michael Clarage (SAFIRE), **"Fractal Toroids — Part 1, Geometry"** — the
  most quantitative description of the geometry in Greenyer's orbit:
  <https://michaelclarage.substack.com/p/fractal-toroids-part-1-geometry>

### 2.1 The claim set, distilled

1. **Morphology.** Across LENR experiments, ball-lightning traces, discharge
   residues and "strange radiation" tracks, Greenyer identifies a recurring
   signature: self-similar nested tori ("shell forms of tori, spindle tori,
   spheres and their aggregates") — the FTM. He treats it as the same object
   as Winston Bostick's *plasmoids* and Ken Shoulders' *EVOs / charge
   clusters*.
2. **Geometry** (per Clarage): a torus whose tube is itself a ring of smaller
   tori, recursively. Working numbers: tube-to-hole aspect ratio ~1:4, ~48
   sub-tori per level, Hausdorff dimension log 48 / log 4 ≈ 2.8. Claim: with
   **≥ 3 levels of nesting the E and B fields become self-contained** and do
   not propagate to the outside world (this is essentially the *anapole /
   non-radiating configuration* idea).
3. **Toroidal (anapole) moment.** Greenyer explicitly anchors on the toroidal
   dipole moment of Zel'dovich (1957) and Dubovik & Tugushev (Phys. Rep. 187,
   1990), and cites Papasimakis et al., *"Electromagnetic toroidal excitations
   in matter and free space"*, Nature Materials 15, 263 (2016), plus Nemkov et
   al. on non-radiating sources. A poloidal current on a torus has zero charge,
   zero electric/magnetic dipole — its leading moment is the toroidal moment
   **T = (1/10) ∫ [(r·J) r − 2 r² J] d³r**.
4. **Monopole-like behavior.** Via Ken Shoulders he calls the EVO an "ideal
   monopole oscillator"; claims of Möbius-strip-like winding producing
   topological monopole signatures. (Note: ∇·B = 0 forbids a true monopole
   moment; at most a configuration can *locally mimic* a monopole's field over
   some shell, or radiate anapole-type potentials. Quantifying what a
   magnetometer would actually read outside such an object is something this
   project can do.)
5. **Charge separation and carriers.** The process is "driven in part by true
   charge separation"; Shoulders' EVs were micron-scale clusters of ~10⁸–10¹¹
   electrons with a small ion fraction. Nowhere is a definite current carrier
   specification given — this is the gap our model has to fill.
6. **Relic neutrinos.** Greenyer hypothesizes the EVO axis carries "a coherent
   beam of spinning dark matter, likely relic neutrinos", and that ordinary
   matter "breathes energy from the relic neutrino background". This is not
   simulatable in classical field theory and is **out of scope**; we model the
   charged-matter sector only and treat any such background at most as a
   phenomenological energy source term (off by default).
7. **Persistence.** Ball-lightning-like objects surviving detachment from the
   apparatus; the ~2-day figure is from Greenyer's talks (exact written source
   not yet pinned down — TODO: locate it in the Remote View archive).

### 2.2 Mainstream anchors (real physics to build on)

The FTM narrative overlaps a set of perfectly respectable topics:

- **Toroidal multipoles / anapoles**: Zel'dovich 1957; Dubovik & Tugushev
  1990; Papasimakis et al. 2016 (measured in metamaterials and nuclei; Cs-133
  has a measured nuclear anapole moment).
- **Force-free / Beltrami fields**: ∇×B = λB; Chandrasekhar–Kendall states;
  Woltjer's theorem (minimum energy at fixed magnetic helicity); Taylor
  relaxation in laboratory plasmas; spheromaks and compact toroids.
- **MHD topological solitons**: Kamchatnov, *"Topological solitons in
  magnetohydrodynamics"*, Sov. Phys. JETP 55, 69 (1982) — an exact Hopf-fibered
  force-free field.
- **Electromagnetic knots / hopfions**: Rañada 1989; Rañada & Trueba, *"Ball
  lightning an electromagnetic knot?"*, Nature 383, 32 (1996); Rañada, Soler &
  Trueba, *"Ball lightning as force-free magnetic knots"*, Phys. Rev. E 62,
  7181 (2000).
- **Plasmoids**: Bostick, Phys. Rev. 106, 404 (1957).

What is *not* mainstream: monopole moments, days-long persistence of an
isolated plasma/EM object in vacuum, relic-neutrino coupling, transmutation.
The known obstacle to persistence is the **virial theorem (Shafranov)**: a
finite blob of classical charged matter + EM field in vacuum, with nothing
pushing in from outside, has no static equilibrium — it must expand, radiate,
or fly apart. (The textbook derivation integrates the Maxwell stress tensor
— field-side bookkeeping — but the conclusion survives charge-side
reformulation: the virial of the interparticle Coulomb/Ampère forces gives
the same no-equilibrium result, so the obstacle does not hinge on where one
locates the energy; cf. neutrinos.md §3.) Any persistence in our model must therefore be *dynamic*
(inertially/rotationally supported, or a slowly-decaying self-organized
state). Measuring the decay time — and whether fractal nesting extends it — is
the central scientific question of this project, and it is falsifiable both
ways.

## 3. Proposed minimal continuum model

### 3.1 Modeling decisions

- **Carrier**: a single cold charged fluid of electrons (charge −e, mass mₑ),
  optionally with a static neutralizing ion background of adjustable fraction
  f ∈ [0, 1] (f = 0: pure electron cluster à la Shoulders; f = 1:
  quasineutral). A second mobile ion fluid is a later extension.
- **Environment**: vacuum, open (absorbing) boundaries. Justified by the
  survives-removal claim.
- **Closure**: pressureless ("cold") to start. Cold fluids wave-break
  (density caustics); if that dominates, add a small isothermal pressure term
  or switch to particle-in-cell.
- **Relativity**: optional γ-factor on the momentum equation; Shoulders-type
  electron clusters are plausibly relativistic, but v ≪ c keeps v1 simple.

### 3.2 Equations (SI)

Maxwell:

$$
\partial_t \mathbf{B} = -\nabla\times\mathbf{E}, \qquad
\partial_t \mathbf{E} = c^2\,\nabla\times\mathbf{B} - \mathbf{J}/\varepsilon_0
$$

$$
\nabla\cdot\mathbf{E} = \rho_q/\varepsilon_0, \qquad
\nabla\cdot\mathbf{B} = 0
$$

Cold electron fluid (number density n, velocity **u**):

$$
\partial_t n + \nabla\cdot(n\mathbf{u}) = 0
$$

$$
\partial_t \mathbf{u} + (\mathbf{u}\cdot\nabla)\mathbf{u}
  = -\frac{e}{m_e}\left(\mathbf{E} + \mathbf{u}\times\mathbf{B}\right)
$$

Coupling:

$$
\rho_q = e\,(f\,n_b - n), \qquad \mathbf{J} = -e\,n\,\mathbf{u}
$$

with $n_b(\mathbf{x})$ the frozen ion background (equal to the initial $n$,
scaled by f). This is the **Euler–Maxwell (cold plasma) system** — arguably
the simplest self-consistent instantiation of "Maxwell's equations where the
current is massive charge accelerated by the EM field".

### 3.3 Normalization

Normalize to plasma units: $\omega_p^2 = n_0 e^2/(\varepsilon_0 m_e)$, skin
depth $d_e = c/\omega_p$, with $t \to \omega_p t$, $\mathbf{x} \to \mathbf{x}/d_e$,
$\mathbf{u} \to \mathbf{u}/c$, $\mathbf{E} \to e\mathbf{E}/(m_e c\,\omega_p)$,
$\mathbf{B} \to e\mathbf{B}/(m_e\omega_p)$,
$n \to n/n_0$. The dimensionless system has **no free constants**; all physics
is in the geometry, the initial speed $u_0 = v_0/c$, and the neutralization
fraction f. Results rescale to any physical size afterwards.

### 3.4 Initial condition: the fractal coil

Recursive space curve ("coil of coils"), built with a parallel-transport frame
(not Frenet, which degenerates on straight segments):

- Level 0: circle of radius $R_0$:
  $\mathbf{C}_0(\theta) = R_0(\cos\theta, \sin\theta, 0)$.
- Level k: wind $w_k$ turns around the level-(k−1) tube at radius $a_k$:

$$
\mathbf{C}_k(\theta) = \mathbf{C}_{k-1}(\theta)
 + a_k\left[\mathbf{N}_{k-1}(\theta)\cos(w_k\theta+\varphi_k)
 + \mathbf{B}\!\mathbf{N}_{k-1}(\theta)\sin(w_k\theta+\varphi_k)\right]
$$

where $\mathbf{N}, \mathbf{BN}$ are the transported normal/binormal of the
previous level. Clarage's numbers suggest $a_k/a_{k-1} \approx 1/4$,
$w_k \approx 48$; both are parameters. The deepest curve is smeared into a
Gaussian tube of radius σ to define the initial $n(\mathbf{x})$ and
$\mathbf{u}(\mathbf{x})$ (hence $\mathbf{J}$), and the initial **E**, **B**
are solved from the static (Biot–Savart / Poisson) fields of that
distribution so the run starts near self-consistency.

### 3.5 Diagnostics

- Energies: field $\tfrac{\varepsilon_0}{2}\!\int (E^2 + c^2B^2)$, kinetic
  $\tfrac{m_e}{2}\!\int n u^2$; Poynting flux through the boundary
  (radiated/lost power). *(Accounting note, cf. neutrinos.md §3: the
  field-energy density is one bookkeeping convention; the charge-side
  convention ½∫(ρφ + J·A) gives identical totals wherever sources exist.
  Our τ metrics localize energy to a core box — a "where" statement that
  is convention-dependent in principle, though numerically the two agree
  closely here since B concentrates at the currents. Notably, our
  radiated-power measurement — energy tallied when the sponge absorbs it
  — is already absorber-side accounting in the Wheeler–Feynman sense.)*
- **Magnetic helicity** $H = \int \mathbf{A}\cdot\mathbf{B}\,dV$ — the
  conserved quantity that drives Taylor relaxation toward force-free states.
- Moments: magnetic dipole $\mathbf{m} = \tfrac12\int \mathbf{r}\times\mathbf{J}\,dV$,
  toroidal/anapole moment
  $\mathbf{T} = \tfrac1{10}\int [(\mathbf{r}\cdot\mathbf{J})\mathbf{r} - 2r^2\mathbf{J}]\,dV$,
  and their ratio vs. nesting depth k.
- **Persistence time** τ(k): e-folding time of the confined field/kinetic
  energy vs. number of fractal levels k. Tests the "3 levels self-contain the
  field" claim directly.
- Synthetic magnetometer: $|\mathbf{B}|$ and its radial profile on shells
  around the object — does anything monopole-*looking* ($B_r \sim 1/r^2$ over
  a shell) appear?

## 4. Numerical plan

Discretization: Yee-grid FDTD for Maxwell (exactly preserves ∇·B = 0),
finite-volume/upwind for the fluid on the same grid, CPML absorbing
boundaries. Memory per cell ≈ 10 floats (E, B, n, u) → a 512³ FP32 grid is
~5.4 GB: **fits the GTX 1080's 8 GB**, with 384³ comfortable. Pascal (sm_61)
has 1:32 FP64 throughput, so we stay FP32 (with FP64 reductions for
diagnostics), and pin CUDA ≤ 12.x (CUDA 13 dropped Pascal).

Phases:

1. **Geometry + magnetostatics** (CPU): fractal-curve generator, Biot–Savart
   field of the frozen coil; reproduce the FTM field morphology; compute
   dipole vs. anapole content vs. depth k.
2. **Vacuum FDTD**: prescribed (non-self-consistent) current pulse through the
   fractal coil; measure how radiation escape varies with k.
3. **Self-consistent Euler–Maxwell**: release the fluid; watch relaxation
   (helicity-conserving?), measure τ(k), f-dependence, u₀-dependence.
4. **CUDA port** of phase 3 for the GTX 1080; parameter scans.

Language: **Julia** is the recommended starting point — phases 1–3 in plain
Julia, then the same kernels moved to GPU via CUDA.jl without a rewrite.
Fallback/alternative for maximum control: C++ with CUDA. (Rust is viable but
its GPU story adds friction for stencil codes.)

## 4a. Implementation status

Phases 1–3 are implemented as the Julia package `FractalToroid` in this
repository (stdlib-only, no external dependencies; multi-threaded). Layout:

- `src/geometry.jl` — recursive coil-of-coils curve (parallel-transport frame
  with holonomy correction so every level closes)
- `src/fields.jl` — Biot–Savart (cell-centered and Yee-staggered), vector
  potential, dipole + anapole moments, helicity, Gaussian-tube current
  deposition, shell-averaged |B| profiles
- `src/yee.jl` — Maxwell on a periodic Yee grid, SSP-RK3 method of lines,
  sponge absorbing layer with removed-energy bookkeeping
- `src/fluid.jl` — cold electron fluid (Rusanov finite volume) coupled to the
  Yee grid; Marder divergence cleaning keeps Gauss's law
- `src/png.jl` — dependency-free PNG writer for field slices
- `scripts/phase{1,2,3}_*.jl` — the three experiments; outputs in `out/`

Run everything (from the repo root):

```sh
julia --project=. test/runtests.jl                          # validation suite
julia -t auto --project=. scripts/phase1_geometry.jl        # [K] [ngrid]
julia -t auto --project=. scripts/phase2_fdtd.jl            # [K] [ngrid]
julia -t auto --project=. scripts/phase3_euler_maxwell.jl   # [K] [ngrid] [t_end]
```

Validation: Biot–Savart against the analytic loop field; m = IπR² and T = 0
for the bare loop; exact div B preservation and energy conservation in vacuum
FDTD; cold-plasma (Langmuir) oscillation at ω_p in the coupled fluid.

First results (quick-mode grids, windings 16/8 at ratio 1/4):

- **Phase 1**: nesting switches on the anapole moment and magnetic helicity
  (k=0: |T| = 0, H = 0 exactly; k=1: |T|/|m| ≈ 0.5, H > 0). The far field
  stays dipole-dominated (log-log slope ≈ −3): the fractal coil is *not*
  monopole-like in its static far field.
- **Phase 2**: at drive wavelength ~ coil size, deeper nesting radiates
  *more*, not less (longer wire at fixed peak current; the anapole
  non-radiating property is a quasi-static/point-source statement). After
  switch-off, the core field energy drops ~3 orders of magnitude within one
  light-crossing regardless of k — prescribed-current vacuum fields do not
  persist, as expected.
- **Phase 3** (48³, u₀ = 0.05c, t = 60/ω_p; `out/phase3/summary_*.csv`):
  - *uniform mode* (coil in a uniform quasineutral plasma): τ ≈ 51, 43, 41/ω_p
    for k = 0, 1, 2. Radiated energy is ≲0.1% — the medium traps radiation —
    and the torus fully phase-mixes away by t = 60.
  - *ball mode* (isolated quasineutral plasma torus in near-vacuum, the
    configuration the persistence claim is about): at 96³, τ ≈ 276, 125,
    74/ω_p for k = 0, 1, 2 (k = 0 agrees with the 48³ value 274 — well
    converged), with only a few % of the energy radiated; the tube structure
    is still recognizable at t = 60, with the density hollowing into a shell
    around each tube axis. An isolated quasineutral current ring is thus
    *quasi-stable* on ~10² plasma periods — but **each level of fractal
    nesting roughly halves the lifetime**, the opposite of the
    self-containment claim. For scale: 2-day persistence at solid-state-like
    densities would need τ ~ 10¹⁵–10²⁰/ω_p.

Known numerical caveats: SSP-RK3 weakly damps grid-scale modes (keep
structures ≳ 3 cells); the Gaussian splat's ∇·J ≠ 0 residue is handled by
Marder cleaning in phase 3 but uncorrected in phase 2; quick-mode grids are
coarse — treat trends, not absolute numbers, as meaningful.

### Videos

Phase-3 ball-mode relaxation at 96³ (t = 0 → 60/ω_p, 24 fps ≈ 10 s each;
left panel: log₁₀|B| on a fixed 3.5-decade scale, right panel: electron
density; xz-slice through the torus). Rendered inline during the simulation,
stitched with `scripts/make_videos.sh`:

Full-quality mp4s: [k = 0, bare current ring — τ ≈ 276/ω_p](out/videos/phase3_k0_ball.mp4)
· [k = 1, coil (16 turns) — τ ≈ 125/ω_p](out/videos/phase3_k1_ball.mp4)
· [k = 2, coil of coils (16×8) — τ ≈ 74/ω_p](out/videos/phase3_k2_ball.mp4)

**k = 0** (τ ≈ 276/ω_p):

![k=0 ball-mode relaxation](out/videos/phase3_k0_ball.gif)

**k = 1** (τ ≈ 125/ω_p):

![k=1 ball-mode relaxation](out/videos/phase3_k1_ball.gif)

**k = 2** (τ ≈ 74/ω_p):

![k=2 ball-mode relaxation](out/videos/phase3_k2_ball.gif)

## 5. Implementation provenance

What each component of the model and code is based on, step by step. (The
original project brief is preserved verbatim in `initial-prompt.txt`.)

**Model choice (cold Euler–Maxwell).** The "simplest possible" reading of the
brief: Maxwell + a single cold charged fluid is the standard *cold plasma*
model of textbook plasma physics (e.g. Nicholson, *Introduction to Plasma
Theory*; Krall & Trivelpiece). The electron-fluid-plus-static-ion-background
split follows Ken Shoulders' picture of EVs as electron clusters with a small
ion fraction. The plasma-unit normalization (ω_p, skin depth d_e) is the
standard nondimensionalization used throughout particle-in-cell literature.

**Fractal coil geometry.** The recursive helix-on-helix construction
formalizes Greenyer's "wheel within a wheel"; default aspect ratio 1/4 and the
idea of a fixed sub-coil count per level come from Michael Clarage's "Fractal
Toroids — Part 1, Geometry". The frame transported along the parent curve is
the *rotation-minimizing (Bishop) frame* from computer-aided geometric design
(Bishop 1975, "There is more than one way to frame a curve"; computed here by
tangent-projection, with the closure holonomy distributed as a uniform
compensating twist — the holonomy itself is the geometric phase familiar from
the Călugăreanu–White–Fuller twist/writhe decomposition).

**Magnetostatics.** Biot–Savart summation over segments with a softened
denominator is the standard *regularized filament* method from vortex-filament
and coil-design codes (the softening length plays the role of the tube/blob
radius). Dipole moment m = ½∮ r×I dl is textbook (Jackson §5.6); the toroidal
(anapole) moment formula T = (1/10)∫[(r·J)r − 2r²J]dV is from the toroidal
multipole literature (Dubovik & Tugushev 1990; Radescu & Vaman 2002;
Papasimakis et al. 2016). Magnetic helicity ∫A·B dV and its role as the
conserved driver of relaxation: Woltjer 1958, Taylor 1974.

**Maxwell solver.** Staggered-grid spatial discretization: Yee 1966 (the FDTD
scheme; see Taflove & Hagness, *Computational Electrodynamics*). Instead of
Yee's leapfrog we use method-of-lines with SSP-RK3 (Shu & Osher 1988) so the
fluid and field share one integrator; the staggering still gives exact
∇·B = 0. The absorbing boundary is a graded "sponge"/masking layer — the
pre-PML absorber family (Israeli & Orszag 1981 damping layers), chosen over
CPML for simplicity, with the absorbed energy tallied as the radiation
diagnostic.

**Fluid solver.** Pressureless cold advection is weakly hyperbolic and forms
δ-shocks (pressureless gas dynamics literature, e.g. Bouchut 1994), so the
flux uses the Rusanov / local Lax–Friedrichs scheme (Rusanov 1961) —
first-order dissipation regularizes the caustics. Gauss's law maintenance by
diffusing ∇·E − ρ errors: Marder 1987, refined by Langdon 1992 —
standard practice in electromagnetic PIC codes.

**Validation tests.** Loop field against the textbook on-axis formula; cold
Langmuir oscillation at ω_p as the canonical coupled-system test (any plasma
textbook); vacuum energy conservation and div-B preservation as the standard
FDTD sanity checks.

**Experiment design.** Phase 2's radiated-power-vs-nesting question
operationalizes the anapole "non-radiating configuration" claims (Zel'dovich
1957; Nemkov et al. on non-radiating sources). Phase 3's
isolated-ball-vs-uniform-medium contrast and the τ(k) persistence measurement
are framed by Shafranov's virial theorem (no static self-confined
equilibrium) and by the ball-lightning-as-relaxed-state proposals (Rañada &
Trueba 1996/2000; Kamchatnov's MHD hopfion 1982).

**Infrastructure.** The dependency-free PNG writer implements the PNG spec
(W3C/RFC 2083) with zlib "stored" deflate blocks (RFC 1950/1951); the
colormap is a linear interpolation of viridis anchors (Smith & van der Walt,
matplotlib). The README math preview uses marked.js + MathJax + GitHub's
markdown CSS; videos are stitched with ffmpeg (libx264, yuv420p for player
compatibility).

## 6. Source list

- Remote View (Greenyer's Substack): <https://remoteview.substack.com>
  - The Fractal Toroidal Moment: <https://remoteview.substack.com/p/the-fractal-toroidal-moment-7b5>
  - Practical Applications of the FTM (ICCF-25): <https://remoteview.substack.com/p/practical-applications-of-the-fractal>
  - O-Day: <https://remoteview.substack.com/p/o-day>
- ICCF-25 abstract announcement: <https://e-catworld.com/2023/04/17/abstract-practical-applications-of-the-fractal-toroidal-moment-bob-greenyer/>
- ICCF-25 talk video: <https://www.youtube.com/watch?v=lnOeQdUaslE>
- Decoding EVOs (altpropulsion.com): <https://www.altpropulsion.com/decoding-evos-a-deep-dive-into-exotic-vacuum-objects/>
- EVOs, Transmutation & Anomalies (altpropulsion.com): <https://www.altpropulsion.com/exotic-vacuum-objects-evos-transmutation-anomalies/>
- Michael Clarage, Fractal Toroids Part 1 — Geometry: <https://michaelclarage.substack.com/p/fractal-toroids-part-1-geometry>
- LENR-forum video thread (index of Greenyer's talks): <https://www.lenr-forum.com/forum/thread/7244-bob-greenyer-mfmp-video-thread/>
- Papasimakis, Fedotov, Savinov, Raybould & Zheludev, Nature Materials 15, 263 (2016).
- Dubovik & Tugushev, Phys. Rep. 187, 145 (1990).
- Kamchatnov, Sov. Phys. JETP 55, 69 (1982).
- Rañada & Trueba, Nature 383, 32 (1996); Rañada, Soler & Trueba, Phys. Rev. E 62, 7181 (2000).
- Bostick, Phys. Rev. 106, 404 (1957).
- Ken Shoulders, "EV — A Tale of Discovery" (1987); US Patent 5,018,180.

## 7. v2 plan: fluid MHD and self-assembly from counter-rotating tori

v1 established what the cold, collisionless, frozen-ion model does with a
*hand-built* fractal coil: anapole + helicity switch on, far field stays
dipolar, and every nesting level shortens the life of the object (τ ≈ 249,
127, 37/ω_p at 144³, k = 2 unconverged and falling). v2 changes two things,
both motivated by a closer reading of the claims.

### 7.1 Motivation

1. **The setting should be fluid, not near-solid.** The source vocabulary is
   magneto-*hydro*-dynamics and vortices, and days-long persistence after
   removal of the apparatus implies a self-carrying lump of matter, not
   charges electrostatically clamped to a frozen background as in v1 (whose
   ion lattice made it behave like a charged solid). v2 moves to a
   single-fluid quasineutral visco-resistive MHD model: the matter itself
   flows, vortices are first-class objects, and decay happens by
   reconnection and resistive diffusion instead of collisionless phase
   mixing.
2. **The fractal should emerge, not be imposed.** Verified against the
   sources (see 7.2): Greenyer's central formation claim is *self-assembly*
   from simple beginnings — a collapsing ionized bubble forming **two
   counter-rotating vortices** with a "zero point" equilibrium plane between
   them, and EVOs then linking "into filaments and larger self-similar
   rings", which "further cluster into rings of rings". v1 built the coil of
   coils by hand; v2 asks whether anything like it self-assembles from a
   two-torus initial condition.

### 7.2 What the sources actually claim about formation

- Two counter-rotating vortices from symmetric cavitation collapse, with a
  "zero point" plane between them (Thunderstorm-Generator orbit of MFMP:
  <https://alchemicalscience.org/thunderstorm-generator-q-a-more-proof-of-cavitation-from-mfmp-plasmoid-tech-updates/>).
- The EVO core as a toroidal "vortex/counter-vortex" configuration
  (<https://www.altpropulsion.com/decoding-evos-a-deep-dive-into-exotic-vacuum-objects/>).
- Hierarchical clustering: EVOs link along their axis into filaments and
  self-similar rings, then "rings of rings" ("THOR — Outside of the Inside":
  <https://remoteview.substack.com/p/thor-outside-of-the-inside>; the axial
  link is attributed to a relic-neutrino beam — out of classical scope, we
  test what plain MHD does).
- "Charge separation and turbulence is critical for formation", with
  acoustic/ion-acoustic resonance as an organizing factor (same THOR post).
  Charge separation is inaccessible to single-fluid MHD — noted as a v3
  direction (Hall or two-fluid MHD).

Mainstream anchors that make self-assembly a fair, testable question:
head-on collision of counter-rotating vortex rings breaks up, via azimuthal
instability, into a **ring of secondary small rings** (Lim & Nickels, Nature
357, 225 (1992)) — literally "ring → ring of rings" in an ordinary fluid;
laboratory merging of two spheromaks self-organizes into new objects, with
outcome controlled by relative helicity sign (counter-helicity → FRC:
Yamada, Ono et al., TS-3 experiments, PRL 65, 721 (1990); Ono et al. 1993);
kink instability converts twist into writhe, i.e. turns an over-twisted flux
ring into a helix — a generic mechanism for spontaneous "coiling of a coil";
Taylor relaxation constrains what any of this can settle into.

### 7.3 v2 model equations

Single-fluid, quasineutral, compressible, visco-resistive MHD in Alfvén
units (B₀, ρ₀, v_A = B₀/√(μ₀ρ₀); lengths in R₀, time in R₀/v_A):

$$
\partial_t \rho + \nabla\cdot(\rho\mathbf{v}) = 0
$$

$$
\partial_t(\rho\mathbf{v}) + \nabla\cdot\left[\rho\mathbf{v}\mathbf{v}
 + \left(p + \tfrac{B^2}{2}\right)\mathbf{I} - \mathbf{B}\mathbf{B}\right]
 = \mu\,\nabla^2\mathbf{v}
$$

$$
\partial_t\mathbf{B} = \nabla\times(\mathbf{v}\times\mathbf{B})
 + \eta\,\nabla^2\mathbf{B}, \qquad p = c_s^2\,\rho \;\text{(isothermal)}
$$

with hyperbolic (GLM/Dedner) divergence cleaning. Dimensionless knobs:
Lundquist number S = 1/η (resistive lifetime = S in code units), sound ratio
c_s/v_A (plasma β), optional explicit viscosity μ (Rusanov dissipation
otherwise), ambient density floor ρ_amb ≈ 10⁻² (standard stand-in for
vacuum). Honest persistence calibration up front: τ_resistive = μ₀σL², so
2 days at L = 1 cm needs σ ≈ 10¹⁵ S/m ≈ 2×10⁷ × copper — i.e. effectively
superconducting coherence; classical MHD can only tell us the *shape* of the
decay and its S-scaling, which we extrapolate.

### 7.4 The two experiment paths

**Path A — fluid re-run of v1 (fractal imposed, fluid physics).** Fractal
coil field from the v1 Biot–Savart generator, density blob following the
tube, v = 0 (in MHD the current is carried by ∇×B; no seeded flow needed).
Measures τ(k) again where the decay channels are reconnection + resistivity.
Sharp question: does the nested coil reconnect internally into a plain torus
(Taylor relaxation predicts helicity-preserving simplification), and does
τ(k) still fall with k?

**Path B — self-assembly from two tori (the genesis claim).** A
two-parameter family of double-ring initial conditions, each ring built from
a vector potential (∇·B = 0 by construction) with core radius a, major
radius R, separation d, and a twist parameter; plus optional hydrodynamic
vortex-ring velocity of either sense:

1. *Counter-helicity magnetic ring pair* (co-directed ring currents,
   opposite toroidal field): laboratory analogue merges into an FRC — also
   our code-validation case against known phenomenology.
2. *Opposed ring currents* pushed together by initial inflow: forced
   reconnection between "two toruses with opposing currents".
3. *Counter-rotating hydrodynamic vortex rings in head-on collision* with a
   weak frozen-in seed field: the Lim–Nickels configuration. The v2
   headline question: do the secondary rings inherit twisted flux — i.e.
   does a **ring of magnetized sub-rings (a 2-level FTM) self-assemble**?

### 7.5 Diagnostics and falsifiable outcomes

Magnetic, kinetic, and cross helicity budgets; energy partition; dipole
m(t) and anapole T(t) — self-assembly of an FTM would announce itself as
growth of |T| from an initial condition with |T| ≈ 0; azimuthal mode
spectrum of current/vorticity density on the mid-plane (counts the number of
self-assembled sub-rings, comparable to Lim–Nickels mode numbers); field-line
tracing for twist/writhe; lifetime vs Lundquist number τ(S) at fixed
geometry, extrapolated toward the 2-day requirement. Outcomes that would
*support* the claim: spontaneous |T| growth, sub-ring formation with
inherited twist, τ growing faster than linearly in S. Outcomes that would
*refute* it in this model family: azimuthal breakup without magnetic
nesting, monotone helicity-preserving simplification to a single torus,
τ ∝ S or slower.

### 7.6 First v2 results

Implementation notes that turned out to matter: plain first-order Rusanov
destroyed the few-cell ring cores within an Alfvén time, so the fluxes use
MUSCL/minmod reconstruction (second-order); and the toroidal (twist)
component of the flux rings must also be built from a vector potential or
its *discrete* divergence pollutes the run.

64³ scenario suite (S = 500, t = 15 Alfvén times, videos in `out/videos/`):

- **No spontaneous anapole moment in any scenario**: |T| stays at machine
  zero throughout. Cause identified: the two-ring initial conditions are
  axisymmetric, and a clean grid gives azimuthal instabilities nothing
  physical to grow from — the azimuthal mode spectrum reads back the
  Cartesian grid's own m = 4 anisotropy (m = 8 harmonic in the vortex-ring
  collision) at the few-percent level. Lesson: the Lim–Nickels breakup is
  noise-seeded in real fluids, so the simulation must seed symmetry-breaking
  noise explicitly (added: 2% random velocity noise, deterministic seed).
- *counterhel* reproduces merging phenomenology: the rings attract, drive a
  current layer at the mid-plane, and reconnect into a single object with
  surviving net dipole (m_z ≈ 0.26); magnetic energy drops ~100× by t = 15
  (reconnection + residual numerical diffusion at 3-cell cores).
- *opposed*: the anti-parallel ring currents annihilate; the net dipole
  cancels to grid zero.
- *limnickels*: kinetic energy decays ~60× with the seed field passively
  advected; no magnetized sub-rings without noise seeding.

96³ noise-seeded runs (2% velocity noise, t = 18, videos
`v2_*_N96*.mp4`): **reconnection does spontaneously generate anapole
current structure — the first self-assembly signal of the project — but at
the sub-percent level.**

- *limnickels*: |T| grows smoothly from machine zero (7×10⁻¹⁷) to
  3.6×10⁻⁴, saturating by t ≈ 10 and holding; the weak seed field is wound
  into a net dipole (m_z: 4×10⁻⁶ → 0.07); the collision annulus breaks up
  at azimuthal mode m = 4 with relative amplitudes up to ~0.5 — well above
  the unseeded grid-anisotropy floor, though m = 4 is also the grid's
  preferred symmetry, so the mode *number* needs a rotated-IC or
  higher-resolution cross-check before it is trusted.
- *counterhel*: a transient anapole peaks at ~1.2×10⁻⁴ around t ≈ 8
  (mid-merge) and relaxes to ~1×10⁻⁴; the merged object keeps a slowly
  decaying dipole m_z ≈ 0.28. First 3D volume-rendered video stream
  (`v2_counterhel_N96_3d.mp4`).
- Time-series plots (SVG, from `scripts/plot_anapole.jl`):
  ![anapole growth](out/plots/anapole_T.svg)
  ![dipole winding](out/plots/dipole_mz.svg)
- Scale honesty: the self-assembled anapole fraction is |T|/|m| ≈ 0.6%
  (limnickels) and ~0.04% (counterhel), versus ~50% for the hand-built v1
  fractal coil. Self-assembly of *weak* FTM character from anapole-free
  initial conditions: observed. Self-assembly of an actual fractal
  toroid: not at these parameters — the follow-up knobs are stronger seed
  field, higher resolution, and longer runs to test whether the saturated
  |T| is a plateau or a slow-growth phase.

**192³/256³ campaign (2026-07-19, GPU).** The higher-resolution, longer-time
follow-up called for just above.

- **Persistence — the plateau is real.** At 192³, limnickels |T| saturates
  at ~9.4×10⁻⁴ and holds flat from t ≈ 24 to t = 36 while kinetic energy
  falls 3.5× and magnetic energy 4.4×: a plateau, not a slow-growth phase, on
  this timescale. Since |T| = (1/10)∫[(r·J)r − 2r²J]dV is a functional of the
  current J = ∇×B, a constant |T| means a *current* persists — it is the
  fluid flow (E_kin) that dies, not the current. The field energy ∫B² also
  decays 4.4×, yet under a uniform decay B→fB one would have E_mag/|T|²
  constant; measured it drops ~5× (t = 20→36), so the decay is
  scale-selective. That is what resistive dissipation (∝ η k²) does — small
  scales first, leaving the largest-scale (lowest-k) toroidal current, which
  is exactly what the anapole measures. The anapole is the long-lived,
  large-scale survivor; it must eventually decay resistively (the longest
  timescale in the box), which t = 36 does not reach.
- **Mechanism — strength tracks the kinetic drive.** Peak |T| at 192³ orders
  with the vortex-ring circulation P0: limnickels (0.40) 9.4×10⁻⁴ > opposed
  (0.30) 5.9×10⁻⁴ (still rising at t = 18) > counterhel (0.10, magnetically
  dominated) 1.6×10⁻⁴. All three self-assemble a persisting anapole, so it is
  generic to the two-ring geometry — but its magnitude is set by the kinetic
  collision, not the field configuration (counterhel has the strongest fields
  and the weakest anapole).
- **Convergence — not there yet.** The anapole magnitude is not
  resolution-converged: 256³ |T| runs 2.4–5× below 192³ through the growth
  phase (the gap narrows with time). Bulk energetics agree to ~14%, so the
  disagreement is specific to the anapole. The rotated-IC / higher-resolution
  cross-check flagged for the 96³ mode number comes back negative: the
  azimuthal mode is grid-influenced — m = 4 at 192³, flickering m = 4/m = 8
  at 256³, m = 8/12 in opposed, all harmonics of the cubic grid's 4-fold
  symmetry — so the sub-ring *count* is not yet physical. A resolution ladder
  (96–256³ to saturation) and a 192³ noise-seed ensemble are running to test
  whether the *saturated* |T| converges and is seed-independent.

Bottom line: "a spontaneous anapole forms and outlives the driving flow" is
robust across resolution and initial condition; any specific value for its
strength or the sub-ring count is not yet converged.

![v2 scenarios at t = 12, 192³ (3D volume render: opacity = |B|, colour =
|ω|). counterhel keeps two distinct rings merging toward an FRC; limnickels
has collapsed to a compact core; opposed stays diffuse — the same two-ring
geometry gives visibly different outcomes.](out/figures/ic_comparison_N192.png)

### 7.7 Numerical plan and reuse

New module `src/mhd.jl` (state: ρ, ρv, B, ψ — 8 cell-centered fields;
Rusanov/HLL fluxes with fast-magnetosonic wave speed; GLM cleaning; explicit
resistive/viscous terms), driven by the existing SSP-RK3 harness. Reused
unchanged: `Box`, sponge layer, PNG writer, inline frame rendering, video
stitcher, moment/helicity diagnostics, fractal-coil generator (path A and
the ring construction). v1 code stays untouched; v2 scripts are
`scripts/v2_*.jl`. Validation tests before experiments: circularly polarized
Alfvén wave (exact solution, speed + resistive damping rate γ = ηk²),
magnetic-loop advection (∇·B control), hydrodynamic vortex-ring
self-propagation, and the counter-helicity merging → FRC phenomenology
check. Same GPU trajectory as v1 (the 8-field cell-centered layout is, if
anything, more CUDA-friendly than the staggered v1 grid).

## 8. v3 direction (planned): potentials-first, coherent matter

Decision (2026-07-18): v3 will be formulated **closer to the source
claims** — at minimum with the electromagnetic **potentials (φ, A) as the
dynamical variables**, per the analysis in neutrinos.md §3–§3.1. The
natural and honest way to do this is one construction, not two:

- **Matter**: a Gross–Pitaevskii condensate ψ (the continuum form of
  "coherent matter", neutrinos.md §4) with minimal coupling — the GPE
  couples to φ and **A** directly, `(−i∇ − qA)²ψ/2m + qφψ + g|ψ|²ψ`;
  fields never enter the matter equation. Potentials-first is not an
  add-on here; it is how a condensate couples to electromagnetism at all.
- **EM sector**: Lorenz-gauge wave equations □A = J, □φ = ρ (same
  hyperbolic machinery as our existing solvers), with E, B derived
  quantities used only for comparison diagnostics.
- **Native observables**: condensate phase (an Aharonov–Bohm-sensitive
  readout), quantized circulation, persistent toroidal currents —
  persistence is *generic* in this model family rather than the puzzle.
- **Experiments this unlocks**: (a) the anapole null, quantified — how
  completely a fractal toroid's external *fields* cancel while its
  external *potentials* survive; (b) the coherent-receiver test — an
  oscillating anapole next to a condensate probe: does the probe's phase
  register what no classical field carries? (the simulatable core of
  Greenyer's "potential waves" + "coherent matter" pairing); (c) whether
  a coil-of-coils of quantized vortex lines is stable, metastable, or
  reconnects away.

v1 (kinetic, field-based) and v2 (MHD, field-based) remain as the
classical control group against which the potentials/coherence
formulation is compared.
