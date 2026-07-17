# Relic neutrinos, coherent matter waves, and Matsumoto
### Research notes toward simulatable physics

Companion to README §2/§7 and `reconn.md`. Three threads that recur in Bob
Greenyer's writings, each examined for (a) what the real physics says and
(b) what part of it could become a continuum-mechanics simulation.

---

## 1. Relic neutrinos: what they actually are

The **cosmic neutrino background (CνB)** is a firm prediction of standard
cosmology: neutrinos decoupled from matter ~1 second after the Big Bang
and have streamed freely since. Today:

- Temperature ≈ **1.95 K** (cooler than the CMB by the (4/11)^{1/3}
  factor from electron–positron annihilation).
- Number density ≈ **56 /cm³ per species** → ~**336 /cm³** total for six
  species (ν and ν̄ of three flavors). They are everywhere, including in
  every experiment's apparatus.
- Typical momenta ~10⁻⁴ eV/c: at least two mass states are
  **non-relativistic today**, moving at ~10⁻³ c.
- **Never directly detected** — cross-sections at these energies are of
  order 10⁻⁵⁶ cm² or below; the proposed PTOLEMY-style capture
  experiments and [nuclear-spin approaches (arXiv 2508.20357)](https://arxiv.org/pdf/2508.20357)
  are still at the concept stage. Gravitational clustering enhances the
  local density only mildly
  ([Ringwald & Wong, hep-ph/0408241](https://arxiv.org/pdf/hep-ph/0408241)).

One genuinely striking property, and probably why the CνB attracts
speculation like Greenyer's: their **de Broglie wavelengths are
macroscopic** (mm-scale), and at least part of the background plausibly
exists as spatially extended coherent superpositions
([*Quantum Coherence of Relic Neutrinos*, arXiv 0811.4370](https://arxiv.org/pdf/0811.4370)).
So "a macroscopically coherent neutrino background permeating everything"
is, in a limited technical sense, real physics.

## 2. Greenyer's claim, and the energy-budget arithmetic

Greenyer proposes ([THOR post](https://remoteview.substack.com/p/thor-outside-of-the-inside))
that the EVO axis carries "a coherent gravitational wave beam of spinning
dark matter, likely relic neutrinos", through which structures link and —
elsewhere — that stable matter "breathes energy from the relic neutrino
background". Could that sustain a fractal toroid? The honest arithmetic:

- Kinetic energy density of the CνB: ~336 cm⁻³ × ~10⁻⁴ eV ≈ **10⁻¹⁴
  J/m³** — about 13 orders of magnitude below the energy density of
  Earth's magnetic field (~10⁻¹ J/m³ at the surface... itself modest).
- Even granting total absorption of the incident flux, the power through
  a 1 cm² object is ~10⁻¹² W; with realistic weak-interaction
  cross-sections it drops ~40 further orders of magnitude.
- The only route around this is coherent scattering off macroscopic
  domains (N² enhancement) — the physics behind Weber's 1980s claims of
  crystal-coherent neutrino detection, which did not survive scrutiny.

Conclusion: no known mechanism lets the CνB power anything macroscopic.
But the claim can still be made *quantitative and falsifiable within our
simulations* — see §5.

## 3. Coherent matter waves: what the term means

Real physics ladder:

1. **De Broglie matter waves** — any particle has wavelength λ = h/p;
   interference observed for electrons up to macromolecules.
2. **Coherent matter** — many particles occupying a single quantum state
   with one macroscopic wavefunction: Bose–Einstein condensates,
   superfluid He, superconducting Cooper pairs, atom lasers. (The
   matter-wave device patent that surfaces in searches,
   [US 6476383](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/6476383),
   is mainstream atom optics of this kind.)
3. **Preparata's "QED coherence in matter"**
   ([book](https://books.google.com/books/about/QED_Coherence_in_Matter.html?id=u-MvobTFGLEC)) —
   the contested claim that ordinary condensed matter contains
   *room-temperature* "coherence domains" oscillating in phase with a
   trapped EM field. This is almost certainly what Greenyer's "coherent
   matter waves at any temperature" refers to
   ([his interview on coherent matter](https://archive.org/details/untitled_20211019_2152));
   it is the standard theoretical scaffolding in LENR circles, and it is
   not accepted mainstream physics.

**The continuum-mechanics door.** A coherent matter wave has an exact
continuum description: the **Gross–Pitaevskii equation** (nonlinear
Schrödinger), which under the Madelung transform *is* a barotropic fluid
with an extra "quantum pressure" term — zero viscosity, and vorticity
locked into **quantized vortex filaments**. Two properties make it the
natural v3 model family for this project:

- **Persistent currents are generic**: a superfluid flow in a torus does
  not decay (no viscosity, no resistivity) — the *only* continuum setting
  where days-long persistence is the default rather than the puzzle.
- Quantized vortex rings/knots reconnect and cascade much like their
  classical cousins (cf. `reconn.md`; Kleckner–Irvine and the BEC
  vortex-knot literature), so "nested vortex structure in a coherent
  medium" is directly expressible.

## 4. Matsumoto: the observational corpus Greenyer leans on

Takaaki Matsumoto (Hokkaido University) ran electrolysis and discharge
experiments through the 1990s and recorded products on **nuclear
emulsions**, reporting ring-shaped, spiral, and **mesh-like traces**
([*Observation of meshlike traces on nuclear emulsions during cold
fusion*, OSTI](https://www.osti.gov/biblio/6364269);
[ordinary-water experiments, lenr-canr PDF](https://lenr-canr.org/acrobat/MatsumotoTcoldfusionb.pdf);
[Fusion Technology 24, 296 (1993)](https://www.tandfonline.com/doi/abs/10.13182/FST93-A30205);
[collected papers 1989–1999](https://www.amazon.co.jp/-/en/Dr-Takaaki-MATSUMOTO/dp/B0B6XS3M7D)).
He interpreted them via his **"Nattoh model"** — degenerate atom clusters
("itons"), "quad-neutrons", "gravity decay", micro black/white holes —
and described **"micro ball lightning"** emerging from his cells.
Greenyer regards him as a key precursor
([LENR-forum thread](https://www.lenr-forum.com/forum/thread/6483-the-outstanding-work-of-takaaki-matsumoto/)),
reading his ring traces as EVO/FTM imprints.

The theoretical superstructure (quad-neutrons, gravity decay) is not
physics we can instantiate. What *is* engageable is the **morphology**:
rings, rings-with-satellites, and meshes are exactly the shapes that
vortex-ring collisions and instabilities produce (see `reconn.md`), and
their multiplicities are countable in both his emulsions and our runs.

## 5. Points of connection → concrete simulation subjects

| Claim/observation | Simulatable counterpart | Status |
|---|---|---|
| CνB "feeds" the toroid | Add a phenomenological volumetric power/momentum source on the torus axis in the v2 MHD runs; measure the **drive power P(k, S) required for τ → ∞**, then compare with the CνB budget of §2 (≤10⁻¹² W/cm²). If required P exceeds the budget by tens of orders of magnitude at any plausible scale, the claim is quantitatively closed in this model family. | Not yet run; small change to `mhd.jl` (source term) |
| "Coherent matter waves at any temperature" | **Gross–Pitaevskii (v3)**: toroidal persistent current carrying nested/knotted quantized vortices on the existing grid infrastructure; measure whether a "coil of coils" of quantized vortices is stable, metastable, or decays by reconnection. GPE is a continuum PDE — the honest home of "indefinite persistence". | v3 candidate, natural next model family |
| Matsumoto's ring & mesh traces | Azimuthal mode multiplicities from our `azimuthal_spectrum` diagnostic (ring → necklace mode number) compared against multiplicities countable in his emulsion figures; mesh traces vs late-stage cascade patterns. | Diagnostic exists; comparison awaits a trace census |
| "Azimuthal beam" linking EVOs into filaments | Two toroids + axial momentum source in MHD; does a stable filament/bridge form? | Variant of path B, cheap |

**Reading order note.** §1 is solid physics; §3's rungs 1–2 are solid;
everything else above is either speculative (Greenyer, Preparata,
Matsumoto's models) or our own proposed tests. The table is the part
intended to survive contact with the simulations.
