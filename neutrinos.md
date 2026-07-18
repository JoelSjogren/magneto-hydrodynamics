# Relic neutrinos, coherent matter waves, and Matsumoto
### Research notes toward simulatable physics

Companion to README §2/§7 and `reconn.md`. Three threads that recur in Bob
Greenyer's writings, each examined for (a) what the real physics says and
(b) what part of it could become a continuum-mechanics simulation.

---

## 1. Relic neutrinos: what they actually are

The **cosmic neutrino background (CνB)** is a prediction of standard
cosmology: neutrinos decoupled from matter ~1 second after the Big Bang
and have streamed freely since. Epistemic status, consistent with §3's
weighting: the CνB has **never been directly observed** — its existence,
temperature, and density are inferred through the same long cosmological
model chain this project treats with caution. That cuts against
over-confident use of the numbers below, but it cuts *harder* against any
proposal that leans on the CνB as a power source: the entity being
invoked is itself known only through the model stack. Standard values:

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
  J/m³** — about 11 orders of magnitude below the energy conventionally
  ascribed to Earth's magnetic field at the surface (B²/2μ₀ ≈ 10⁻³ J/m³
  for ~50 μT; on the charge-based accounting of §3 the same total sits in
  the core currents — either way the comparison stands).
- Even granting total absorption of the incident flux, the power through
  a 1 cm² object is ~10⁻¹² W; with realistic weak-interaction
  cross-sections it drops ~40 further orders of magnitude.
- The only route around this is coherent scattering off macroscopic
  domains (N² enhancement) — the physics behind Weber's 1980s claims of
  crystal-coherent neutrino detection, which did not survive scrutiny.

Note the logical form of this argument: it is an **internal-consistency
check**, not an appeal to consensus. *Granting* the claimant's own
presupposition (the model-derived CνB, with its model-derived density and
energy), the budget still fails by ~40 orders of magnitude. One does not
need to trust cosmology to run this argument — only to observe that the
proposal borrows its central entity from cosmology and then violates that
same framework's arithmetic.

Conclusion: no known mechanism lets the CνB power anything macroscopic.
But the claim can still be made *quantitative and falsifiable within our
simulations* — see §6.

## 3. Interlude: where does electromagnetic energy live? (the Mead objection)

A foundational objection to phrases like "the field's energy density",
raised in this project's context via Carver Mead's appearance on the same
podcast circuit as Greenyer
([DemystifySci episode](https://goodpods.com/podcasts/the-demystifysci-podcast-209266/gravitational-variation-in-the-speed-of-light-dr-carver-mead-caltech-30826283)).
Mead — whose *Collective Electrodynamics* rebuilds EM from superconducting
loops, flux quantization, and the potentials — formulates energy **not**
as ∫(E² + B²)/2 dV but as the coupling of charge-current to the
potentials: his
[*Collective Electrodynamics I* (PNAS 1997, open access)](https://pmc.ncbi.nlm.nih.gov/articles/PMC20992/)
writes W = ∫ **J·A** dV and ascribes the energy density to where the
*current* is; "it is the vector potential A, rather than the magnetic
field B, that has a natural connection with the quantum nature of matter"
(Aharonov–Bohm, flux quantization). On this accounting, potentials
multiplied by zero charge give zero energy: **no charge, no
electromagnetic energy** — the position recalled from the podcast.

How seriously to take it, in three layers:

1. **Wherever sources exist, the dispute is about *where*, not *how
   much*.** In magnetostatics the two accountings are mathematically
   identical in total: ∫ B²/2μ₀ dV = ½∫ J·A dV (integration by parts,
   fields decaying at infinity), and likewise ∫ ε₀E²/2 = ½∫ ρφ. An
   inductor's ½LI² is the same number computed either way. This
   ambiguity of *location* is acknowledged mainstream physics — Feynman
   Lectures II §27 concedes classical EM cannot say where field energy
   resides, only what the totals and flows are. So Mead's bookkeeping is
   a legitimate reformulation, not a crank position (though see
   [Parrott's critical review](https://www.math.umb.edu/~sp/mead.pdf)
   for a careful mainstream dissent on how much of EM it rederives).
2. **The accountings genuinely diverge for source-free radiation — and
   the divergence may be empirically undecidable.** Light in transit has
   energy located in vacuum (Poynting view) or is a direct
   emitter-to-absorber transaction with no energy "in flight"
   (Wheeler–Feynman direct-action lineage, where Mead's sympathies lie).
   The decisive observation is that **every detector is made of charge**:
   there is no charge-free way to sample "energy in the vacuum", so all
   local experiments — radiation pressure, calorimetry, interferometry —
   measure charge-side bookkeeping that both accountings reproduce
   identically. Precedent that "vacuum energy" claims evaporate under
   direct accounting: the Casimir force, long advertised as zero-point
   vacuum energy made manifest, is fully computable as interactions
   among the charges in the plates
   ([Jaffe 2005](https://arxiv.org/abs/hep-th/0503158)). The only
   argument that even *attempts* to locate energy in source-free space
   is gravitational/cosmological (radiation-era expansion driven by
   photon energy density), and this project deliberately assigns it low
   weight: it tests a long conjunction of model assumptions (GR + FRW
   homogeneity + the ΛCDM component inventory + Gyr-scale
   extrapolation), with residuals historically absorbed by postulating
   new components — an inference chain too indirect to settle a
   foundational bookkeeping question. Working position: treat the two
   accountings as empirically equivalent for anything we can actually
   measure or simulate, and prefer the one that stays closest to what
   is observed — charges moving.
3. **For this project the dispute is almost invisible — and where it
   isn't, our results side with Mead's intuition.** Every v2 MHD object
   carries its currents with it (J = ∇×B in the matter), so all energy
   budgets in this repo are numerically identical under either
   accounting. And the one place the views differ — can a *pure field*
   structure with no charge persist? — is precisely what v1 phase 2
   tested dynamically: prescribed-current vacuum fields disperse within a
   light-crossing, and only charge/matter-bearing configurations persist.
   "No charge, no enduring electromagnetic object" is both Mead's
   ontology and our simulation result. It also lands on our modeling
   choice from day one (README §1: the current must be carried by massive
   charge) and strengthens the superconductor-inspired GPE direction of
   §4, since Mead's entire construction starts from persistent currents
   in superconductors.

A cheap diagnostic can make the bookkeeping dispute *visible* in our
runs: track ∫B²/2 dV and ½∫J·A dV side by side — equal during
quasi-static evolution, diverging by exactly the boundary Poynting flux
during radiative episodes (added to the table in §6).

### 3.1 Where Greenyer falls on this divide

Emphatically on the potentials side — but on its *other wing*. Via
Shoulders he describes the EVO as an "ideal monopole oscillator"
generating **"vector and scalar potential waves"** with the fields
themselves absent, and FTMs that "self-cancel to create
non-electromagnetic boundaries similar to anapoles"
([Decoding EVOs](https://www.altpropulsion.com/decoding-evos-a-deep-dive-into-exotic-vacuum-objects/);
[THOR](https://remoteview.substack.com/p/thor-outside-of-the-inside)).
His mainstream anchors (Zel'dovich, Dubovik, non-radiating sources) are
precisely the corner of orthodox EM where external fields vanish but
external potentials do not.

The divide-within-the-divide: **Mead is a conservative
potentials-fundamentalist** — potentials are the real quantities but only
matter multiplied by charge (W = ∫J·A), so nothing is delivered where
nothing charged or coherent sits. **Greenyer treats potential waves as an
autonomous channel** (energy/structure propagating where E = B = 0) — the
Tesla–Bearden "scalar wave" lineage. On Mead's own accounting, a
potential wave arriving at an uncharged, incoherent receiver delivers
exactly nothing: the disciplined wing of the potentials camp is *hostile*
to the promiscuous wing.

The legitimate kernel connecting them is **quantum phase**. In field-free
regions potentials demonstrably do one thing: shift the phase of coherent
matter waves — the Aharonov–Bohm effect, measured in both vector and
scalar forms
([neutron interferometry](https://inis.iaea.org/records/p1g3f-dbz68)).
Mainstream work has even proposed exactly Greenyer's source as the
testbed: a **time-dependent anapole** is a non-radiating source whose
external potentials could imprint a time-dependent AB phase
([Nemkov, Basharin & Fedotov, arXiv 1605.09033](https://arxiv.org/pdf/1605.09033)) —
while extended "AB electrodynamics" analyses still find a scalar wave
[cannot deposit energy in a medium](https://arxiv.org/pdf/2302.10224).
The steelman of Greenyer's position is therefore internally elegant: his
two recurring themes — anapole sources and "coherent matter" — are
exactly the pair AB physics requires. A potential-wave story is coherent
*only if* both emitter and receiver are quantum-coherent; classically,
both ends vanish. This is the direct motivation for the v3 formulation
(README §8): potentials as the formulation variables, a condensate as the
matter — the one continuum setting in which Greenyer-style claims are
even expressible.

## 4. Coherent matter waves: what the term means

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

## 5. Matsumoto: the observational corpus Greenyer leans on

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

## 6. Points of connection → concrete simulation subjects

| Claim/observation | Simulatable counterpart | Status |
|---|---|---|
| CνB "feeds" the toroid | Add a phenomenological volumetric power/momentum source on the torus axis in the v2 MHD runs; measure the **drive power P(k, S) required for τ → ∞**, then compare with the CνB budget of §2 (≤10⁻¹² W/cm²). If required P exceeds the budget by tens of orders of magnitude at any plausible scale, the claim is quantitatively closed in this model family. | Not yet run; small change to `mhd.jl` (source term) |
| "Coherent matter waves at any temperature" | **Gross–Pitaevskii (v3)**: toroidal persistent current carrying nested/knotted quantized vortices on the existing grid infrastructure; measure whether a "coil of coils" of quantized vortices is stable, metastable, or decays by reconnection. GPE is a continuum PDE — the honest home of "indefinite persistence". | v3 candidate, natural next model family |
| Matsumoto's ring & mesh traces | Azimuthal mode multiplicities from our `azimuthal_spectrum` diagnostic (ring → necklace mode number) compared against multiplicities countable in his emulsion figures; mesh traces vs late-stage cascade patterns. | Diagnostic exists; comparison awaits a trace census |
| "Azimuthal beam" linking EVOs into filaments | Two toroids + axial momentum source in MHD; does a stable filament/bridge form? | Variant of path B, cheap |
| Mead: energy at currents vs. in fields (§3) | Track ∫B²/2 dV and ½∫J·A dV side by side each frame: equal in quasi-static phases, diverging by exactly the boundary Poynting flux during radiative episodes — the bookkeeping dispute made visible in one plot. | Diagnostic not yet implemented; A available from a Poisson solve of ∇²A = −J |

**Reading order note.** §1 is solid physics; §3 is a live foundational
debate — this project's working position is charge-side accounting, with
the cosmological counterargument noted and deliberately down-weighted;
§4's rungs 1–2 are solid; everything else above is either
speculative (Greenyer, Preparata, Matsumoto's models) or our own proposed
tests. The table is the part intended to survive contact with the
simulations.
