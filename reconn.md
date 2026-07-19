# Vortex reconnection → nested toroidal structures: a visual catalogue

Curated graphical/animated resources on vortex reconnection as it pertains
to the **self-assembly of nested ("ring of rings") toroidal structures** —
the process the v2 path-B experiments (README §7) try to reproduce in MHD.
Inline video embeds play in the local preview; each entry also carries
plain links.

---

## 1. Head-on ring collisions → a ring of secondary rings

The canonical demonstration that "ring → ring of rings" self-assembly is
real fluid dynamics: two counter-rotating vortex rings collide head-on,
stretch into a thin annulus, go azimuthally unstable, and **reconnect into
a necklace of small secondary rings** that fly off radially.

- **Simulation of the Lim & Nickels configuration** — watch the necklace
  form at the end ([direct link](https://www.youtube.com/watch?v=Qxr7tsZUy1c)):

  <iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/Qxr7tsZUy1c" title="Head-on vortex ring collision" frameborder="0" allowfullscreen></iframe>

- Original experiment: T. T. Lim & T. B. Nickels, *Instability and
  reconnection in the head-on collision of two vortex rings*, Nature 357,
  225 (1992) — animated-GIF summary at
  [FYFD](https://fyfluiddynamics.com/2018/07/one-of-the-most-enduringly-popular-submissions-i/).
- Slow-motion water-tank realization:
  [SmarterEveryDay, via Hackaday (build details + video)](https://hackaday.com/2018/06/22/when-vortex-rings-collide/).
- Numerical anatomy: [*A Cascade Leading to the Emergence of Small
  Structures in Vortex Ring Collisions* (arXiv 1802.09973)](https://ar5iv.labs.arxiv.org/html/1802.09973)

**Relevance to v2:** this is scenario `limnickels`. The azimuthal mode
counter in our diagnostics (`azimuthal_spectrum`) is meant to detect
exactly this necklace; the open question is whether the secondary rings
inherit *twisted magnetic flux* — which would make the necklace a 2-level
fractal toroidal structure.

## 2. The iterative cascade: rings of rings of rings

McKeown et al. showed the breakup is not one-shot: each generation of
counter-rotating filaments can undergo the **elliptical instability
again**, producing a *hierarchy* of ever-smaller ring structures — the
closest mainstream physics gets to a self-assembled fractal torus cascade.

- **McKeown, Ostilla-Mónico, Pumir, Brenner & Rubinstein**, *Turbulence
  generation through an iterative cascade of the elliptical instability*,
  Science Advances 6, eaaz2717 (2020) — **open access with supplementary
  movies** of the cascade generations:
  [paper + movies](https://www.science.org/doi/10.1126/sciadv.aaz2717) ·
  [arXiv PDF](https://arxiv.org/pdf/1908.01804)
- Follow-up with sustained collisions: [*Turbulence through sustained
  vortex ring collisions*, Phys. Rev. Fluids 8, 110507 (2023)](https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.8.110507)
- Breakup mechanics across regimes: [*Instability and disintegration of
  vortex rings during head-on collisions and wall interactions* (arXiv
  2107.12324)](https://arxiv.org/pdf/2107.12324)

**Relevance to v2:** if Greenyer's "rings cluster into rings of rings" has
a classical mechanism, this iterated elliptical instability is the
candidate. In MHD the question becomes whether magnetic tension arrests
the cascade or organizes it.

## 3. Knotted and linked vortices: reconnection and helicity across scales

Reconnection is what lets topology change — links and knots untie into
separate rings while (approximately) conserving helicity by converting
linking into coiling and twisting: **structure cascading from large
writhe to small-scale twist**, i.e. helical fine structure appearing on a
ring spontaneously.

- **Kleckner & Irvine (U. Chicago)**, *Creation and dynamics of knotted
  vortices*, Nature Physics 9, 253 (2013) — supplementary videos:

  <iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/YCA0VIExVhg" title="Knotted vortices in real fluids" frameborder="0" allowfullscreen></iframe>

  Links: [Knotted vortices in real fluids (YouTube)](https://www.youtube.com/watch?v=YCA0VIExVhg) ·
  [3D version (YouTube)](https://www.youtube.com/watch?v=rcnw8NeJqjU) ·
  [paper PDF](https://faculty.ucmerced.edu/dkleckner/papers/KnottedVortices.pdf) ·
  [UChicago story](https://news.uchicago.edu/story/vortex-loops-could-untie-knotty-physics-problems) ·
  [phys.org with video](https://phys.org/news/2013-03-physics-duo-vortex-fluid-video.html)
- Helicity transfer during reconnection (linking → coiling → twisting):
  [*Helicity conservation by flow across scales in reconnecting vortex
  links and knots* (arXiv 1404.6513)](https://arxiv.org/pdf/1404.6513) ·
  [*The Life of a Vortex Knot* (arXiv 1310.3321)](https://arxiv.org/pdf/1310.3321)
- Detailed reconnection scaling: [*Scaling of Navier–Stokes trefoil
  reconnection* (arXiv 1610.00398)](https://arxiv.org/pdf/1610.00398)

**Relevance to v2:** helicity flowing from writhe (large-scale linking)
into twist (small-scale coiling) is precisely the "coil acquiring
sub-coils" direction of the FTM story — and it is measurable in our runs
via the twist/writhe diagnostics.

## 4. The magnetic counterpart: merging plasmoids

The MHD version of ring reconnection, realized in laboratories: two
spheromaks (magnetized tori) merge through reconnection, and the outcome
is set by relative helicity — counter-helicity merging annihilates the
toroidal field and relaxes into a **field-reversed configuration**.

- Swarthmore SSX, *Three-dimensional MHD simulations of counter-helicity
  spheromak merging* — figures of the full merge sequence:
  [PDF](https://plasma.physics.swarthmore.edu/brownpapers/DoubletCT.pdf)
- Princeton MRX (Magnetic Reconnection Experiment) — experiment site with
  image/movie galleries: <https://mrx.pppl.gov>
- FRC formation by spheromak merging (overview):
  [OSTI record](https://www.osti.gov/pages/biblio/1309540)

**Relevance to v2:** scenario `counterhel` is exactly this configuration;
our 64³ run already reproduces the mid-plane current layer and merge (see
`out/videos/v2_counterhel_N64.mp4`).

## 5. What to look for in our own runs

| Catalogue phenomenon | v2 diagnostic | Seen so far |
|---|---|---|
| Necklace of secondary rings | azimuthal mode spectrum at the collision annulus | noise-seeded collision does break up azimuthally, but the dominant mode is grid-dependent — m=4 (192³), m=4/8 flicker (256³), m=8/12 (opposed), all harmonics of the grid's 4-fold symmetry; not yet a physical count (§6) |
| Iterated cascade (rings³) | mode spectrum vs time + frame videos | — |
| Writhe→twist helicity flow | anapole (toroidal) moment |T|(t) — zero for a plain ring, nonzero for coil-of-coils winding | machine zero pre-noise; **noise-seeded collision grows a persisting anapole — 192³ saturates at ~9.4×10⁻⁴ and outlives the flow (E_kin/E_mag decay 3.5–4.4×); strength ∝ vortex drive; magnitude not yet resolution-converged (§6)** |
| Counter-helicity merge → FRC | mid-plane current sheet, m_z survival | reproduced at 64³ |

## 6. Results — spontaneous-anapole campaign (2026-07-19)

GPU runs at 96³–256³ on the `limnickels`, `counterhel`, and `opposed`
scenarios (full numbers and discussion in README §7.6). Headline:

- **The anapole persists and saturates.** `limnickels` |T| holds flat at
  ~9.4×10⁻⁴ (192³) from t≈24 to t=36 while kinetic and magnetic energy decay
  3.5–4.4×. |T| is a functional of the current J=∇×B, so this is a persisting
  *current* — the largest-scale toroidal one, the slow survivor of
  scale-selective resistive decay (∝ηk²). It is the fluid *flow* that dies,
  not the current.
- **Strength ∝ kinetic drive.** Peak |T| orders with vortex-ring circulation:
  limnickels (P0=0.40) 9.4×10⁻⁴ > opposed (0.30) 5.9×10⁻⁴ > counterhel (0.10)
  1.6×10⁻⁴. Generic to the two-ring geometry, but a vortex-collision effect,
  not a magnetic one (counterhel has the strongest fields, the weakest
  anapole).
- **Not yet converged.** 256³ |T| runs 2.4–5× below 192³ (gap narrowing with
  time); the breakup mode is grid-harmonic (m=4/8/12), not physical. A
  resolution ladder (96–256³ to saturation) and a 192³ noise-seed ensemble
  are running to test whether the *saturated* |T| converges and is
  seed-independent.

![v2 scenarios, 192³ — initial conditions (top, t = 0) vs evolved (bottom,
t = 12); 3D volume render (opacity = |B|, colour = |ω|). Similar two-ring
starts diverge: counterhel keeps two distinct rings (FRC-forming merge),
limnickels collapses to a compact core, opposed stays diffuse.](out/figures/ic_comparison_N192.png)
