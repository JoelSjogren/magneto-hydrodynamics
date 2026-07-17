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

- **Simulation of the Lim & Nickels configuration** (watch the necklace
  form at the end):

  <iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/Qxr7tsZUy1c" title="Head-on vortex ring collision" frameborder="0" allowfullscreen></iframe>

  Link: [Head-on vortex ring collision (YouTube)](https://www.youtube.com/watch?v=Qxr7tsZUy1c)
- Original experiment: T. T. Lim & T. B. Nickels, *Instability and
  reconnection in the head-on collision of two vortex rings*, Nature 357,
  225 (1992).
- Animated-GIF treatments at FYFD:
  [Vortex Ring Collisions (2018)](https://fyfluiddynamics.com/2018/07/one-of-the-most-enduringly-popular-submissions-i/) ·
  [Vortex Ring Collision (2012)](https://fyfluiddynamics.com/2012/01/two-vortex-rings-collide-head-on-in-this-video-if/)
- The famous SmarterEveryDay slow-motion water-tank realization, with
  build details: [When Vortex Rings Collide (Hackaday writeup + embedded video)](https://hackaday.com/2018/06/22/when-vortex-rings-collide/) ·
  [Boing Boing writeup](https://boingboing.net/2018/06/22/watch-this-absolutely-glorious.html)
- Numerical anatomy of the same event: [*A Cascade Leading to the
  Emergence of Small Structures in Vortex Ring Collisions* (arXiv
  1802.09973, figures viewable in-browser)](https://ar5iv.labs.arxiv.org/html/1802.09973)

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
| Necklace of secondary rings | azimuthal mode spectrum at the collision annulus | grid m=4 only before noise seeding; 96³ noise-seeded run in progress |
| Iterated cascade (rings³) | mode spectrum vs time + frame videos | — |
| Writhe→twist helicity flow | anapole moment |T|(t), field-line twist | |T| ≈ 0 so far |
| Counter-helicity merge → FRC | mid-plane current sheet, m_z survival | reproduced at 64³ |
