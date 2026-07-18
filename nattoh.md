# Matsumoto's Nattoh corpus: an image gallery

Companion to `neutrinos.md` §5 and `reconn.md`. Takaaki Matsumoto (Dept.
of Nuclear Engineering, Hokkaido University) spent the 1990s recording
what he interpreted as "micro ball lightning" and its imprints — on
nuclear emulsions, cathode surfaces, and acrylite sheets — and explaining
them with his **"Nattoh" model** (*Fusion Technology* 16, 532 (1989)):
hydrogen clusters condensing into degenerate "itons", with
"electro-nuclear collapse", "quad-neutrons", and "gravity decay" as the
downstream vocabulary. The theory is not something we can instantiate;
the **morphology is** — rings, rings-with-satellites, meshes, strings,
and trails are exactly the shapes vortex-ring dynamics produces (see
`reconn.md`), and their multiplicities are countable in both his scans
and our simulations.

Images below were extracted from scanned papers hosted by
[lenr-canr.org](https://lenr-canr.org) and live (deliberately
uncommitted) under `out/nattoh/`; the PDFs themselves are cached in
`out/nattoh/pdf/`. View this page through the local preview so the
relative image paths resolve.

---

## 1. Experimental scans

### Mesh-like trace on nuclear emulsion

![mesh trace montage](out/nattoh/01_mesh_trace_montage.png)

A ~300 μm mesh-like etch region (50 μm scale bar), stitched by Matsumoto
from micrographs. He classified these meshes into ring-associated and
ring-free types and attributed them to "gravity decay of quad-neutrons"
(*[Cold fusion experiments with ordinary water and thin nickel foil,
technical note](https://lenr-canr.org/acrobat/MatsumotoTcoldfusionb.pdf)*;
cf. [*Observation of meshlike traces on nuclear emulsions*, OSTI
6364269](https://www.osti.gov/biblio/6364269)). To a fluid-dynamical eye
the mesh reads as a late-stage tangle — compare the post-breakup frames
of our path-B runs.

### Ring trace with film fragments

![ring trace and films](out/nattoh/02_ring_trace_and_films.png)

The canonical **ring trace** (dark circle, right of center; ~30 μm) amid
polygonal film-like fragments. Rings are the signature Matsumoto (and
after him, Greenyer) read as the footprint of a toroidal object — the
same annulus-with-substructure morphology our azimuthal-mode diagnostic
counts. Same source paper.

### Branched "string" trace

![string tree trace](out/nattoh/03_string_tree_trace.png)

A branched, tree-like filament with an adjacent closed loop — captioned
in the source among "gravity decay products with strings". In vortex
vocabulary: filaments, linking, and a detached ring; in Greenyer's:
EVOs linking into filaments (`neutrinos.md` §5's "azimuthal beam" row).

### Ring-clusters on cathode surfaces (SEM)

![ring clusters SEM page](out/nattoh/04_ring_clusters_sem_page.png)

A full page from his *Journal of New Energy* work: deposited material on
a Ni cathode with a "tiny hole" (Fig. 7), **ring-clusters on an Fe
cathode** (Fig. 8) — the structures he proposed as the *generators* of
micro ball lightning — and hexagonal "decay products" (Fig. 9)
(*[experiments PDF](https://lenr-canr.org/acrobat/MatsumotoTexperimenta.pdf)*).

### Tiny sparks in flight

![tiny sparks page](out/nattoh/05_tiny_sparks_page.png)

Voltage-dependent "tiny sparks" on a Pd cathode and a 33 ms/frame
sequence of sparks *separated from the electrode* — his claimed direct
observation of free-flying micro ball lightning, the laboratory-scale
version of the detachment-and-persistence claim this project simulates.

### Trail marks

![trail marks](out/nattoh/06_trail_marks.png)

Trail marks (100 μm scale bar) of the kind described as ball lightning
that "slid and hopped" across the recording medium — reproduced in E. H.
Lewis's survey
(*[Microscopic Ball Lightning](https://www.lenr-canr.org/acrobat/LewisEmicroscopi.pdf)*),
which collects and compares Matsumoto's markings with similar reports.

## 2. Theoretical drawings — a gap

The Nattoh model's schematic drawings (iton structure, the hydrogen
"natto-like" cluster condensation, quad-neutron decay chains) are in
papers not freely mirrored: *"Nattoh" Model for Cold Fusion* (Fusion
Technology 16, 532 (1989)), *Observation of Quad-Neutrons and Gravity
Decay* (Fusion Technology 19, 2125 (1991)), and the collected volume
*[Steps to the Discovery of Electro-Nuclear
Collapse](https://www.amazon.co.jp/-/en/Dr-Takaaki-MATSUMOTO/dp/B0B6XS3M7D)*
(1989–1999). TODO: source these scans; until then this gallery is
experimental-only.

## 3. How the MFMP frames Matsumoto

- Greenyer curates his memory directly: a Remote View post reading the
  preface and biography of the collected papers
  ([Matsumoto — Steps to the Discovery of Electro-Nuclear
  Collapse](https://remoteview.substack.com/p/matsumoto-steps-to-the-discovery),
  audio; includes a photograph of Matsumoto in Sapporo, 1999), framing
  him as the great under-recognized experimentalist of the field.
- The LENR-forum thread [The outstanding work of Takaaki
  Matsumoto](https://www.lenr-forum.com/forum/thread/6483-the-outstanding-work-of-takaaki-matsumoto/)
  (Greenyer-adjacent) collects his emulsion imagery and connects it to
  EVOs/charge clusters.
- In Greenyer's synthesis (README §2), Matsumoto's rings/meshes are FTM
  imprints, his micro ball lightning is the EVO, and his ring-clusters
  are self-assembly sites — i.e. the MFMP reads this corpus as *the
  observational archive* for exactly the structures our v2/v3 runs try
  to grow.

## 4. What our simulations can say about these images

| Image | Simulatable counterpart | Where |
|---|---|---|
| Ring traces | footprint of a toroidal object; annulus diameter distribution | path-B objects, any run |
| Rings with satellites / mesh | azimuthal breakup mode number; late-stage tangle | `azimuthal_spectrum`, path-B videos |
| String/filament traces | filament linking between toroids | "azimuthal beam" variant, planned |
| Trail marks | propagating ring in a dissipative medium | vortex-ring propagation (already validated) |
| Sparks separated from electrode | detachment + persistence τ | the project's central measurement |

The honest caveat, always: etch marks on emulsions have mundane
candidate explanations (chemical damage, mechanical scratching,
radiochemistry) that Matsumoto's papers dismiss quickly. The gallery
documents what the claims *look like*, so that morphological comparison
with simulations is possible at all — resemblance alone will never carry
the argument.
