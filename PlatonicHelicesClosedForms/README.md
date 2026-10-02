# Closed Forms for Platonic Helices

Author-review draft for Robert L. Read and Public Invention, October 2, 2026.

This is a standalone companion to Read's 2022 paper. It proves a common law
for the principal subunit rotation `theta`, segment–axis angle `phi`, signed
axial advance `d`, and radius `r` of the helix through the joining-face
centroids. It gives unit-edge exact expressions for all 28 original catalog
inputs and three corrected icosahedral matching configurations.

For the repository, a separate directory such as
`MathematicJournal/PlatonicHelicesClosedForms/` would keep this companion and
its verification files together. The PDF and its single LaTeX source are
sufficient for reading and compilation. Nothing in this package requires an
external image or bibliography file.

## Compile the paper

Use a TeX installation with the standard packages listed in the source:

```sh
latexmk -pdf PlatonicHelicesClosedForms.tex
```

Alternatively, run `pdflatex PlatonicHelicesClosedForms.tex` twice. The exact
tables are embedded in this source. Generated table fragments are supplied
for checking and later editing; the manuscript does not load them.

## Verify the mathematics

The independent geometric verifier needs only Python 3:

```sh
python3 verify_closed_forms.py
```

It checks the unit-edge coordinate solids, every printed parameter in the
28-row original catalog after normalization, all 162 corrected distinct-face
matching constructions and their representative multiplicities, 3,948
continuous-twist cases, and the three corrected icosahedral rows. It constructs
rotations with matrices and measures radii with projected circumcircles.

The optional symbolic generator uses SymPy:

```sh
python3 -m pip install -r requirements.txt
python3 build_exact_catalog.py
```

This checks the displayed radicals against the theorem and checks the exact
vertex calculation of the corrected icosahedral base angle. It rewrites
`ExactCatalog.json` and the five `tables_*.tex` fragments. It does not rewrite
the manuscript body or its embedded tables. Both accompanying run logs are
included.

## Review priorities

1. The formulas for `r`, `d`, and `phi` assume two joint centers of the form
   `B = rho*n_B`, `C = rho*n_C`, with a common positive `rho`. The normal-only
   rotation law has broader scope. `L` is the joint-centroid distance within
   one solid; `d` is the advance along the screw axis.
2. The twist is measured from the shortest rotation sending the incoming
   outward normal to the negative receiving outward normal. The paper proves
   this agrees with the repository baseline for the Platonic cases. Signed
   advances follow the chosen directed axis; straight rows use its forward
   direction.
3. The original constructor lengths have circumradius one for four solids
   and edge length one for the cube. The paper gives exact conversion factors
   and reproduces all four original printed columns to their three-decimal
   precision. It does not attribute an explicit unit-edge declaration to the
   original table.
4. Straight joint orbits hide rotations of the subunit. The new `theta` column
   records the actual principal rigid-motion rotation; the original zero
   markers and their interpretation are stated separately.
5. For icosahedral face 8, the exact matching base is
   `beta = acos((3*sqrt(5)-1)/8) = 44.4775121859299...` degrees. The corrected
   twists are `beta-120`, `beta`, `beta+120`. The six faces in this class have
   three positive and three negative bases. Original `-75`, `45`, `165` degree
   inputs remain explicitly tabulated as literal configurations, alongside
   the corrected matching configurations.

The 28 names and their multiplicities are inherited catalog labels. Different
twists can produce the same straight joint orbit. Equal radii and advances
across different solids do not alone imply the same ordered helix. The paper
does not establish collision-free infinite solid chains or priority over
existing screw-motion and tetrahelix results.

Diagrams are left for the next author-guided revision. The face-normal
geometry, projected chord calculation, and corrected icosahedral vertex match
are natural subjects for those illustrations.

## Files and sources

- `PlatonicHelicesClosedForms.pdf` and `.tex`: paper and self-contained source.
- `catalog_data.py`: source catalog transcription and numbered unit-edge faces.
- `ExactCatalog.json`: exact expressions, fixed constants, and numerical values.
- `verify_closed_forms.py`: independent standard-library geometric verifier.
- `build_exact_catalog.py`: optional SymPy generator and symbolic checks.
- `tables_Tet.tex`, `tables_Cub.tex`, `tables_Oct.tex`, `tables_Dod.tex`,
  `tables_Ico.tex`: generated exact tables.
- `verification-results.txt`, `symbolic-results.txt`: accompanying check logs.
- `requirements.txt`: optional symbolic dependency.

Primary catalog: R. L. Read, *Mathematics* **10** (2022), 2533,
[doi:10.3390/math10142533](https://doi.org/10.3390/math10142533).
Implementation audited:
[PubInv/segmented-helixes](https://github.com/PubInv/segmented-helixes),
revision `8be3839de2204b966c8aedb7c28f6a6f5c18076e`.

Prepared with ChatGPT assistance in derivation, calculation, verification
code, and writing. This is a working draft for detailed author review.
