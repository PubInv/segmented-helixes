# Joint-face rotation addendum

This is a standalone working note for colleagues, dated October 1, 2026.
It supplements the Joint Face Normal Method in Robert L. Read's 2022
segmented-helix paper. It is suitable for placing beside the original source
in `MathematicJournal/`, while keeping the published manuscript intact.

## Contents

- `JointFaceRotationAddendum.tex`: editable, single-file LaTeX source.
- `JointFaceRotationAddendum.pdf`: compiled companion note.
- `verify_joint_face_rotation.py`: reproducible numerical checks, using only
  the Python 3 standard library.
- `verification-results.txt`: the output from those checks.

The note includes the theorem and proof, the precise twist reference,
face-normal classes, exact values for all 23 nonzero-radius catalog entries,
and the five straight-joint degeneracies. It records the correction that
zero rotation in a straight-joint catalog row does not determine the rotation
of a labeled subunit. It does not change the original paper or software.

## Build the PDF

In the directory containing the source, run:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error JointFaceRotationAddendum.tex
```

Alternatively, run this command twice:

```sh
pdflatex -interaction=nonstopmode -halt-on-error JointFaceRotationAddendum.tex
```

The source uses ordinary LaTeX packages, including AMS mathematics, TikZ,
booktabs, and hyperref. There are no external images, bibliography files,
custom fonts, or MDPI-class dependencies. It can also be uploaded as one
`.tex` file to Overleaf and compiled with pdfLaTeX.

## Run the checks

```sh
python3 verify_joint_face_rotation.py
```

The checks compare radical expressions with the original rounded catalog
values, enumerate normal-dot-product classes, test the shortest-alignment
baseline for centered regular solids, test matrix composition and reference
phase, and verify that all five straight-joint cases retain the stated full
rotation while their joint points remain collinear. No extra Python packages
are required. The proof in the note establishes the identity; the checks
provide independently reproducible numerical validation.

## Source provenance

Repository: https://github.com/PubInv/segmented-helixes

Source revision: `8be3839de2204b966c8aedb7c28f6a6f5c18076e`

Reviewed files:

- `MathematicJournal/StackingHelix.tex`
- `js/segment_helix_math.js`, especially `adjoinPrism`, `PointAxis`, and
  `computeThetaU`
- `js/main.js`, especially the regular-solid face-center construction

Original paper: https://doi.org/10.3390/math10142533

## Possible repository placement

Copy these files into `MathematicJournal/` or a new companion-note directory.
Add a link to the PDF from the repository's existing README. Build outputs
such as `.aux`, `.log`, `.out`, `.fls`, and `.fdb_latexmk` need not be committed.
If the note is later incorporated into a revised manuscript, its mathematics
can be moved into the original body after the normal-method construction;
the standalone preamble and bibliography should remain separate.

The draft credits assistance from ChatGPT in a preparation note. Authorship,
wording, and any release decision remain available for Robert to review.
