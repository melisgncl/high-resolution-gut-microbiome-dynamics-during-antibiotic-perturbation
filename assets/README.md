# assets/

The figures **as submitted**, one image per figure:
`figure_{1,2,3,4,5}_as_submitted.png` and `figure_{S1..S5}_as_submitted.png`.

They are here for one purpose: so a reader can put a panel rebuilt by
`make figures` next to the published version and compare them directly.

## These are whole figures, not panels

Three published panels have no code — **1A** (experimental schematic), **5B**
(the *rpsE* deletion mapped onto the 70S ribosome) and **5C** (colony
morphology on spectinomycin). They were made outside the analysis: a drawing, a
structural render, and a photograph.

They are **not** shipped as separate files. They exist here only inside
`figure_1_as_submitted.png` and `figure_5_as_submitted.png`. So the automated
build reproduces every panel that is derived from data, and those three are
reproduced only in the sense that the submitted composite is included for
reference.

## Where the rebuilt panels differ from the submitted ones

Two deliberate differences, both documented in the root README:

- **Figure 3A** — the submitted panel labels x as Im and y as Re; the code
  computes the opposite. Rebuilt as computed.
- **Figure 4B/4C** — the submitted panels read control indices with the
  colonised clock. Rebuilt with the control clock by default; run
  `make figures-published` to regenerate the submitted versions as
  `*_published.png` and compare.

Everything else should match in content. Exact pixel agreement is not expected:
these were rendered at a different size and assembled into composites for
submission.

## Supplementary S6 has no submitted counterpart

The submitted supplementary set runs S1–S5. This repository also builds **S6**
(co-clustering pooled across all mice), which was not part of the submission —
it is included because the pooled distance matrix it draws is shipped in
`data/coclustering/` and supports Figure 1I. Treat it as an addition, not a
reproduction.
