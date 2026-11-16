# figure-regen — manuscript figure scripts

Plotting scripts for the manuscript figures. All emit vector **PDF + PNG**
(zoom-proof), use Figure-2 fonts (Liberation Sans / Arial-compatible) and a
shared soothing palette. Run each from any directory; data paths are hardcoded
to the analysis TSVs.

## Figure 2 — benchmarking (scatter / bar / violin × organism)

- `plot_fig2_panels.py` — the 3×3 benchmarking grid (concordance scatter,
  off-target category bar with CFD retention overlay, global off-target violin)
  for human / mouse / yeast. Emits one PDF+PNG per panel.

## Figure 3 — individual–reference genomic discordance

The discordance figure is assembled in LaTeX from **atomic** panels: one PDF per
(organism × class × block), plus two shared legends.

- `panel_hh.py ORG CLASS` — one head-to-head panel: 3 tools (CRISPOR /
  GuideScan2 / kRISP-meR), 3-category stacked bars (generated&accepted /
  generated&rejected / not-generated). `ORG` = `rice|human`,
  `CLASS` = `multiple_v_single|present_v_absent`. → `hh_<org>_<class>.pdf/.png`
- `panel_cfd.py ORG CLASS` — one assembly-CFD panel (accepted guides): box +
  strip + mean diamond, reference lines at 1 (single-copy) and 2 (duplicated).
  → `cfd_<org>_<class>.pdf/.png`
- `panel_legends.py` — the two shared legends (`legend_hh`, `legend_cfd`) placed
  once in the LaTeX table.

Superseded per-organism combined panels (kept for reference; the atomic panels
above are used in the manuscript):

- `fig_hh_panel.py ORG [rowtag]` — 2 classes × 3 tools head-to-head, one column.
- `fig_cfd_accepted.py ORG [rowtag]` — 2 classes × 3 tools assembly-CFD, one column.

## Data

- Figure 3 CFD panels read `data/asm_cfd_dist_{rice_n22,human_chm13}.tsv`
  (cols: `cfd_cuts class group guide tool truth verdict`).
- Head-to-head counts are hardcoded from the both-strand class-scoring runs
  (`scripts/fresh_n22/*_bothstrand.*`).
