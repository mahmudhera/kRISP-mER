# fresh_n22 — kRISP-meR manuscript experiment pipeline

Reproducible pipeline for the head-to-head experiments (kRISP-meR vs CRISPOR vs
GuideScan2, plus AlleleAnalyzer), run on rice (aus-N22 vs japonica) and human
(CHM13 vs GRCh38). Scripts are numbered in execution order. Rice uses the bare
number / `R` prefix (japonica = reference build); human uses the `h` suffix.

## Class definitions
Four discordance classes, defined by guide copy-number in the INDIVIDUAL assembly
(c_I) vs the REFERENCE assembly (c_R):

| class | c_I | c_R | truth |
|---|---|---|---|
| present_v_absent | ≥1 | 0 | accept |
| multiple_v_single | ≥2 | 1 | reject |
| ref_present_v_absent | 0 | ≥1 | reject |
| ref_multiple_v_single | 1 | ≥2 | accept |

## Pipeline order

1. **`01`/`R01h` build_exons** — spliced transcripts from the individual (N22 / CHM13)
   assembly + GFF → the one exon namespace for all target windows.
2. **`02`/`R02*` build_targets** — align class bucket guides to exons (`bowtie -v 0`),
   cut 150 bp windows → `targets/<class>/target_N.fasta`. `R02*` builds the
   reference-anchored twin (`ref_targets/`).
3. **`03` verify_targets** — check each window contains its class guide.
4. **`04`/`R04*` class_scoring** — kRISP-meR Expected_cuts per guide per target.
5. **`05`/`R05*` run_krispmer** — both-strand kRISP-meR run (`-n`).
6. **`06`/`R06*` run_crispor**, **`07`/`R07*` run_guidescan** — the other two tools.
7. **`08`/`R08*` headtohead_classN** — accept/reject/not-generated counts per class per tool.
8. **`09`,`10` plot_headtohead** — the combined head-to-head figures.

## Task scripts
- **`T2*` / `plot_*cfd*`** — Task 2: CFD off-target distributions (assembly-CFD `T2c`,
  reads-CFD `T2b`, CRISPOR cfdSpec `T2_cfd_distributions`).
- **`T1_*`** — Task 1: reference-tool-mis-accepted dangerous guides landing in famous /
  trait genes (human + rice).
- **`plot_anchor_comparison_*`** — reference-anchored vs individual-anchored comparison.

## AlleleAnalyzer (AA*) — human only
AlleleAnalyzer was added as a 4th, personalized competitor (human run only; rice not run).
- **`AA1_discriminability_human`** — per-class locus discriminability (loci AA can build ≥1
  allele-specific guide for).
- **`AA2_plot_headtohead_human_with_aa`** — head-to-head with AA as a discriminability bar.
- **`AA3_variant_reality_human`** — ref/alt allele × CHM13/GRCh38 2×2: does the allele AA
  targets actually exist in the individual (CHM13 k-mers)? Handles AA's homopolymer
  "no-guide" sentinels; SNP-only.
- **`AA4_plot_variant_reality`** — the 2×2 figure.
- **`AA5_groundtruth_recall_human`** — THE main AA result. Ground truth = NGG guides in the
  CHM13 assembly but absent from GRCh38, within the same target windows both tools used.
  Recall = of that set, how many each tool designed (matched on 20-mer protospacer).
  kRISP-meR ~97%, AlleleAnalyzer ~0% (pooled, two single-in-reference classes).
- **`AA6_plot_gt_recall`** — the recall figure.

Note: the two `ref_*` classes are excluded from the AA5 recall figure (reference-copy
ambiguity makes the GT window→CHM13 mapping unreliable there); only the two
single-in-reference classes are reported.
