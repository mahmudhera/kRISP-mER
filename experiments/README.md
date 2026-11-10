# kRISP-meR experiments

Scripts used to produce the results in the kRISP-meR manuscript. Code and
documentation only: no sequencing data, assemblies, intermediate files or
result tables are stored here. Paths inside the scripts point at the
directories they were run from on the lab machines and will need adjusting.

Installation, including the Rule Set 2 setup, is in [INSTALL.md](INSTALL.md).

## Layout

| Directory | Contents |
|---|---|
| `pipeline/` | The numbered head-to-head pipeline. Start here. See its own README. |
| `figures/` | Manuscript figure scripts. See its own README. |
| `crispor/` | Batch driver and scraper for running CRISPOR. |
| `kmer-counting/` | Jellyfish counting and lexicographic sorting. |
| `difference-with-ref/` | `hawk`, which splits k-mer counts into the four discordance classes. Source only; compile before use. |
| `analysis/` | Head-to-head comparison and plotting scripts used during development. |

Top level also holds `runBowTie.sh`, `runBowTie_unique.sh`,
`getgRNAs_unique_case.sh` and `getGeneList_multiple_v_single_all.sh`, used by
the k-mer class analysis, and `stageB_full_pipeline.md`, which documents the
assembly-vs-assembly counting and classification procedure end to end.

## The pipeline

`pipeline/` is numbered in execution order and covers three organism pairs:

| Prefix | Individual | Reference |
|---|---|---|
| bare number | rice aus-N22 | rice japonica |
| `R` | rice japonica | (reference build) |
| `h` suffix | human CHM13 | human GRCh38 |

The stages are: build exons, build 150 bp target windows, verify them, score
with each tool, then compare. Steps `05`, `06` and `07` run kRISP-meR, CRISPOR
and GuideScan2 respectively on the same target set, which is what makes the
comparison head to head.

`pipeline/README.md` documents the four discordance classes and each stage.

## How the comparison tools were run

**kRISP-meR**: `05_run_krispmer_*.sh`, both strands (`-n`).

**CRISPOR**: `crispor/run.sh` drives `crispor_batch.py` over a directory of
target FASTAs, with the genome and PAM given on the command line, for example
`--genome "Oryza sativa Japonica" --pam NGG`. `crispor_scrape.py` collects the
results. `pipeline/06_run_crispor_*.sh` wraps this per organism.

**GuideScan2**: `pipeline/07_run_guidescan_*.sh`. Guides are enumerated from
each target FASTA with GuideScan2's own `generate_kmers.py` (`--pam NGG`),
concatenated per class, then scored against a prebuilt genome index with
`guidescan enumerate -m 2 -a NGG --format csv --mode complete`.

The GuideScan2 genome index is built once per reference and reused across
runs; the scripts reference an existing index path rather than rebuilding it.

## Target generation

Targets are not arbitrary windows. The procedure is:

1. Count canonical 23-mers in both assemblies with Jellyfish (`-m 23 -C`) and
   keep those with a PAM (`^CC` or `GG$`).
2. Compare the two count sets with `hawk` to produce four discordance classes,
   defined by copy number in the individual assembly against the reference.
3. Align the guides in each class back to the exon set with `bowtie -v 0`.
4. Cut a 150 bp window centred on each hit. That window is the target.
5. Verify each window still contains the guide that generated it (step `03`).

Step 2 is documented in detail in `stageB_full_pipeline.md`.
