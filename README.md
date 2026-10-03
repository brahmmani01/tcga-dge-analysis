# TCGA-BRCA Differential Gene Expression Analysis

A differential gene expression (DGE) analysis pipeline comparing
breast tumor vs. normal tissue RNA-seq data, built around the method
described in ["The One Biotech Project"](https://x.com/soumalyabhandari)
guide — adapted here for a health informatics / data analyst context
rather than a biotech MSc application, since that's the actual skill
this demonstrates: working with real public genomic data, applying
the field-standard statistical method, and interpreting results
critically.

## Why this belongs in a health informatics portfolio

Precision/genomic medicine is an increasingly central part of health
data work — understanding how RNA-seq differential expression
analysis works, what a p-value/log2FC cutoff actually means, and how
to validate a pipeline against known biology are directly transferable
skills to clinical genomics, pharmacogenomics, and any role touching
molecular data pipelines.

## The honest split: synthetic demo vs. real analysis

**I built this without internet access** (the sandbox that assembled
this repo can't reach the NCI's GDC data portal), so there are two
parallel tracks here, clearly separated:

| | Synthetic demo (runs now) | Real analysis (you run) |
|---|---|---|
| Data | `synthetic_rnaseq_generator.py` — simulated counts, 31 real curated breast-cancer genes with literature-documented direction + ~1970 background genes | `real_tcga_analysis.R` — actual GDC download via `TCGAbiolinks::GDCquery()` |
| Method | `dge_analysis.py` — Welch's t-test + BH-FDR on log2(CPM+1), a simplified stand-in | `real_tcga_analysis.R` — DESeq2's actual negative-binomial GLM (the real field-standard method) |
| Pathway analysis | `make_figures.py` Figure 3 — hand-curated gene→pathway mapping, illustrative only | `real_tcga_analysis.R` — real `clusterProfiler::enrichGO()` hypergeometric test |
| Validates | Pipeline logic: recovers 27/31 known real genes correctly, correct direction 95%+ of the time | Actual TCGA-BRCA biology |

**Do not present the synthetic numbers as real TCGA results.** The
synthetic track exists to prove the pipeline works end-to-end (figures
render, stats are computed correctly, known genes come out on top)
before you spend 10-30 minutes downloading real data. Once you run
`real_tcga_analysis.R`, the `results/real_*.csv` and `figures/real_*.png`
outputs are the genuine article — use those in your actual report.

## A real bug this found, worth knowing about

The first version of `make_figures.py`'s heatmap sorted samples
alphabetically ("normal" < "tumor"), then labeled the first block
"Tumor" — silently mislabeling every column. Caught by checking the
known direction of TOP2A/FABP4 numerically against the heatmap colors
rather than trusting the plot visually. Fixed by sorting explicitly
instead of relying on alphabetical order. Worth remembering: a
figure that *looks* right isn't the same as being right — check
against a number you already know.

## Running it

**Synthetic demo (works right now, no setup beyond Python):**
```bash
pip install numpy pandas scipy matplotlib
python3 synthetic_rnaseq_generator.py
python3 dge_analysis.py
python3 make_figures.py
python3 tests/test_dge_analysis.py
```

**Real analysis (requires R, RStudio, internet, ~30-60 min runtime):**
1. Install R (r-project.org) and RStudio (posit.co/download/rstudio-desktop)
2. Open `real_tcga_analysis.R` in RStudio, run it top to bottom
3. Real results land in `results/real_dge_results.csv` and
   `figures/real_fig*.png`

## Files

```
synthetic_rnaseq_generator.py   synthetic count matrix + ground truth
dge_analysis.py                 Welch's t-test + BH-FDR (demo method)
make_figures.py                 volcano plot, heatmap, pathway chart
real_tcga_analysis.R             THE REAL PIPELINE - run this for genuine results
data/                           synthetic counts, coldata, ground truth
results/                        DGE results tables (synthetic + real, once run)
figures/                        all figures (synthetic + real, once run)
report/project_report_template.md   filled-in example report structure
tests/test_dge_analysis.py      4 tests validating the demo pipeline
```

## What would make this a genuinely strong portfolio piece

- Run `real_tcga_analysis.R` and replace every synthetic figure/number
  with the real output before showing this to anyone
- Write the actual project report (`report/` has a filled-in template
  using synthetic numbers as a placeholder — replace with real ones)
- Pick the cancer type that connects to something you actually care
  about, not just the default (breast cancer here was chosen because
  it's the best-documented example in the literature, not because
  it's the most relevant to your interests)
- Consider Extension 1 from the original guide (survival analysis via
  TCGA clinical data) if you want to add a translational/outcomes
  dimension, which connects more directly to health informatics than
  pure molecular biology does

## Attribution

Method and project structure adapted from "The One Biotech Project
That Strengthens Any University Application" by Soumalya Bhandari
(@soumalyabhandari, MSc Biotechnology, University of Glasgow). Curated
gene list and direction based on general TCGA-BRCA literature (Perou
et al. 2000, *Nature*; TCGA Research Network 2012, *Nature*) —
illustrative of known biology, not quoted from a specific results
table.
