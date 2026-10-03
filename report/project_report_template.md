# Differential Gene Expression Analysis of Breast Tumor vs. Normal Tissue

**⚠️ EXAMPLE DRAFT — filled in with SYNTHETIC data to show the structure.
Replace every number/figure below with output from `real_tcga_analysis.R`
before using this anywhere real.**

## 1. Background

Breast cancer remains the most commonly diagnosed cancer in women
worldwide. Prior work using TCGA data has established molecular
subtypes with distinct expression signatures (Perou et al. 2000) and
catalogued the mutational and expression landscape across over a
thousand tumors (TCGA Research Network, 2012). This analysis asks
which genes are significantly differentially expressed between tumor
and adjacent normal breast tissue, and what biological processes those
genes implicate.

*(Real version: cite 3-5 actual papers from PubMed relevant to your
chosen cancer type.)*

## 2. Methods

Data: [SYNTHETIC — replace with: TCGA-BRCA RNA-seq, STAR-counts
workflow, downloaded via TCGAbiolinks from the NCI GDC portal, N=XXX
tumor / XXX normal samples]. Differential expression: [SYNTHETIC —
replace with: DESeq2 negative-binomial GLM]. Significance threshold:
adjusted p < 0.05 (Benjamini-Hochberg), |log2 fold change| > 1. Pathway
analysis: [SYNTHETIC — replace with: clusterProfiler GO Biological
Process enrichment, BH-corrected].

## 3. Results

*(Example numbers from the synthetic demo — replace with real output)*

42 genes met the significance threshold (23 upregulated, 19
downregulated in tumor vs. normal). The most significantly altered
gene was TOP2A (log2FC = 2.94, padj = 5.5×10⁻²¹), consistent with its
known role as a proliferation marker. The most significantly
downregulated gene was ADIPOQ (log2FC = -3.14), consistent with the
loss of normal adipose tissue signal in tumor samples. [Insert Figure
1: volcano plot] [Insert Figure 2: heatmap of top 30 genes]

Pathway analysis of upregulated genes showed strong representation of
cell-cycle/proliferation genes (MKI67, TOP2A, CCNB1, CDK1, PLK1,
AURKA) and ERBB2/growth-factor signaling genes (ERBB2, EGFR, MYC).
Downregulated genes were dominated by adipocyte/PPAR signaling genes
(ADIPOQ, FABP4, PPARG), consistent with the known histological
difference between fat-rich normal breast tissue and cellular tumor
tissue. [Insert Figure 3: pathway chart]

## 4. Discussion

*(This section is the one to spend the most time on — see the guide's
note that this is what separates a strong applicant from an average
one.)*

The recovery of well-established proliferation markers (MKI67, TOP2A)
and the known ERBB2 signaling axis aligns with decades of breast
cancer molecular biology — a sanity check that the pipeline is working
correctly, not a novel finding on its own. [Real version: discuss
whether your actual results recovered expected genes, and name any
unexpected findings.]

**Limitations**: [SYNTHETIC demo used a Welch's t-test rather than
DESeq2's negative-binomial model — a real limitation of the demo
pipeline, not of the real analysis once run]. A real analysis should
discuss: sample size and tumor/normal imbalance, potential batch
effects across TCGA sequencing centers, and that differential
expression alone does not establish causality in cancer biology.

**Next steps**: survival analysis to test whether expression of the
top genes correlates with patient outcomes (see Extension 1 in the
source guide); comparison against a second cancer type to distinguish
shared vs. cancer-type-specific dysregulation.

## 5. References

[Minimum 5 peer-reviewed papers — e.g.:]
1. Perou CM, et al. Molecular portraits of human breast tumours.
   *Nature*. 2000;406:747-752.
2. The Cancer Genome Atlas Network. Comprehensive molecular portraits
   of human breast tumours. *Nature*. 2012;490:61-70.
3. Love MI, Huber W, Anders S. Moderated estimation of fold change and
   dispersion for RNA-seq data with DESeq2. *Genome Biol*.
   2014;15:550.
