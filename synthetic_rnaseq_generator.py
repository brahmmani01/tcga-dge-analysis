"""
Generates a synthetic RNA-seq count matrix standing in for a real
TCGA download (TCGAbiolinks::GDCdownload requires internet access to
the NCI's GDC portal, which this sandbox doesn't have - see
real_tcga_analysis.R for the actual pipeline to run yourself).

Design: 40 tumor + 40 "normal" synthetic samples, ~2000 genes. A
curated set of ~30 genes are real, well-documented breast cancer
differentially-expressed genes (from published TCGA-BRCA literature:
Perou et al. 2000 intrinsic subtypes; TCGA Network 2012 comprehensive
analysis), assigned real up/down direction and realistic effect
sizes. The remaining ~1970 genes are non-differential background
noise, labeled generically (not fabricated real gene symbols) so
there's no ambiguity about which part of this dataset is literature-
grounded and which is synthetic noise.

Counts are drawn from a negative binomial distribution, matching the
actual statistical structure of real RNA-seq count data (this is the
same distributional assumption DESeq2 itself uses) - this isn't just
a Gaussian stand-in, it's the right distribution family for the
method to be a faithful demo of.
"""
import numpy as np
import pandas as pd

rng = np.random.default_rng(42)

N_TUMOR = 40
N_NORMAL = 40
N_BACKGROUND_GENES = 1970

# Real breast cancer genes with literature-documented direction in
# tumor vs. normal breast tissue (TCGA-BRCA, Perou et al. 2000, TCGA
# Research Network 2012, Nature). Values are illustrative log2FC
# targets consistent with the general literature, not pulled from a
# specific cited table - the point is these are REAL genes with REAL
# known direction, not invented symbols.
CURATED_GENES = {
    # upregulated in tumor - proliferation / oncogenes
    "MKI67":  2.8, "TOP2A": 2.5, "BIRC5": 2.3, "AURKA": 2.1, "CDK1": 2.4,
    "ERBB2":  1.9, "MYC": 1.6, "CCND1": 1.4, "EGFR": 1.2, "MMP9": 1.8,
    "VEGFA":  1.3, "FOXM1": 2.0, "CCNB1": 2.2, "PLK1": 1.9, "RRM2": 1.7,
    # upregulated, hormone/luminal markers (direction varies by subtype
    # in reality; simplified here as modestly up for demo purposes)
    "ESR1":   0.8, "GATA3": 0.6, "FOXA1": 0.7, "PGR": 0.5,
    # downregulated in tumor - normal breast/adipose tissue markers
    "ADIPOQ": -3.2, "FABP4": -3.0, "PPARG": -2.1, "LEP": -2.8,
    "PLIN1":  -2.6, "CIDEC": -2.4, "ADIPOR2": -1.5, "PCK1": -1.8,
    "GPD1":   -2.0, "AQP7": -2.3, "CFD": -1.9, "LPL": -2.2,
}

all_gene_names = list(CURATED_GENES.keys()) + [f"GENE{i:04d}" for i in range(N_BACKGROUND_GENES)]
n_genes = len(all_gene_names)
n_samples = N_TUMOR + N_NORMAL

base_mean = rng.lognormal(mean=5, sigma=1.5, size=n_genes)
dispersion = 0.3 + rng.exponential(0.2, size=n_genes)

log2fc = np.zeros(n_genes)
for i, gene in enumerate(all_gene_names):
    if gene in CURATED_GENES:
        log2fc[i] = CURATED_GENES[gene] + rng.normal(0, 0.15)  # small realistic noise
    else:
        # background: ~97% null, ~3% small random DE (realistic - not
        # every background gene is perfectly null in real data either)
        log2fc[i] = rng.normal(0, 0.1) if rng.random() > 0.03 else rng.normal(0, 0.9)

counts = np.zeros((n_genes, n_samples), dtype=int)
condition = ["tumor"] * N_TUMOR + ["normal"] * N_NORMAL

for i in range(n_genes):
    mu_normal = base_mean[i]
    mu_tumor = base_mean[i] * (2 ** log2fc[i])
    r = 1 / dispersion[i]  # negative binomial size parameter

    p_tumor = r / (r + mu_tumor)
    p_normal = r / (r + mu_normal)

    counts[i, :N_TUMOR] = rng.negative_binomial(r, p_tumor, size=N_TUMOR)
    counts[i, N_TUMOR:] = rng.negative_binomial(r, p_normal, size=N_NORMAL)

sample_names = [f"TCGA-synthetic-T{i+1:02d}" for i in range(N_TUMOR)] + \
               [f"TCGA-synthetic-N{i+1:02d}" for i in range(N_NORMAL)]

counts_df = pd.DataFrame(counts, index=all_gene_names, columns=sample_names)
coldata_df = pd.DataFrame({"sample": sample_names, "condition": condition})

counts_df.to_csv("data/synthetic_counts.csv")
coldata_df.to_csv("data/synthetic_coldata.csv", index=False)

# also save which genes are "ground truth" DE for validating the pipeline
truth_df = pd.DataFrame({
    "gene": all_gene_names,
    "true_log2fc": log2fc,
    "is_curated_real_gene": [g in CURATED_GENES for g in all_gene_names],
})
truth_df.to_csv("data/ground_truth.csv", index=False)

print(f"Generated {n_genes} genes x {n_samples} samples "
      f"({N_TUMOR} tumor, {N_NORMAL} normal)")
print(f"Curated real breast-cancer genes included: {len(CURATED_GENES)}")
print("Saved to data/synthetic_counts.csv, data/synthetic_coldata.csv, data/ground_truth.csv")
