"""
Differential gene expression analysis on the synthetic count matrix -
a Python approximation of what real_tcga_analysis.R's DESeq2 call
does, since DESeq2 (Bioconductor/R) isn't available in this sandbox.

IMPORTANT METHODOLOGICAL NOTE: DESeq2 fits a negative-binomial GLM per
gene with shrinkage-estimated dispersion - the actual field-standard
method (Love, Huber & Anders 2014, Genome Biology). This script uses
Welch's t-test on log2(CPM+1) values with Benjamini-Hochberg FDR
correction instead - simpler, runnable without Bioconductor, and fine
for validating the pipeline logic end-to-end, but it is NOT the same
statistical method and will not give identical numbers to DESeq2 on
the same data. Do not present this script's output as "DESeq2
results" - it's a stand-in that proves the pipeline shape (counts in,
volcano plot + heatmap + pathway hits out), not a substitute for
running the real method in real_tcga_analysis.R.
"""
import numpy as np
import pandas as pd
from scipy import stats

PADJ_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.0


def benjamini_hochberg(pvals):
    """Standard BH FDR correction - same procedure DESeq2 uses for padj."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted = ranked * n / (np.arange(n) + 1)
    # enforce monotonicity (standard BH step)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result = np.empty(n)
    result[order] = adjusted
    return result


def run_dge(counts_path="data/synthetic_counts.csv", coldata_path="data/synthetic_coldata.csv"):
    counts = pd.read_csv(counts_path, index_col=0)
    coldata = pd.read_csv(coldata_path)

    tumor_cols = coldata.loc[coldata.condition == "tumor", "sample"].tolist()
    normal_cols = coldata.loc[coldata.condition == "normal", "sample"].tolist()

    # library-size normalization (CPM) - a simplified stand-in for
    # DESeq2's median-of-ratios size factor normalization
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes, axis=1) * 1e6
    log_cpm = np.log2(cpm + 1)

    # filter low-count genes (same filter the guide applies)
    keep = counts.sum(axis=1) >= 10
    log_cpm = log_cpm[keep]

    tumor_vals = log_cpm[tumor_cols].to_numpy()
    normal_vals = log_cpm[normal_cols].to_numpy()

    t_stat, pval = stats.ttest_ind(tumor_vals, normal_vals, axis=1, equal_var=False)
    log2fc = tumor_vals.mean(axis=1) - normal_vals.mean(axis=1)
    padj = benjamini_hochberg(pval)

    res_df = pd.DataFrame({
        "gene": log_cpm.index,
        "log2FoldChange": log2fc,
        "pvalue": pval,
        "padj": padj,
    }).sort_values("padj")

    res_df["significant"] = (res_df.padj < PADJ_CUTOFF) & (res_df.log2FoldChange.abs() > LOG2FC_CUTOFF)
    return res_df, log_cpm, coldata


if __name__ == "__main__":
    res_df, log_cpm, coldata = run_dge()
    res_df.to_csv("results/dge_results.csv", index=False)

    n_up = ((res_df.significant) & (res_df.log2FoldChange > 0)).sum()
    n_down = ((res_df.significant) & (res_df.log2FoldChange < 0)).sum()
    print(f"Significant genes (padj<{PADJ_CUTOFF}, |log2FC|>{LOG2FC_CUTOFF}): "
          f"{n_up + n_down} total ({n_up} up, {n_down} down)")
    print("\nTop 10 by adjusted p-value:")
    print(res_df.head(10)[["gene", "log2FoldChange", "padj"]].to_string(index=False))

    # validate against ground truth - what fraction of the known
    # curated real genes did this pipeline actually recover?
    truth = pd.read_csv("data/ground_truth.csv")
    curated = truth[truth.is_curated_real_gene]
    merged = res_df.merge(curated, on="gene")
    recovered = merged.significant.sum()
    print(f"\nValidation: recovered {recovered}/{len(curated)} curated real "
          f"breast-cancer genes as significant (sanity check on pipeline logic).")
