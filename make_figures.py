"""
Figures 1-3 for the project report: volcano plot, heatmap of top DE
genes, and a pathway-enrichment-style bar chart.

The pathway chart here uses a small, hand-curated gene->pathway
mapping for the specific genes in this demo (real, well-documented
associations - e.g. MKI67/TOP2A/CCNB1/PLK1 -> cell cycle; ADIPOQ/
FABP4/PPARG -> adipocyte/PPAR signaling), NOT a real GO enrichment
statistical test. clusterProfiler + org.Hs.eg.db (used in
real_tcga_analysis.R) run a proper hypergeometric enrichment test
against the full Gene Ontology database - that requires Bioconductor
annotation packages not available in this sandbox. This chart is
illustrative of what the output looks like, not a substitute for
running the real enrichment test.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PADJ_CUTOFF = 0.05
LOG2FC_CUTOFF = 1.0

res_df = pd.read_csv("results/dge_results.csv")
log_cpm = pd.read_csv("data/synthetic_counts.csv", index_col=0)  # reload raw for heatmap z-scoring
coldata = pd.read_csv("data/synthetic_coldata.csv")

# ---------------------------------------------------------------------
# Figure 1: Volcano plot
# ---------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))
neg_log10_padj = -np.log10(res_df.padj.clip(lower=1e-300))

colors = np.where(
    (res_df.padj < PADJ_CUTOFF) & (res_df.log2FoldChange > LOG2FC_CUTOFF), "#C0392B",
    np.where((res_df.padj < PADJ_CUTOFF) & (res_df.log2FoldChange < -LOG2FC_CUTOFF), "#2980B9", "#B0B0B0")
)
ax.scatter(res_df.log2FoldChange, neg_log10_padj, c=colors, s=14, alpha=0.6, edgecolors="none")

# label the curated real genes that are significant
truth = pd.read_csv("data/ground_truth.csv")
curated_sig = res_df.merge(truth[truth.is_curated_real_gene], on="gene")
curated_sig = curated_sig[curated_sig.significant]
for _, row in curated_sig.iterrows():
    y = -np.log10(max(row.padj, 1e-300))
    ax.annotate(row.gene, (row.log2FoldChange, y), fontsize=7.5,
                xytext=(3, 3), textcoords="offset points")

ax.axvline(LOG2FC_CUTOFF, color="gray", linestyle="--", linewidth=0.8)
ax.axvline(-LOG2FC_CUTOFF, color="gray", linestyle="--", linewidth=0.8)
ax.axhline(-np.log10(PADJ_CUTOFF), color="gray", linestyle="--", linewidth=0.8)
ax.set_xlabel("log2 Fold Change (tumor vs normal)")
ax.set_ylabel("-log10(adjusted p-value)")
ax.set_title("Figure 1: Differential Expression, Synthetic TCGA-BRCA-style Data")
fig.tight_layout()
fig.savefig("figures/fig1_volcano_plot.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------
# Figure 2: Heatmap of top 30 significant genes
# ---------------------------------------------------------------------
top_genes = res_df[res_df.significant].reindex(
    res_df[res_df.significant].padj.sort_values().index
).head(30).gene.tolist()

lib_sizes = log_cpm.sum(axis=0)
cpm = log_cpm.div(lib_sizes, axis=1) * 1e6
log2_cpm_top = np.log2(cpm.loc[top_genes] + 1)
z = log2_cpm_top.sub(log2_cpm_top.mean(axis=1), axis=0).div(log2_cpm_top.std(axis=1), axis=0)

# explicit tumor-then-normal order - do NOT rely on alphabetical sort
# ("normal" < "tumor" lexicographically, which silently puts normal
# first and would mislabel the two blocks drawn below)
condition_order = (
    coldata.loc[coldata.condition == "tumor", "sample"].tolist() +
    coldata.loc[coldata.condition == "normal", "sample"].tolist()
)
z = z[condition_order]

fig, ax = plt.subplots(figsize=(10, 8))
im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
ax.set_yticks(range(len(top_genes)))
ax.set_yticklabels(top_genes, fontsize=7)
ax.set_xticks([])
n_tumor = (coldata.condition == "tumor").sum()
ax.axvline(n_tumor - 0.5, color="black", linewidth=1.2)
ax.text(n_tumor / 2, -1.5, "Tumor", ha="center", fontsize=10, fontweight="bold")
ax.text(n_tumor + (len(condition_order) - n_tumor) / 2, -1.5, "Normal", ha="center", fontsize=10, fontweight="bold")
ax.set_title("Figure 2: Top 30 Differentially Expressed Genes (row z-score)")
fig.colorbar(im, ax=ax, shrink=0.6, label="z-score (log2 CPM)")
fig.tight_layout()
fig.savefig("figures/fig2_heatmap.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------
# Figure 3: Pathway "enrichment" bar chart (illustrative, see docstring)
# ---------------------------------------------------------------------
PATHWAY_MAP = {
    "Cell cycle / proliferation": ["MKI67", "TOP2A", "CCNB1", "CDK1", "PLK1", "AURKA", "RRM2", "FOXM1", "BIRC5"],
    "ERBB2/growth factor signaling": ["ERBB2", "EGFR", "MYC", "CCND1"],
    "Angiogenesis / invasion": ["VEGFA", "MMP9"],
    "Hormone receptor signaling": ["ESR1", "GATA3", "FOXA1", "PGR"],
    "Adipocyte / PPAR signaling (down in tumor)": ["ADIPOQ", "FABP4", "PPARG", "PLIN1", "CIDEC", "ADIPOR2", "LEP", "AQP7", "CFD", "LPL"],
    "Lipid/energy metabolism (down in tumor)": ["PCK1", "GPD1"],
}

sig_genes = set(res_df[res_df.significant].gene)
pathway_counts = {p: len(sig_genes & set(genes)) for p, genes in PATHWAY_MAP.items()}
pathway_counts = {k: v for k, v in pathway_counts.items() if v > 0}

fig, ax = plt.subplots(figsize=(8, 4.5))
pathways = list(pathway_counts.keys())
counts = list(pathway_counts.values())
order = np.argsort(counts)
ax.barh([pathways[i] for i in order], [counts[i] for i in order], color="#2E86AB")
ax.set_xlabel("Number of significant genes in pathway (out of curated set)")
ax.set_title("Figure 3: Pathway representation among significant genes\n(illustrative - see script docstring)")
fig.tight_layout()
fig.savefig("figures/fig3_pathway_bars.png", dpi=150)
plt.close(fig)

print("Saved figures/fig1_volcano_plot.png, fig2_heatmap.png, fig3_pathway_bars.png")
