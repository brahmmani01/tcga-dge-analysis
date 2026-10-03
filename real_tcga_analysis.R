# =====================================================================
# REAL TCGA-BRCA differential gene expression analysis
#
# This is the actual pipeline - run this yourself in R/RStudio with an
# internet connection to get REAL results on REAL patient data from
# the NCI's Genomic Data Commons. Everything in dge_analysis.py and
# make_figures.py in this repo is a synthetic stand-in built because
# the sandbox that assembled this repo has no internet access; this
# script is the real thing.
#
# Swap TCGA-BRCA for any of the 33 TCGA project codes to analyze a
# different cancer type (TCGA-LUAD, TCGA-COAD, TCGA-LIHC, TCGA-KIRC,
# etc.)
# =====================================================================

# --- Step 0: install packages (run once) ------------------------------
if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')
BiocManager::install(c('TCGAbiolinks', 'DESeq2', 'EnhancedVolcano',
                        'clusterProfiler', 'org.Hs.eg.db', 'ComplexHeatmap'))
install.packages(c('ggplot2', 'dplyr', 'pheatmap'))

library(TCGAbiolinks)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# --- Step 1: download real TCGA-BRCA RNA-seq data ----------------------
query <- GDCquery(
  project = 'TCGA-BRCA',
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts'
)
GDCdownload(query)
data <- GDCprepare(query)

# --- Step 2: prepare for DESeq2 ----------------------------------------
counts <- assay(data, 'unstranded')
coldata <- colData(data)
coldata$condition <- ifelse(substr(coldata$sample_type, 1, 6) == 'Primar', 'tumour', 'normal')
keep_samples <- !is.na(coldata$condition)
counts <- counts[, keep_samples]
coldata <- coldata[keep_samples, ]

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# --- Step 3: run DESeq2 (the real negative-binomial GLM method) -------
dds <- DESeq(dds)
res <- results(dds, contrast = c('condition', 'tumour', 'normal'))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, 'results/real_dge_results.csv', row.names = FALSE)

cat('Significant genes (padj<0.05, |log2FC|>1):',
    sum(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1), '\n')

# --- Step 4: figures -----------------------------------------------------
# Volcano plot
EnhancedVolcano(res_df, lab = res_df$gene, x = 'log2FoldChange', y = 'padj',
                pCutoff = 0.05, FCcutoff = 1,
                title = 'Tumour vs Normal - TCGA-BRCA',
                subtitle = 'Real GDC data, DESeq2')
ggplot2::ggsave('figures/real_fig1_volcano.png', width = 9, height = 7, dpi = 300)

# Heatmap of top 50
sig_genes <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj) %>% head(50)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[sig_genes$gene, ]
mat_scaled <- t(scale(t(mat)))
annotation <- data.frame(condition = coldata$condition, row.names = colnames(mat))
pheatmap(mat_scaled, annotation_col = annotation, show_rownames = TRUE,
         show_colnames = FALSE, main = 'Top 50 DEGs: Tumour vs Normal',
         filename = 'figures/real_fig2_heatmap.png', width = 10, height = 10)

# --- Step 5: real GO enrichment (proper hypergeometric test) -----------
up_genes <- res_df %>% filter(!is.na(padj), padj < 0.05, log2FoldChange > 1) %>% pull(gene)
gene_ids <- bitr(up_genes, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = 'BP',
                 pAdjustMethod = 'BH', pvalueCutoff = 0.05, readable = TRUE)
dotplot(ego, showCategory = 20, title = 'GO Biological Process Enrichment')
ggplot2::ggsave('figures/real_fig3_go_enrichment.png', width = 9, height = 8, dpi = 300)

write.csv(as.data.frame(ego), 'results/real_go_enrichment.csv', row.names = FALSE)

cat('\nDone. Real results in results/, real figures in figures/.\n')
cat('Compare these against dge_analysis.py + make_figures.py output on the\n')
cat('synthetic data - the pipeline SHAPE should match; the exact numbers will not,\n')
cat('since DESeq2 (negative-binomial GLM) and the synthetic Welch-t-test stand-in\n')
cat('are genuinely different statistical methods.\n')
