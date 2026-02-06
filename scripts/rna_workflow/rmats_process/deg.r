# Packages
library(tximport)
library(tximeta)
library(GenomicFeatures)
library(data.table)
library(dplyr)
library(tibble)
library(DESeq2)
library(DEGreport)
library(ggplot2)
library(org.Hs.eg.db)
library(biomaRt)
library(plyranges)
library(readr)
library(ggtext)
library(formattable)
library(flextable)
library(officer)
library(readxl)
library(openxlsx)

# paths
base_dir <- "/proj/phanstiel_lab/Data/processed/CCSP/rna/seq_pipelines"
rna_dir <- file.path(base_dir,"rna_output")
differential_dir <- file.path(rna_dir,"00_differential_expression")
diff_data_dir <- file.path(differential_dir,"data")
diff_plot_dir <- file.path(differential_dir,"plots")

if(!dir.exists(differential_dir)){dir.create(differential_dir)}
if(!dir.exists(diff_data_dir)){dir.create(diff_data_dir)}
if(!dir.exists(diff_plot_dir)){dir.create(diff_plot_dir)}

# ---- Sample Info ---------------------
sampleInfo <- fread(file.path(base_dir, "samplesheet.txt"))
sampleInfo <- sampleInfo |> mutate(names = paste0(Proj,"_",Condition,"_",Bio_Rep,"_",Tech_Rep)) |> 
                    dplyr::select("names", "Condition","Bio_Rep","Tech_Rep","RIN","Cong(ng/ml)", "DV 200 (for RNA)") |>
                    rename( "DV_200" =  "DV 200 (for RNA)",
                            "Cong_ng_ml" = "Cong(ng/ml)")
sample_nm <- sampleInfo$names

# Add quant path for names
coldata <- sampleInfo 
coldata$files <- file.path(rna_dir, "quant",sample_nm,"quant.sf")
file.exists(coldata$files)

se <- tximeta(coldata)
gse <- summarizeToGene(se)

# Gene level counts
colData(gse)[] <- lapply(colData(gse), as.factor)
colData(gse)$Condition <- factor(colData(gse)$Condition, levels = c("hESC","NPC","Neuron"))
colData(gse)$RIN <- as.numeric(colData(gse)$RIN)
colData(gse)$Cong_ng_ml <- as.numeric(colData(gse)$Cong_ng_ml)
colData(gse)$DV_200 <- as.numeric(colData(gse)$DV_200)

# DeSeq2
dds <- DESeqDataSet(gse, design = ~ Condition)
# expressed genes for 10 FOR at least 34% of samples approx 3 samples
keep <- rowSums(counts(dds) >= 10) >= ceiling(nrow(colData(gse)) * 0.34)
dds_filt <- dds[keep,]

# Fit model
dds_fit <- DESeq(dds_filt, test ="LRT", reduced = ~ 1)

# Get results
res_lrt <- results(dds_fit)
sig_genes_lrt <- res_lrt |> as.data.frame() |> 
    rownames_to_column(var = "gene_id") |> 
    filter(padj < 0.05)

vst_dds <- vst(dds_fit)
rlog_dds <- rlog(dds_fit)
norm_counts <- normTransform(dds_fit)

# PCA 
pca_vst <- prcomp(t(assay(vst_dds)))
pca_df <- as.data.frame(pca_vst$x)
variance_explained <- summary(pca_vst)$importance[2, ]

# Add metadata
pca_df$Condition <- colData(vst_dds)$Condition
pca_df$Tech_Rep <- colData(vst_dds)$Tech_Rep
pca_df$Sample <- rownames(pca_df)

# Condition colors
cond_cols <- c("hESC" = "#419164", "NPC" = "#F9C555", "Neuron" = "#876EC4")

# PCA plot
RNA_PCAPlot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 2, alpha = 0.8) +
    #geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    scale_color_manual(values = cond_cols) +
    labs(
        title = paste0("Gene Expression PCA (n=", nrow(vst_dds), " genes)"),
        x = paste0("PC1: ", round(variance_explained[1] * 100, 2), "%"),
        y = paste0("PC2: ", round(variance_explained[2] * 100, 2), "%")
    ) +
    theme_classic() +
        theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5)
)

save(RNA_PCAPlot, file = file.path(diff_plot_dir, "RNA_PCAPlot.RData"))

#------------------------------------------------------------------------------------------------------
# pairwise DESeq2 results
dds_pair <- DESeq(dds_filt)
res_hesc_npc <- results(dds_pair, contrast = c("Condition", "NPC", "hESC"))
res_npc_neuron <- results(dds_pair, contrast = c("Condition", "Neuron", "NPC"))
res_hesc_neuron <- results(dds_pair, contrast = c("Condition", "Neuron", "hESC"))

# Shrink LFC
resShrink_hesc_npc <- lfcShrink(dds_pair, coef = "Condition_NPC_vs_hESC", type = "apeglm",format =  "GRanges") |>
  plyranges::names_to_column("gene_id")

resShrink_hesc_neuron <- lfcShrink(dds_pair, coef = "Condition_Neuron_vs_hESC", type = "apeglm", format =  "GRanges") |>
  plyranges::names_to_column("gene_id")


# Due to  hesc was the reference that there is no neuron  vs npc coef I will re-level to th NPC 
dds_pair_npc <- dds_filt
dds_pair_npc$Condition <- relevel(dds_pair_npc$Condition, ref = "NPC")
dds_pair_npc <- DESeq(dds_pair_npc)
#resultsNames(dds_pair_npc)

resShrink_npc_neuron <- lfcShrink(dds_pair_npc, coef = "Condition_Neuron_vs_NPC", type = "apeglm", format =  "GRanges") |>
  plyranges::names_to_column("gene_id")

save(dds_fit, dds_filt, dds_pair, dds_pair_npc, file = file.path(diff_data_dir, "deseq2_dds_objects.RData"))
save (resShrink_hesc_npc, resShrink_npc_neuron, resShrink_hesc_neuron,
      file = file.path(diff_data_dir, "pairwise_deseq2_results_shrunken_lfc.RData"))
save(vst_dds, rlog_dds, norm_counts,file = file.path(diff_data_dir, "deseq2_normalized_lrt.RData"))
load(file.path(diff_data_dir, "deseq2_dds_objects.RData"))
load(file.path(diff_data_dir, "pairwise_deseq2_results_shrunken_lfc.RData"))
load(file.path(diff_data_dir, "deseq2_normalized_lrt.RData"))
#------------------------------------------------------------------------------------------------------
# Finding siginificant genes from pairwise comparisons
padj_cutoff <- 0.05
lfc_cutoff <- 1


sig_genes_lrt <- res_lrt |> as.data.frame() |> 
    rownames_to_column(var = "gene_id") |> 
    filter(padj < 0.01)
dim(sig_genes_lrt)

# Cluster as tutorial
top_genes <- sig_genes_lrt |> arrange(padj) |>  pull(gene_id) # 3000 genes out of  15679

cluster_vst <- assay(vst_dds)[top_genes, ]

clusters  <- degPatterns(cluster_vst, 
                         metadata = as.data.frame(colData(vst_dds)),
                         time = "Condition")

save(clusters, file = file.path(diff_data_dir, "degpatterns_clusters.RData"))
load(file.path(diff_data_dir, "degpatterns_clusters.RData"))
