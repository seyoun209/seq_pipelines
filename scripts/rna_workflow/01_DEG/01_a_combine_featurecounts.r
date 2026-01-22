# Combining Featurecounts
  library(data.table)
  library(dplyr)
  library(DESeq2)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(biomaRt)
  library(tibble)

base_dir <- "/work/users/s/e/seyoun/colabs/LCSP/seq_pipelines"
rna_dir <- file.path(base_dir,"rna_output")
differential_dir <- file.path(rna_dir,"00_differential_expression")
diff_data_dir <- file.path(differential_dir,"data")
diff_plot_dir <- file.path(differential_dir,"plots")

if(!dir.exists(differential_dir)){dir.create(differential_dir)}
if(!dir.exists(diff_data_dir)){dir.create(diff_data_dir)}
if(!dir.exists(diff_plot_dir)){dir.create(diff_plot_dir)}

count_files <- list.files(file.path(rna_dir,"featurecounts"), pattern = "_counts\\.txt$", full.names = TRUE)
summary_files <- list.files(file.path(rna_dir,"featurecounts"), pattern = "_counts\\.txt\\.summary$", full.names = TRUE)


# Function to read the counts files ------------------
counts_df <- function(file){
    counts_df <- fread(file)
    sampleID <- sub("_counts\\.txt$","", basename(file))
    counts_col <- names(counts_df)[ncol(counts_df)]
     newcol_df <- counts_df[, .(Geneid, counts = get(counts_col))]
    setnames(newcol_df, "counts", sampleID)
}


#--------------------------------------------------------------------------
counts_li <- lapply(count_files, counts_df)
merged_counts <- Reduce(function(x,y) merge(x,y, by="Geneid"), counts_li)

# Samplesheet
sampleInfo <- fread(file.path(base_dir,"samplesheet.txt"))
DonorInfo <- fread(file.path(base_dir,"donorinfo.txt"))
donorinfo_sub <- DonorInfo |> dplyr::select(Donor, Sex, Age,Date_of_collection)

sampleInfo_sub <- sampleInfo |> mutate(names = paste0(Proj, "_", Donor, "_", Condition, "_", Treatment)) |>
                    dplyr::select(names,Donor, Condition,Treatment,RIN,RNA_Isolation_Date) |> 
                    left_join(donorinfo_sub, by="Donor")

sampleNMs <- sampleInfo_sub$names

sample_norm  <- sampleInfo_sub %>% filter(Condition == "norm")
sample_hypo  <- sampleInfo_sub %>% filter(Condition == "hypo")

sampleNMs_norm <- sample_norm$names
sampleNMs_hypo <- sample_hypo$names

# Subset merged counts for normoxia and hypoxia
merged_counts_norm <- merged_counts[, c("Geneid", intersect(sampleNMs_norm, colnames(merged_counts))), with=FALSE]
merged_counts_hypo <- merged_counts[, c("Geneid", intersect(sampleNMs_hypo, colnames(merged_counts))), with=FALSE]


#--------------------
# Normoxia
# Function to keep which rows
keep_filter <- function(dds) {
  keep <- rowSums(counts(dds) >= 10) >= ceiling(ncol(dds) * 0.4)
  dds[keep, ]
}

normoxiaNMs <- intersect(colnames(merged_counts_norm), sample_norm$names)
counts_norm_dt <- merged_counts_norm[, c("Geneid", normoxiaNMs), with=FALSE]
counts_norm_mat <- counts_norm_dt |> dplyr::select(-Geneid) |> as.data.frame()
rownames(counts_norm_mat) <- counts_norm_dt$Geneid
coldata_norm <- sample_norm |>  filter(names %in% normoxiaNMs) |> arrange(match(names, normoxiaNMs))
rownames(coldata_norm) <- coldata_norm$names

# factors + level
coldata_norm$Treatment <- factor(coldata_norm$Treatment, levels = c("PBS","DOXO"))
coldata_norm$Donor <- factor(coldata_norm$Donor)

dds_norm <- DESeqDataSetFromMatrix(countData = counts_norm_mat, colData = as.data.frame(coldata_norm), 
                    design = ~ Donor + Treatment)

dds_norm_sub <- keep_filter(dds_norm)
dds_norm_fit <- DESeq(dds_norm_sub)

resRNA_norm <- results(dds_norm_fit)
res_Shrink_norm_df <- lfcShrink(dds_norm_fit, coef="Treatment_DOXO_vs_PBS", type="apeglm") |> as.data.frame()
res_Shrink_norm_df$gene_id <- rownames(res_Shrink_norm_df)

save(res_Shrink_norm_df, dds_norm_fit, resRNA_norm, coldata_norm,
     file = file.path(diff_data_dir, "deseq_results_norm_featureCounts.Rdata"))

load(file.path(diff_data_dir, "deseq_results_norm_featureCounts.Rdata")) #res_Shrink_norm_df, dds_norm_fit, resRNA_norm, coldata_norm

res_Shrink_norm_df <- res_Shrink_norm_df[!is.na(res_Shrink_norm_df$padj),]

norm_dds_norm <- normTransform(dds_norm_fit)
rlog_dds_norm <- rlog(dds_norm_fit)
vst_dds_norm <- vst(dds_norm_fit)

save(norm_dds_norm, rlog_dds_norm, vst_dds_norm, file = file.path(diff_data_dir,"deseq_norm_featurecounts.Rdata"))


pca_vst_re <- prcomp(t(assay(vst_dds_norm)))
pca_df <- as.data.frame(pca_vst_re$x)
variance_explained <- summary(pca_vst_re)$importance[2, ]
pca_df$Donor <- colData(vst_dds_norm)$Donor
pca_df$Treatment <- colData(vst_dds_norm)$Treatment
gene_count <- nrow(assay(vst_dds_norm))

normPCAplot <- ggplot(pca_df, aes(x=PC1, y=PC2, label=Donor, color=Treatment)) +
  geom_point(size=2) +
  scale_color_manual(values = c("PBS" = "#0C7BDC", "DOXO" = "#FFC20A")) +
  labs(
    title = paste0("Total genes (n=", gene_count, ") with all samples"),
    x = paste0("PC1 Variance: ", round(variance_explained[1]*100, 2), "%"),
    y = paste0("PC2 Variance: ", round(variance_explained[2]*100, 2), "%")
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        plot.background = element_rect(fill='transparent'))

save(normPCAplot, file = file.path(diff_plot_dir, "PCA_normoxia_deseq_featureCounts.rda"))

# hypoxia ----------------------------------------------------------------------------------------------------
hypoxiaNMs <- intersect(colnames(merged_counts_hypo), sample_hypo$names)
counts_hypo_dt <- merged_counts_hypo[, c("Geneid", hypoxiaNMs), with=FALSE]
counts_hypo_mat <- counts_hypo_dt |> dplyr::select(-Geneid) |> as.data.frame()
rownames(counts_hypo_mat) <- counts_hypo_dt$Geneid
coldata_hypo <- sample_hypo |>  filter(names %in% hypoxiaNMs) |> arrange(match(names, hypoxiaNMs))
rownames(coldata_hypo) <- coldata_hypo$names
# factors + level

coldata_hypo$Treatment <- factor(coldata_hypo$Treatment, levels = c("PBS","DOXO"))
coldata_hypo$Donor <- factor(coldata_hypo$Donor)

dds_hypo <- DESeqDataSetFromMatrix(countData = counts_hypo_mat, colData = as.data.frame(coldata_hypo), 
                    design = ~ Donor + Treatment)

dds_hypo_sub <- keep_filter(dds_hypo)
dds_hypo_fit <- DESeq(dds_hypo_sub)

resRNA_hypo <- results(dds_hypo_fit)
res_Shrink_hypo_df <- lfcShrink(dds_hypo_fit, coef="Treatment_DOXO_vs_PBS", type="apeglm") |> as.data.frame()
res_Shrink_hypo_df$gene_id <- rownames(res_Shrink_hypo_df)

save(res_Shrink_hypo_df, dds_hypo_fit, resRNA_hypo, coldata_hypo,
     file = file.path(diff_data_dir, "deseq_results_hypoxia_featureCounts.Rdata"))
load(file.path(diff_data_dir, "deseq_results_hypoxia_featureCounts.Rdata")) #res_Shrink_hypo_df, dds_hypo_fit, resRNA_hypo, coldata_hypo

res_Shrink_hypo_df <- res_Shrink_hypo_df[!is.na(res_Shrink_hypo_df$padj),]

res_Shrink_hypo_df <- as.data.frame(res_Shrink_hypo_gr)
res_Shrink_hypo_df$gene_id <- rownames(res_Shrink_hypo_df)

norm_dds_hypo <- normTransform(dds_hypo_fit)
rlog_dds_hypo <- rlog(dds_hypo_fit)
vst_dds_hypo <- vst(dds_hypo_fit)
save(norm_dds_hypo, rlog_dds_hypo, vst_dds_hypo, file = file.path(diff_data_dir,"deseq_hypo_featurecounts.Rdata"))


pca_vst_re <- prcomp(t(assay(vst_dds_hypo)))
pca_df <- as.data.frame(pca_vst_re$x)
variance_explained <- summary(pca_vst_re)$importance[2, ]
pca_df$Donor <- colData(vst_dds_hypo)$Donor
pca_df$Treatment <- colData(vst_dds_hypo)$Treatment
gene_count <- nrow(assay(vst_dds_hypo))

hypoPCAplot <- ggplot(pca_df, aes(x=PC1, y=PC2, label=Donor, color=Treatment)) +
  geom_point(size=2) +
  scale_color_manual(values = c("PBS" = "#0C7BDC", "DOXO" = "#FFC20A")) +
  labs(
    title = paste0("Total genes (n=", gene_count, ") with all samples"),
    x = paste0("PC1 Variance: ", round(variance_explained[1]*100, 2), "%"),
    y = paste0("PC2 Variance: ", round(variance_explained[2]*100, 2), "%")
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        plot.background = element_rect(fill='transparent'))

save(hypoPCAplot, file = file.path(diff_plot_dir, "PCA_hypoxia_deseq_featureCounts.rda"))


#------------------------------------------------------------------------------------------------------------------------------------
# Interaction between condition:treatment 
#------------------------------------------------------------------------------------------------------------------------------------

sampleNMs <- intersect(colnames(merged_counts), sampleInfo_sub$names)
merged_counts_dt <- merged_counts[, c("Geneid", sampleNMs), with=FALSE]
counts_dt <- merged_counts_dt |> dplyr::select(-Geneid) |> as.data.frame()
rownames(counts_dt) <- merged_counts_dt$Geneid
coldata <- sampleInfo_sub |>  filter(names %in% sampleNMs) |> arrange(match(names, sampleNMs))

rownames(coldata) <- coldata$names

# factors + level
coldata$Treatment <- factor(coldata$Treatment, levels = c("PBS","DOXO"))
coldata$Donor <- factor(coldata$Donor)
coldata$Condition <- factor(coldata$Condition, levels = c("norm","hypo"))

dds <- DESeqDataSetFromMatrix(countData = counts_dt, colData = as.data.frame(coldata), 
                    design =~ Donor + Condition + Treatment + Condition:Treatment)

dds_sub <- keep_filter(dds)
dds_fit <- DESeq(dds_sub)
res <- results(dds_fit)

res_Shrink_df <- lfcShrink(dds_fit, coef="Conditionhypo.TreatmentDOXO", type="apeglm") |> as.data.frame()
res_Shrink_df$gene_id <- rownames(res_Shrink_df)

save(res_Shrink_df, dds_fit, res, coldata,
     file = file.path(diff_data_dir, "deseq_results_interaction_featureCounts.Rdata"))

load(file.path(diff_data_dir, "deseq_results_interaction_featureCounts.Rdata")) #res_Shrink_df, dds_fit, res, coldata
res_Shrink_df <- res_Shrink_df[!is.na(res_Shrink_df$padj),]

norm_dds <- normTransform(dds_fit)
rlog_dds <- rlog(dds_fit)
vst_dds <- vst(dds_fit)

save(norm_dds, rlog_dds, vst_dds, file = file.path(diff_data_dir,"deseq_interaction_featurecounts.Rdata"))


pca_vst_re <- prcomp(t(assay(vst_dds)))
pca_df <- as.data.frame(pca_vst_re$x)
variance_explained <- summary(pca_vst_re)$importance[2, ]
pca_df$Donor <- colData(vst_dds)$Donor
pca_df$Treatment <- colData(vst_dds)$Treatment
pca_df$Condition <- colData(vst_dds)$Condition
gene_count <- nrow(assay(vst_dds))

interactionPCAplot <- ggplot(pca_df, aes(x=PC1, y=PC2, label=Donor, color=Treatment,shape = Condition)) +
  geom_point(size=2) +
  scale_color_manual(values = c("PBS" = "#0C7BDC", "DOXO" = "#FFC20A")) +
  scale_shape_manual(values = c("norm" = 16, "hypo" = 17)) +
  labs(
    title = paste0("Total genes (n=", gene_count, ") with all samples"),
    x = paste0("PC1 Variance: ", round(variance_explained[1]*100, 2), "%"),
    y = paste0("PC2 Variance: ", round(variance_explained[2]*100, 2), "%")
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        plot.background = element_rect(fill='transparent'))

save(interactionPCAplot, file = file.path(diff_plot_dir, "PCA_interaction_deseq_featureCounts.rda"))

#--------------------------------------------------------------------------------------------------------
# Compare with the tximeta vs featurecounts output from deseq2
#--------------------------------------------------------------------------------------------------------
load(file.path(diff_data_dir, "deseq_results_interaction_featureCounts.Rdata")) #res_Shrink_df, dds_fit, res, coldata
load(file = file.path(diff_data_dir,"deseq_results_interaction_tximeta.Rdata")) # res_Shrink_gr,dds_sub_fit,resRNA_inter, sampleInfo_sub
res_Shrink_inter_gr <- res_Shrink_gr[!is.na(res_Shrink_gr$padj),]
res_Shrink_tximeta_df <- as.data.frame(res_Shrink_inter_gr)

log2fc_compare_df <- res_Shrink_tximeta_df |>
  semi_join(res_Shrink_df, by = "gene_id") |>
  dplyr::select(gene_id, log2FoldChange) |>
  dplyr::rename(log2FC_tximeta = log2FoldChange) |>
  left_join(res_Shrink_df |> dplyr::select(gene_id, log2FoldChange) |> dplyr::rename(log2FC_featureCounts = log2FoldChange),
                by = "gene_id") 

cor_test <- cor.test(log2fc_compare_df$log2FC_tximeta, log2fc_compare_df$log2FC_featureCounts, method = "pearson")
pearson_r  <- cor_test$estimate
pearson_r2 <- pearson_r^2

scatter_plot <- ggplot(log2fc_compare_df,aes(x = log2FC_tximeta, y = log2FC_featureCounts)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  annotate( "text",
    x = min(log2fc_compare_df$log2FC_tximeta),
    y = max(log2fc_compare_df$log2FC_featureCounts),
    hjust = 0,
    label = paste0(
      "Pearson r = ", round(pearson_r, 3),
      "\nR² = ", round(pearson_r2, 3)
    )) +
  theme_classic()



#--------------------------------------------------------------------------------------------------------
# Significant genes 
#--------------------------------------------------------------------------------------------------------
load(file.path(diff_data_dir, "deseq_results_interaction_featureCounts.Rdata")) #res_Shrink_df, dds_fit, res, coldata
load(file.path(diff_data_dir, "deseq_results_norm_featureCounts.Rdata")) #res_Shrink_norm_df, dds_norm_fit, resRNA_norm, coldata_norm
load(file.path(diff_data_dir, "deseq_results_hypoxia_featureCounts.Rdata")) #res_Shrink_hypo_df, dds_hypo_fit, resRNA_hypo, coldata_hypo

# normoxia significant genes -------------------------------------------------------------------------------

diff_norm_sig_p05log2fc1 <- res_Shrink_norm_df |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) # 1774
diff_norm_sig_p01log2fc1 <- res_Shrink_norm_df |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.01) # 1531
diff_norm_sig_p05log2fc15 <- res_Shrink_norm_df |> dplyr::filter(abs(log2FoldChange) > 1.5 & padj < 0.05) # 883


diff_norm_all  <- diff_norm_sig_p05log2fc1
#finiding up sig and down sig
up_norm <- diff_norm_all[diff_norm_all$log2FoldChange > 0 & diff_norm_all$padj < 0.05,]
down_norm <- diff_norm_all[diff_norm_all$log2FoldChange < 0 & diff_norm_all$padj < 0.05,]
static_norm <- res_Shrink_norm_df[!(res_Shrink_norm_df$gene_id %in% c(up_norm$gene_id, down_norm$gene_id)),]


# Making all the background genes and adding the class
res_Shrink_norm_df <- as.data.frame(res_Shrink_norm_gr)
res_Shrink_norm_df$class <- "static"
res_Shrink_norm_df$class[res_Shrink_norm_df$gene_id %in% up_norm$gene_id] <- "gained"
res_Shrink_norm_df$class[res_Shrink_norm_df$gene_id %in% down_norm$gene_id] <- "lost"

# Adding the HGNC symbol
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
res_Shrink_norm_df_cleanVer <- res_Shrink_norm_df |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id))

annots <- AnnotationDbi::select(org.Hs.eg.db, res_Shrink_norm_df_cleanVer$gene_id,
                                columns=c("SYMBOL"), keytype="ENSEMBL")
annots <- dplyr::rename(annots, gene_id = ENSEMBL)
annots_unique <- annots %>% distinct(gene_id, .keep_all = TRUE)
res_Shrink_norm_df_hgnc <- left_join(res_Shrink_norm_df_cleanVer, annots_unique, by = "gene_id") # 2865 of them couldn't find the HGNC

#  Save to run GO and pathway analysis
# up genes
res_Shrink_norm_df_hgnc |> filter(class == "gained") |> dplyr::select(gene_id) |> 
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_normoxia_upregulated_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

# down genes
res_Shrink_norm_df_hgnc |> filter(class == "lost") |> dplyr::select(gene_id) |> 
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_normoxia_downregulated_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

#background genes
res_Shrink_norm_df_hgnc  |> dplyr::select(gene_id) |> 
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_normoxia_background_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

homer_rna <- file.path(base_dir,"scripts/rna_workflow/01_DEG/run_homer.sh")

up_norm_cmd <- paste(
  homer_rna,
  file.path(diff_data_dir, "diff_normoxia_upregulated_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "diff_normoxia_background_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "deseq_homer", "normoxia_up_LFC1_padj05")
)

down_norm_cmd <- paste(
  homer_rna,
  file.path(diff_data_dir, "diff_normoxia_downregulated_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "diff_normoxia_background_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "deseq_homer", "normoxia_down_LFC1_padj05")
)

#system(up_norm_cmd)
#system(down_norm_cmd)
