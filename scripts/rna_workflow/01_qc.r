library(tidyverse)
library(data.table)
library(stringr)
library(dplyr)
library(DESeq2)
library(DEGreport)
library(colorspace)
library(plotgardener)
library(ggtext)

# QC for RNA-seq
# paths
base_dir <- "/proj/phanstiel_lab/Data/processed/CCSP/rna/seq_pipelines"
rna_dir <- file.path(base_dir,"rna_output")
differential_dir <- file.path(rna_dir,"00_differential_expression")
diff_data_dir <- file.path(differential_dir,"data")
diff_plot_dir <- file.path(differential_dir,"plots")


# align directory
log_dir <- "/proj/phanstiel_lab/Data/processed/CCSP/rna/seq_pipelines/rna_output/align"
log_files <- list.files(log_dir, pattern = "Log.final.out$", full.names = TRUE)

parse_star_log <- function(f) {
  x <- readLines(f)
  tibble(
    sample = str_remove(basename(f), "_Log.final.out"),
    uniquely_mapped_reads = as.numeric(str_extract(x[grep("Uniquely mapped reads number", x)], "[0-9]+")),
    uniquely_mapped_pct = as.numeric(str_extract(x[grep("Uniquely mapped reads %", x)], "[0-9.]+")),
    input_reads = as.numeric(str_extract(x[grep("Number of input reads", x)], "[0-9]+"))
  )
}

star_qc_df <- map_dfr(log_files, parse_star_log) |>
  mutate(
    condition = case_when(
      str_detect(sample, "hESC") ~ "hESC",
      str_detect(sample, "NPC") ~ "NPC",
      str_detect(sample, "Neuron") ~ "Neuron"
    )) |>
  mutate(
    condition = factor(condition, levels = c("hESC","NPC","Neuron")),
    sample = str_replace(sample, "CCSP_", "") |> str_replace("_1_", "_"),
    rep = str_extract(sample, "[0-9]+$"),
    rep_label = paste0("rep", rep)
)

cond_cols <- c(
  "hESC" = "#419164",
  "NPC" = "#F9C555",
  "Neuron" = "#876EC4"
)
# reordering
rep_order <- star_qc_df |>
  distinct(rep_label) |>
  arrange(as.numeric(str_extract(rep_label, "[0-9]+"))) |>
  pull(rep_label)

star_qc_long <- star_qc_df |>
  mutate(rep_label = factor(rep_label, levels = rep_order)) |>
  pivot_longer(cols = c(uniquely_mapped_reads, uniquely_mapped_pct), 
               names_to = "metric", values_to = "value") |>
  mutate(metric = case_when(
    metric == "uniquely_mapped_reads" ~ "Uniquely mapped reads",
    metric == "uniquely_mapped_pct" ~ "Uniquely mapped reads (%)")) 

readsBarplot <- ggplot(star_qc_long, aes(x = rep_label, y = value, fill = condition)) +
  geom_col() +
  labs(y = "Uniquely mapped reads", x = "") +
  scale_y_continuous(breaks = seq(0, 200e6, 25e6), labels = function(x) paste0(x/1e6, "M"), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("hESC" = "#8DBDA2","NPC" = "#FBDC99","Neuron" = "#B7A8DC")) +
  facet_grid(. ~ condition, switch = "x") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(0, 0, 0, 0),
    axis.text.x = element_markdown(size = 6),
    axis.text.y = element_markdown(size = 6),
    axis.title.y = element_markdown(size = 7),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 7)
)

mappedpctBarplot <- ggplot(star_qc_long, aes(x = rep_label, y = value, fill = condition)) +
  geom_col() +
  labs(y = "Uniquely mapped reads (%)", x = "") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("hESC" = "#8DBDA2","NPC" = "#FBDC99","Neuron" = "#B7A8DC")) +
  facet_grid(. ~ condition, switch = "x") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
     plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(0, 0, 0, 0),
     axis.text.x = element_markdown(size = 5),
     axis.text.y = element_markdown(size = 5),
     axis.title.y = element_markdown(size = 7),
     strip.placement = "outside",
     strip.background = element_blank(),
     strip.text.x = element_text(size = 7)
)

save(readsBarplot, file=file.path(diff_plot_dir, "RNA_readsBarplot.rda"))
save(mappedpctBarplot, file=file.path(diff_plot_dir, "RNA_mappedpctBarplot.rda"))






#----------------------------------------------------------------
rmats_dir <- "/proj/phanstiel_lab/Data/processed/CCSP/rna/seq_pipelines/rna_output/rmats/hESC_vs_Neuron"

se_jc <- read_tsv(
  file.path(rmats_dir, "SE.MATS.JC.txt"),
  show_col_types = FALSE
)

parse_psi <- function(x) {
  as.numeric(strsplit(x, ",")[[1]])
}

psi1 <- se_jc$IncLevel1 %>%
  map(parse_psi) %>%
  do.call(rbind, .) %>%
  as.data.frame()

psi2 <- se_jc$IncLevel2 %>%
  map(parse_psi) %>%
  do.call(rbind, .) %>%
  as.data.frame()

colnames(psi1) <- paste0("hESC_", 1:ncol(psi1))
colnames(psi2) <- paste0("Neuron_", 1:ncol(psi2))
psi_mat <- cbind(psi1, psi2)
parse_counts <- function(x) {
  mat <- do.call(rbind, strsplit(as.character(x), ","))
  mat <- apply(mat, 2, as.numeric)
  rowSums(mat, na.rm = TRUE)
}

support1 <- parse_counts(se_jc$IJC_SAMPLE_1) +
            parse_counts(se_jc$SJC_SAMPLE_1)
support2 <- parse_counts(se_jc$IJC_SAMPLE_2) +
            parse_counts(se_jc$SJC_SAMPLE_2)
keep <- support1 >= 20 & support2 >= 20

psi_filt <- psi_mat[keep, ]
psi_filt <- psi_filt[rowMeans(is.na(psi_filt)) < 0.3, ]
cor(psi_filt$hESC_1, psi_filt$hESC_2, use = "pairwise.complete.obs")

parse_psi_vec <- function(x) as.numeric(strsplit(x, ",")[[1]])
  make_psi_long <- function(rmats_file, event_type, group1, group2) {
  df <- read_tsv(rmats_file, show_col_types = FALSE)
  psi1 <- do.call(rbind, map(df$IncLevel1, parse_psi_vec)) |> as.data.frame()
  psi2 <- do.call(rbind, map(df$IncLevel2, parse_psi_vec)) |> as.data.frame()

  colnames(psi1) <- paste0(group1, "_rep", seq_len(ncol(psi1)))
  colnames(psi2) <- paste0(group2, "_rep", seq_len(ncol(psi2)))
  psi_mat <- bind_cols(psi1, psi2)

  psi_long <- psi_mat |> mutate(event_id = row_number()) |>
    pivot_longer(cols = -event_id, names_to = "sample",values_to = "PSI") |>
    separate(sample, into = c("celltype", "rep"), sep = "_") |>
    mutate( rep = factor(rep), celltype = factor(celltype, levels = c("hESC", "NPC", "Neuron")),event = event_type)
  psi_long
}

rmats_base <- "/proj/phanstiel_lab/Data/processed/CCSP/rna/seq_pipelines/rna_output/rmats"

psi_all <- bind_rows(
  make_psi_long(file.path(rmats_base, "hESC_vs_Neuron/SE.MATS.JC.txt"),  "SE",  "hESC", "Neuron"),
  make_psi_long(file.path(rmats_base, "hESC_vs_Neuron/RI.MATS.JC.txt"),  "RI",  "hESC", "Neuron"),
  make_psi_long(file.path(rmats_base, "hESC_vs_Neuron/A3SS.MATS.JC.txt"),"A3SS","hESC", "Neuron"),
  make_psi_long(file.path(rmats_base, "hESC_vs_Neuron/A5SS.MATS.JC.txt"),"A5SS","hESC", "Neuron"),
  make_psi_long(file.path(rmats_base, "hESC_vs_Neuron/MXE.MATS.JC.txt"), "MXE", "hESC", "Neuron")
)

get_group_psi <- function(folder_name) {
  path <- file.path(rmats_base, folder_name, "SE.MATS.JC.txt")
  df <- read_tsv(path, show_col_types = FALSE)
  
  df <- df %>% 
    mutate(event_key = paste(chr, exonStart_0base, exonEnd, strand, sep = "_")) %>%
    select(event_key, IncLevel1, IncLevel2)
  
  psi1 <- do.call(rbind, map(df$IncLevel1, parse_psi_vec)) %>% as.data.frame()
  psi2 <- do.call(rbind, map(df$IncLevel2, parse_psi_vec)) %>% as.data.frame()
  
  groups <- strsplit(folder_name, "_vs_")[[1]]
  colnames(psi1) <- paste0(groups[1], "_rep", 1:ncol(psi1))
  colnames(psi2) <- paste0(groups[2], "_rep", 1:ncol(psi2))
  
  bind_cols(event_key = df$event_key, psi1, psi2)
}

h_vs_neu <- get_group_psi("hESC_vs_Neuron")
h_vs_npc <- get_group_psi("hESC_vs_NPC")

full_psi_mat <- h_vs_neu %>%
  inner_join(h_vs_npc %>% select(event_key, contains("NPC")), by = "event_key") %>%
  distinct(event_key, .keep_all = TRUE)
  column_to_rownames("event_key")

# library(pheatmap)
# cor_mat <- cor(full_psi_mat, method = "spearman", use = "complete.obs")
# pheatmap(cor_mat, 
#          display_numbers = TRUE, 
#          color = colorRampPalette(c("white", "#FBDC99", "#d62728"))(100),
#          main = "Spearman Correlation: Replicate Consistency Check")



# # Load Deseq2 results


# load(file = file.path(diff_data_dir, "deseq2_dds_objects.RData")) # dds_fit, dds_filt, dds_pair, dds_pair_npc
# load(file = file.path(diff_data_dir, "pairwise_deseq2_results_shrunken_lfc.RData")) # resShrink_hesc_npc, resShrink_npc_neuron, resShrink_hesc_neuron
# load(file = file.path(diff_data_dir, "deseq2_normalized_lrt.RData")) # vst_dds, rlog_dds, norm_counts

# # deseqreport

# # Size facor QC
# degCheckFactors(norm_counts)


# Paths
base_dir <- "/proj/phanstiel_lab/Data/processed/CCSP/rna/seq_pipelines"
rmats_dir <- file.path(base_dir, "rna_output","rmats")
rmats_data_dir <- file.path(rmats_dir,"data")
rmats_plots_dir <- file.path(rmats_dir,"plots")

# Combine QC plots 
pdf(file=file.path(rmats_plots_dir, "splicingQC.pdf"),width=8.5,height=6,bg="transparent")

pageCreate(width = 8.5, height = 6, showGuides = FALSE)

plotText("a", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file=file.path(diff_plot_dir, "RNA_mappedpctBarplot.rda")) # mappedpctBarplot
load(file=file.path(diff_plot_dir, "RNA_readsBarplot.rda")) # readsBarplot

plotText("RNA-seq QC uniquely mapping rate", x = 1.7, y = 0.5, just = c("center", "top"), fontfamily = "Helvetica",fontsize = 8)
plotGG(plot=mappedpctBarplot, x=0.35, y=0.65, height=2, width=2.5)

plotText("b", x = 3, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotText("RNA-seq QC total mapped reads", x = 4.75, y = 0.5, just = c("center", "top"), fontfamily = "Helvetica",fontsize = 8)
plotGG(plot=readsBarplot, x=3.35, y=0.65, height=2, width=2.5)

# PCA plot
plotText("c", x = 6, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file.path(rmats_plots_dir, "splicePCAplot.Rdata")) # splicePCAplot         
plotGG(plot=splicePCAplot, x=6.35, y=0.5, height=2, width=2)


# Correlation plot 
plotText("d", x = 0.1, y = 2.8, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load(file.path(rmats_plots_dir, "pearsonBarplot.Rdata")) # pearsonBarplot
plotGG(plot=pearsonBarplot, x=0.35, y= 3.2, height=2.5, width=8.2)

dev.off()
