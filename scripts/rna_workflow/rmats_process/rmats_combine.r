# Library packages
library(maser)
library(tidyverse)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(data.table)
library(stringr)
library(dplyr)
library(limma)
library(ggvenn)
library(plotgardener)
library(openxlsx)
library(readr)
library(eulerr)

source("/proj/phanstiel_lab/Data/processed/CCSP/rna/seq_pipelines/scripts/rna_workflow/rmats_process/utils.r")
# Paths
base_dir <- "/proj/phanstiel_lab/Data/processed/CCSP/rna/seq_pipelines"
rmats_dir <- file.path(base_dir, "rna_output","rmats")
rmats_data_dir <- file.path(rmats_dir,"data")
rmats_plots_dir <- file.path(rmats_dir,"plots")
rmats_homer <- file.path(rmats_dir,"homer")

if(!dir.exists(rmats_data_dir)) {dir.create(rmats_data_dir, recursive = TRUE)}
if(!dir.exists(rmats_plots_dir)) {dir.create(rmats_plots_dir, recursive = TRUE)}
if(!dir.exists(rmats_homer)) {dir.create(rmats_homer, recursive = TRUE)}

hesc_npc_dir <- file.path(rmats_dir,"hESC_vs_NPC")
hesc_neuron_dir <- file.path(rmats_dir,"hESC_vs_Neuron")
npc_neuron_dir <- file.path(rmats_dir,"NPC_vs_Neuron")
chrs <- paste0("chr",c(1:22,"X"))

#==============================================================================================
# Junction Only first
#==============================================================================================

# hESC vs NPC----------------------------------------------------------------------------------
hesc_npc_jc <- maser(hesc_npc_dir, c("hESC", "NPC"), ftype = "JC")
hesc_npc_jcec <- maser(hesc_npc_dir, c("hESC", "NPC"), ftype = "JCEC")
save(hesc_npc_jc,file = file.path(rmats_data_dir, "hesc_npc_jc.Rdata"))
save(hesc_npc_jcec,file = file.path(rmats_data_dir, "hesc_npc_jcec.Rdata"))

#load(file.path(rmats_data_dir, "hesc_npc_jc.Rdata"))
head(summary(hesc_npc_jc, type = "SE")[, 1:8])
hesc_npc_filt <- filterByCoverage(hesc_npc_jc, avg_reads = 10)
hesc_npc_filt <- filterByChromosome(hesc_npc_filt, chromosomes = chrs)
sig_hesc_npc <- topEvents(hesc_npc_filt, fdr = 0.05, deltaPSI = 0.1)

# summary(sig_hesc_npc,'SE') |> dim()
# volcano(hesc_npc_filt, fdr = 0.05, deltaPSI = 0.1, "SE")
# splicingDistribution(hesc_npc_filt)

hesc_npc_sig_genes <- c(summary(sig_hesc_npc,'SE')[,"GeneID"],
                    summary(sig_hesc_npc,'RI')[,"GeneID"],
                    summary(sig_hesc_npc,'A3SS')[,"GeneID"],
                    summary(sig_hesc_npc,'A5SS')[,"GeneID"],
                    summary(sig_hesc_npc,'MXE')[,"GeneID"]) |> unique()

hesc_npc_background_genes <- c(summary(hesc_npc_filt,'SE')[,"GeneID"],
                    summary(hesc_npc_filt,'RI')[,"GeneID"],
                    summary(hesc_npc_filt,'A3SS')[,"GeneID"],
                    summary(hesc_npc_filt,'A5SS')[,"GeneID"],
                    summary(hesc_npc_filt,'MXE')[,"GeneID"]) |> unique() 
                    
hesc_npc_background_genes <- gsub("\\..*$", "",hesc_npc_background_genes)

# npc vs Neuron---------------------------------------------------------------------------------    
npc_neuron_jc <- maser(npc_neuron_dir, c("NPC", "Neuron"), ftype = "JC")
npc_neuron_jcec <- maser(npc_neuron_dir, c("NPC", "Neuron"), ftype = "JCEC")
save(npc_neuron_jc,file = file.path(rmats_data_dir, "npc_neuron_jc.Rdata"))
save(npc_neuron_jcec,file = file.path(rmats_data_dir, "npc_neuron_jcec.Rdata"))

npc_neuron_filt <- filterByCoverage(npc_neuron_jc, avg_reads = 10)
npc_neuron_filt <- filterByChromosome(npc_neuron_filt, chromosomes = chrs)
sig_npc_neuron <- topEvents(npc_neuron_filt, fdr = 0.05, deltaPSI = 0.1)
npc_neuron_sig_genes <- c(summary(sig_npc_neuron,'SE')[,"GeneID"],
                        summary(sig_npc_neuron,'RI')[,"GeneID"],
                        summary(sig_npc_neuron,'A3SS')[,"GeneID"],
                        summary(sig_npc_neuron,'A5SS')[,"GeneID"],
                        summary(sig_npc_neuron,'MXE')[,"GeneID"]) |> unique()

npc_neuron_background_genes <- c(summary(npc_neuron_filt,'SE')[,"GeneID"],
                    summary(npc_neuron_filt,'RI')[,"GeneID"],
                    summary(npc_neuron_filt,'A3SS')[,"GeneID"],
                    summary(npc_neuron_filt,'A5SS')[,"GeneID"],
                    summary(npc_neuron_filt,'MXE')[,"GeneID"]) |> unique() 
                    
npc_neuron_background_genes <- gsub("\\..*$", "",npc_neuron_background_genes)

# hESC vs Neuron--------------------------------------------------------------------------------
hesc_neuron_jc <- maser(hesc_neuron_dir, c("hESC", "Neuron"), ftype = "JC")
hesc_neuron_jcec <- maser(hesc_neuron_dir, c("hESC", "Neuron"), ftype = "JCEC")
save(hesc_neuron_jc,file = file.path(rmats_data_dir, "hesc_neuron_jc.Rdata"))
save(hesc_neuron_jcec,file = file.path(rmats_data_dir, "hesc_neuron_jcec.Rdata"))

hesc_neuron_filt <- filterByCoverage(hesc_neuron_jc, avg_reads = 10)
hesc_neuron_filt <- filterByChromosome(hesc_neuron_filt, chromosomes = chrs)
sig_hesc_neuron <- topEvents(hesc_neuron_filt, fdr = 0.05, deltaPSI = 0.1)
hesc_neuron_sig_genes <- c(summary(sig_hesc_neuron,'SE')[,"GeneID"],
                        summary(sig_hesc_neuron,'RI')[,"GeneID"],
                        summary(sig_hesc_neuron,'A3SS')[,"GeneID"],
                        summary(sig_hesc_neuron,'A5SS')[,"GeneID"],
                        summary(sig_hesc_neuron,'MXE')[,"GeneID"]) |> unique()

hesc_neuron_background_genes <- c(summary(hesc_neuron_filt,'SE')[,"GeneID"],
                    summary(hesc_neuron_filt,'RI')[,"GeneID"],
                    summary(hesc_neuron_filt,'A3SS')[,"GeneID"],
                    summary(hesc_neuron_filt,'A5SS')[,"GeneID"],
                    summary(hesc_neuron_filt,'MXE')[,"GeneID"]) |> unique() 
                    
hesc_neuron_background_genes <- gsub("\\..*$", "",hesc_neuron_background_genes)


# Save significant splicing genes in excel
supp_table <- file.path(rmats_dir,"supplementary_tables")
if(!dir.exists(supp_table)) dir.create(supp_table)

comparisons <- list(
  "hESC_vs_NPC" = sig_hesc_npc,
  "NPC_vs_Neuron" = sig_npc_neuron,
  "hESC_vs_Neuron" = sig_hesc_neuron
)

event_types <- c("SE", "RI", "A3SS", "A5SS", "MXE")

for (comp_name in names(comparisons)) {
  maser_obj <- comparisons[[comp_name]]
  groups <- strsplit(comp_name, "_vs_")[[1]]
  group1 <- groups[1]
  group2 <- groups[2]

  wb <- createWorkbook()
  
  for (event in event_types) {
    event_data <- summary(maser_obj, type = event)
    if (nrow(event_data) > 0) {
      # change column name
      colnames(event_data)[colnames(event_data) == "PSI_1"] <- paste0("PSI_", group1)
      colnames(event_data)[colnames(event_data) == "PSI_2"] <- paste0("PSI_", group2)
      # Make tab
      addWorksheet(wb, event)
      writeData(wb, sheet = event, x = event_data)
    }
  }
  output_file <- file.path(supp_table, paste0("Sig_diff_AS_fdr_05_deltaPSI_01", comp_name, ".xlsx"))
  saveWorkbook(wb, output_file, overwrite = TRUE)
}


# All background genes
all_background_genes <- c(hesc_neuron_background_genes, npc_neuron_background_genes, hesc_npc_background_genes) |> unique()
write.table(all_background_genes, file.path(rmats_data_dir, "all_rmats_background_genes.txt"),
             sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE)

# Significant genes
# hESC vs NPC
gsub("\\..*$", "", hesc_npc_sig_genes) |> write.table(file.path(rmats_data_dir, "hesc_npc_sig_genes_psi10fdr05.txt"),
             sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE)

# NPC vs Neuron
gsub("\\..*$", "", npc_neuron_sig_genes) |> write.table(file.path(rmats_data_dir, "npc_neuron_sig_genes_psi10fdr05.txt"),
             sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE)

# hESC vs Neuron
gsub("\\..*$", "", hesc_neuron_sig_genes) |> write.table(file.path(rmats_data_dir, "hesc_neuron_sig_genes_psi10fdr05.txt"),
             sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE)


save(hesc_npc_filt,npc_neuron_filt,hesc_neuron_filt,
     sig_hesc_npc,sig_npc_neuron,sig_hesc_neuron,
     hesc_npc_sig_genes,npc_neuron_sig_genes,hesc_neuron_sig_genes,
     file = file.path(rmats_data_dir, "rmats_jc_filtered_results.Rdata"))

load(file.path(rmats_data_dir, "rmats_jc_filtered_results.Rdata"))

#==============================================================================================
# Venn diagram for events
#==============================================================================================
# Function to create unique event IDs
create_event_id <- function(maser_obj) {
    event_types <- c("SE", "RI", "MXE", "A3SS", "A5SS")
    all_events <- map_dfr(event_types, function(type) {
        df <- summary(maser_obj, type)
        if (nrow(df) == 0) return(NULL)
        df <- df |> filter(Chr %in% chrs)
        # Get coordinate columns
        coord_cols <- switch(type,
            "SE"   = c("exon_target","exon_upstream","exon_downstream"),
            "MXE"  = c("exon_1", "exon_2","exon_upstream","exon_downstream"),
            "RI"   = c("exon_ir","exon_upstream","exon_downstream"),
            "A3SS" = c("exon_long", "exon_short","exon_flanking"),
            "A5SS" = c("exon_long", "exon_short","exon_flanking")
        )
        
        # Create unique event ID
        df <- df |> unite("coords", all_of(coord_cols), sep = "_", remove = FALSE) |>
            mutate(event_id = paste(type, Chr, Strand, coords, sep = "_"),EventType = type) |>
            select(event_id, EventType, GeneID,geneSymbol)
        return(df)
    })
    return(all_events)
}

# Extract events from each comparison
hesc_npc_events <- create_event_id(sig_hesc_npc)
npc_neuron_events <- create_event_id(sig_npc_neuron)
hesc_neuron_events <- create_event_id(sig_hesc_neuron)

# events
venn_events_list <- list(
    "hESC vs. NPC" = hesc_npc_events$event_id,
    "NPC vs. Neuron" = npc_neuron_events$event_id,
    "hESC vs. Neuron" = hesc_neuron_events$event_id
)


eventsVenndiagram <- ggvenn(venn_events_list, fill_color = c("hESC vs. NPC" = "#F8A27E","NPC vs. Neuron" = "#FBDC99","hESC vs. Neuron"= "#5BB5D5"),
text_color = "black",
text_size = 2.5,
set_name_size = 2.5,
set_name_color = c("hESC vs. NPC" = "#ee6c34","NPC vs. Neuron" = "#f5bb3f","hESC vs. Neuron"= "#23aadb"),
show_percentage = FALSE,
stroke_size = 0,padding=0) +
  theme_void() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),    
        axis.line = element_blank(),
        legend.position = "none",
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0, 0, 0, 0))

save(eventsVenndiagram, file = file.path(rmats_plots_dir, "eventsVenndiagram.rda"))


geneVenndiagram <- ggvenn(venn_gene_list, fill_color = c("hESC vs. NPC" = "#F8A27E","NPC vs. Neuron" = "#FBDC99","hESC vs. Neuron"= "#5BB5D5"),
text_color = "black",
text_size = 2.5,
set_name_size = 2.5,
set_name_color = c("hESC vs. NPC" = "#ee6c34","NPC vs. Neuron" = "#f5bb3f","hESC vs. Neuron"= "#23aadb"),
show_percentage = FALSE,
stroke_size = 0,padding=0) +
  theme_void() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),    
        axis.line = element_blank(),
        legend.position = "none",
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0, 0, 0, 0))
save(geneVenndiagram, file = file.path(rmats_plots_dir, "geneVenndiagram.rda"))

# barplot to counts each event type 

summarise_event_gene <- function(event_df) {
    event_df |> group_by(EventType) |> 
    summarise(
      n_events = n_distinct(event_id),
      n_genes  = n_distinct(GeneID),
      .groups = "drop"
    )
}

hesc_npc_summary <- summarise_event_gene(hesc_npc_events)
npc_neuron_summary <- summarise_event_gene(npc_neuron_events)
hesc_neuron_summary <- summarise_event_gene(hesc_neuron_events)

event_order <- c("SE", "RI", "A3SS", "A5SS", "MXE")

hesc_npc_summary$EventType <- factor(hesc_npc_summary$EventType, levels = event_order)
npc_neuron_summary$EventType <- factor(npc_neuron_summary$EventType, levels = event_order)
hesc_neuron_summary$EventType <- factor(hesc_neuron_summary$EventType, levels = event_order)

hesc_npc_summary$comparison <- "hESC vs. NPC"
npc_neuron_summary$comparison <- "NPC vs. Neuron"
hesc_neuron_summary$comparison <- "hESC vs. Neuron"

barplot_summary_all <- bind_rows(hesc_npc_summary,npc_neuron_summary,hesc_neuron_summary)
#barplot_event_summary <- barplot_summary_all |> select(EventType, n_events, comparison)
#barplot_gene_summary <- barplot_summary_all |> select(EventType, n_genes, comparison)

barplot_summary_long <- barplot_summary_all |> pivot_longer(cols = c(n_events, n_genes),names_to = "level",values_to = "count") |> 
    mutate(level = recode(level, n_events = "Differential splicing (event types)",n_genes = "Differential splicing (genes)"),
            EventType = factor(EventType, levels = c("SE", "RI", "A3SS", "A5SS", "MXE")),
            comparison = factor(comparison, levels = c("hESC vs. NPC", "NPC vs. Neuron", "hESC vs. Neuron")))

n_annotation <- barplot_summary_long |> group_by(level) |> 
            summarise(n = sum(count),.groups = "drop") |> 
            mutate(label = paste0("n = ", n), x = Inf,y = -Inf)

comp_cols <- c("hESC vs. NPC" = "#F8A27E","NPC vs. Neuron" = "#FBDC99","hESC vs. Neuron"= "#5BB5D5")

countAS_barPlot <- function(df, level_name) {
  ggplot( df |> dplyr::filter(level == level_name),aes(x = EventType, y = count, fill = comparison)) +
        geom_col(position = position_dodge(width = 0.7),width = 0.65) +
        geom_text(aes(label = count),position = position_dodge(width = 0.7),vjust = -0.3,size = 2,color = "black") +
        scale_fill_manual(values = comp_cols) +
        scale_y_continuous(limits = c(0, NA),expand = expansion(mult = c(0, 0.15))) +
        labs(x = "Splicing event type",y = level_name) +
        theme_classic() +
        theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = c(0.9, 0.9),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0, 0, 0, 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 5.5),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(0.25, "cm"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7)
    )
}

AS_eventsBarplot <- countAS_barPlot(barplot_summary_long, "Differential splicing (event types)")
AS_genesBarplot <- countAS_barPlot(barplot_summary_long, "Differential splicing (genes)")

save(AS_eventsBarplot, file = file.path(rmats_plots_dir, "AS_eventsBarplot.rda"))
save(AS_genesBarplot, file = file.path(rmats_plots_dir, "AS_genesBarplot.rda"))

#==============================================================================================
# Pathway analysis for the pairwise comparison
homer_script <- file.path(base_dir,"scripts/rna_workflow/rmats_process/run_homer.sh")

hesc_npc_cmd <- paste(
  homer_script,
  file.path(rmats_data_dir, "hesc_npc_sig_genes_psi10fdr05.txt"),
  file.path(rmats_data_dir, "all_rmats_background_genes.txt"),
  file.path(rmats_homer, "hesc_npc_padj05psi10")
)
system(hesc_npc_cmd)

npc_neuron_cmd <- paste(
  homer_script,
  file.path(rmats_data_dir, "npc_neuron_sig_genes_psi10fdr05.txt"),
  file.path(rmats_data_dir, "all_rmats_background_genes.txt"),
  file.path(rmats_homer, "npc_neuron_padj05psi10")
)
system(npc_neuron_cmd)

hesc_neuron_cmd <- paste(
  homer_script,
  file.path(rmats_data_dir, "hesc_neuron_sig_genes_psi10fdr05.txt"),
  file.path(rmats_data_dir, "all_rmats_background_genes.txt"),
  file.path(rmats_homer, "hesc_neuron_padj05psi10")
)
system(hesc_neuron_cmd)


# Subset only specific to the pairwise comparison genes
hesc_npc_specific <- setdiff(hesc_npc_sig_genes, c(npc_neuron_sig_genes, hesc_neuron_sig_genes))
gsub("\\..*$", "", hesc_npc_specific) |> write.table(file.path(rmats_data_dir, "hesc_npc_specific_psi10fdr05.txt"),
             sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE)

npc_neuron_specific <- setdiff(npc_neuron_sig_genes, c(hesc_npc_sig_genes, hesc_neuron_sig_genes))
gsub("\\..*$", "", npc_neuron_specific) |> write.table(file.path(rmats_data_dir, "npc_neuron_specific_psi10fdr05.txt"),
             sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE)

hesc_neuron_specific <- setdiff(hesc_neuron_sig_genes, c(hesc_npc_sig_genes, npc_neuron_sig_genes))
gsub("\\..*$", "", hesc_neuron_specific) |> write.table(file.path(rmats_data_dir, "hesc_neuron_specific_psi10fdr05.txt"),
             sep="\t",quote = FALSE, row.names = FALSE, col.names = FALSE)

hesc_npc_specific_cmd <- paste(
  homer_script,
  file.path(rmats_data_dir, "hesc_npc_specific_psi10fdr05.txt"),
  file.path(rmats_data_dir, "all_rmats_background_genes.txt"),
  file.path(rmats_homer, "hesc_npc_specific")
)
system(hesc_npc_specific_cmd)


npc_neuron_specific_cmd <- paste(
  homer_script,
  file.path(rmats_data_dir, "npc_neuron_specific_psi10fdr05.txt"),
  file.path(rmats_data_dir, "all_rmats_background_genes.txt"),
  file.path(rmats_homer, "npc_neuron_specific")
)
system(npc_neuron_specific_cmd)

hesc_neuron_specific_cmd <- paste(
  homer_script,
  file.path(rmats_data_dir, "hesc_neuron_specific_psi10fdr05.txt"),
  file.path(rmats_data_dir, "all_rmats_background_genes.txt"),
  file.path(rmats_homer, "hesc_neuron_specific")
)
system(hesc_neuron_specific_cmd)


# Pathway plot --------------------------------------------------------------------------
homer_comparisons <- list(
  "hESC_vs_NPC" = file.path(rmats_homer, "hesc_npc_specific"),
  "NPC_vs_Neuron" = file.path(rmats_homer, "npc_neuron_specific"),
  "hESC_vs_Neuron" = file.path(rmats_homer, "hesc_neuron_specific")
)

pathway_types <- c("kegg.txt")

homer_kegg_list <- list()

for (comp_name in names(homer_comparisons)) {
  homer_dir <- homer_comparisons[[comp_name]]
  homer_kegg_list[[comp_name]] <- list()
  for (pathway_file in pathway_types) {
    input_path <- file.path(homer_dir, pathway_file)
    
    if (file.exists(input_path)) {
      pathway_data <- read_delim(input_path, delim = "\t", show_col_types = FALSE) |>
        mutate(pval = exp(1)^logP) |>
        filter(pval < 0.05) |>  # p-value cutoff
        distinct(Term, .keep_all = TRUE) |>
        mutate(`-log10pval` = -log10(pval))
      
      # fomatting
      pathway_table <- pathway_data |>
        dplyr::select(-logP, -pval) |>
        relocate(`-log10pval`, .after = Enrichment) |>
        arrange(desc(`-log10pval`))
      
       homer_kegg_list[[comp_name]] <- pathway_table

      pathway_name <- gsub(".txt", "", pathway_file)
      output_file <- file.path(rmats_homer, "tables", paste0("specific_", pathway_name, "_", comp_name, ".csv"))
      if (!dir.exists(file.path(rmats_homer,"tables"))) { dir.create(file.path(rmats_homer,"tables"), recursive = TRUE)}
      write_csv(pathway_table, file = output_file)
    }
  }
}

select_kegg_term <- function(df, category_name) {
    df |> as_tibble() |>
        arrange(desc(`-log10pval`)) |> head(5) |> mutate(category = category_name)
}

# combined kegg Data
kegg_plotting_data <- bind_rows(
  select_kegg_term(homer_kegg_list$hESC_vs_NPC, "hESC vs. NPC"),
  select_kegg_term(homer_kegg_list$NPC_vs_Neuron, "NPC vs. Neuron"),
  select_kegg_term(homer_kegg_list$hESC_vs_Neuron, "hESC vs. Neuron")
)

kegg_plotting_data <- kegg_plotting_data |> mutate(category = factor(category, levels = names(comp_cols))) |>
  group_by(category) |> 
  mutate(Term_unique = reorder(Term, `-log10pval`)) |> 
  ungroup()

final_kegg_plot <- kegg_barplot(kegg_plotting_data, fill_colors = comp_cols, title = "KEGG Pathways")

save(final_kegg_plot, file = file.path(rmats_plots_dir, "final_kegg_plot.rda"))
#------------------------------
load(file.path(rmats_plots_dir, "AS_eventsBarplot.rda"))
load(file.path(rmats_plots_dir, "AS_genesBarplot.rda"))
load(file.path(rmats_plots_dir, "eventsVenndiagram.rda"))
load(file.path(rmats_plots_dir, "geneVenndiagram.rda"))
load(file.path(rmats_plots_dir, "final_kegg_plot.rda"))

# Group them as one figure
pdf(file=file.path(rmats_plots_dir, "pairwise_comparison.pdf"),width=10,height=5.6,bg="transparent")

pageCreate(width = 10, height = 5.6, showGuides = FALSE)
# Event Level
plotText("a", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotText("Differential splicing total events", x = 1.2, y = 0.45, just = c("center", "top"), fontfamily = "Helvetica",fontsize = 8)
plotGG(plot=eventsVenndiagram, x=0.1, y=0.65, height=2, width=2)

# Gene level
plotText("b", x = 0.1, y = 3, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotText("Differential splicing total unique genes", x = 1.2, y = 3.35, just = c("center", "top"), fontfamily = "Helvetica",fontsize = 8)
plotGG(plot=geneVenndiagram, x=0.1, y=3.55, height=2, width=2)


plotText("c", x = 2.3, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotText("Differential splicing events type in pairwise comparison", 
    x = 4.75, y = 0.45, just = c("center", "top"), fontfamily = "Helvetica",fontsize = 8)
plotGG(plot=AS_eventsBarplot, x=2.55, y=0.65, height=2, width=4.25)

plotText("d", x = 2.3, y = 3, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotText("Differential splicing unique genes type in pairwise comparison", 
    x = 4.75, y = 3.35, just = c("center", "top"), fontfamily = "Helvetica",fontsize = 8)
plotGG(plot=AS_genesBarplot, x=2.55, y=3.55, height=2, width=4.25)

# Kegg pathway
plotText("e", x = 7.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(plot=final_kegg_plot, x=7.3, y=0.35, height=5.3, width=2.7)

dev.off()

# PCA test --------------------------------------------------------------------------------
chrs  <- paste0("chr",c(1:22,"X"))
load(file.path(rmats_data_dir, "hesc_npc_jc.Rdata"))
load(file.path(rmats_data_dir, "npc_neuron_jc.Rdata"))
load(file.path(rmats_data_dir, "hesc_neuron_jc.Rdata"))
load(file.path(rmats_data_dir, "rmats_jc_filtered_results.Rdata"))

extract_full_data <- function(maser_obj, label1, label2) {
    event_types <- c("SE", "RI", "MXE", "A3SS", "A5SS")
    all_types_df <- map_dfr(event_types, function(type) {
        df <- summary(maser_obj, type)
        if (nrow(df) == 0) return(NULL)
        df <- df |> filter(Chr %in% chrs)
        coord_cols <- switch(type,
            "SE"   = c("exon_target","exon_upstream","exon_downstream"),
            "MXE"  = c("exon_1", "exon_2","exon_upstream","exon_downstream"),
            "RI"   = c("exon_ir","exon_upstream","exon_downstream"),
            "A3SS" = c("exon_long", "exon_short","exon_flanking"),
            "A5SS" = c("exon_long", "exon_short","exon_flanking")
        )

        if (!all(coord_cols %in% colnames(df))) {
            df <- df |> unite("coords", starts_with("exon_"), sep = "_", remove = FALSE)
            } else {
            df <- df |> unite("coords", all_of(coord_cols), sep = "_", remove = FALSE)
        }

        df |> mutate(EventType = type, event_key = paste(type, Chr, Strand, coords, sep = "_")) |>
            select(event_key, EventType, Chr, GeneID,geneSymbol, PSI_1, PSI_2) |> 
            rename(!!label1 := PSI_1, !!label2 := PSI_2)
        })
return(all_types_df)
}

#combined data
full_hesc_npc_df <- extract_full_data(hesc_npc_filt, "hESC", "NPC") # Total 211,494
full_npc_neuron_df <- extract_full_data(npc_neuron_filt, "NPC", "Neuron") # Total 189,275
full_hesc_neuron_df <- extract_full_data(hesc_neuron_filt, "hESC", "Neuron")

full_all_df <- full_hesc_npc_df |> 
    inner_join(full_npc_neuron_df, by = c("event_key","EventType","Chr","GeneID","geneSymbol","NPC"))  # 163, 421

# Change to Average PSI 
all_matrix_psi <- full_all_df |> separate(hESC, into = c("hESC_1", "hESC_2", "hESC_3"), sep = ",") |> 
    separate(NPC, into =c("NPC_1", "NPC_2", "NPC_3"), sep = ",") |>
    separate(Neuron, into =c("Neuron_1", "Neuron_2", "Neuron_3"), sep = ",") |>
    mutate(across(starts_with(c("hESC","NPC","Neuron")), as.numeric)) 

all_matrix_psi |> filter(if_any(starts_with(c("hESC", "NPC", "Neuron")), is.na)) |> nrow()   # NAs included n=1680
all_matrix_psi_noNA <- all_matrix_psi |> drop_na(starts_with(c("hESC", "NPC", "Neuron"))) # n = 161,741

save(all_matrix_psi_noNA, file = file.path(rmats_data_dir, "all_matrix_psi_noNA.Rdata"))
load(file.path(rmats_data_dir, "all_matrix_psi_noNA.Rdata")) # all_matrix_psi_noNA
#Getting only PSI
psi_mat <- all_matrix_psi_noNA |> select(hESC_1, hESC_2, hESC_3,
                                        NPC_1, NPC_2, NPC_3,
                                        Neuron_1, Neuron_2, Neuron_3) |> as.matrix()
rownames(psi_mat) <- all_matrix_psi_noNA$event_key


psi_mat <- psi_mat[apply(psi_mat, 1, var) > 0, ] # if the variance is 0 across samples filter total n= 1792
pca_res <- prcomp(t(psi_mat), scale. = TRUE)

pca_df <- as.data.frame(pca_res$x)
pca_df$Groups <- factor(rep(c("hESC", "NPC", "Neuron"), each = 3), 
                       levels = c("hESC", "NPC", "Neuron"))
variance_explained <- summary(pca_res)$importance[2, ]

splicePCAplot <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Groups)) +
  geom_point(size=2) +
  scale_color_manual(values = c("hESC" = "#419164", "NPC" = "#F9C555", "Neuron" = "#876EC4")) +
  labs(title = paste0("Total spliceing events n=", nrow(psi_mat)),
    x = paste0("PC1: ", round(variance_explained[1]*100, 2), "%"),
    y = paste0("PC2: ", round(variance_explained[2]*100, 2), "%")) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),    
        axis.line = element_line(colour = "black", linewidth = 0.25),
        legend.position = "none",
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0, 0, 0, 0),
        plot.title = element_text(size=7),
        axis.text.x  = element_text(size = 6),
        axis.text.y  = element_text(size = 6),
        axis.title.y = element_text(size = 7),
        axis.title.x = element_text(size = 7)
)

save(splicePCAplot, file = file.path(rmats_plots_dir, "splicePCAplot.Rdata"))

#-------------------------------------------------------------------------------------------------
# Replicate correlation 
#-------------------------------------------------------------------------------------------------
cond_cols <- c("hESC" = "#8DBDA2","NPC" = "#FBDC99","Neuron" = "#B7A8DC")

psi_sources <- list(
  hESC = list(
    SE   = hesc_npc_jc@SE_PSI,
    RI   = hesc_npc_jc@RI_PSI,
    A3SS = hesc_npc_jc@A3SS_PSI,
    A5SS = hesc_npc_jc@A5SS_PSI,
    MXE  = hesc_npc_jc@MXE_PSI
  ),
  NPC = list(
    SE   = hesc_npc_jc@SE_PSI,
    RI   = hesc_npc_jc@RI_PSI,
    A3SS = hesc_npc_jc@A3SS_PSI,
    A5SS = hesc_npc_jc@A5SS_PSI,
    MXE  = hesc_npc_jc@MXE_PSI
  ),
  Neuron = list(
    SE   = hesc_neuron_jc@SE_PSI,
    RI   = hesc_neuron_jc@RI_PSI,
    A3SS = hesc_neuron_jc@A3SS_PSI,
    A5SS = hesc_neuron_jc@A5SS_PSI,
    MXE  = hesc_neuron_jc@MXE_PSI
  )
)

# Function to calcualte correlation
calc_rep_cor_simple <- function(psi_mat, celltype) {
    cols <- grep(paste0("^", celltype), colnames(psi_mat), value = TRUE)
    mat <- psi_mat[, cols]
    combs <- combn(seq_along(cols), 2)

    apply(combs, 2, function(idx) {
        v1 <- mat[, idx[1]]
        v2 <- mat[, idx[2]]
        keep <- complete.cases(v1, v2)
        cor(v1[keep], v2[keep], method = "pearson")
    })
}


event_types <- c("SE", "RI", "A3SS", "A5SS", "MXE")
celltypes <- c("hESC", "NPC", "Neuron")

cor_df <- map_dfr(event_types, function(ev) {
  map_dfr(celltypes, function(ct) {
    psi <- psi_sources[[ct]][[ev]]
    tibble(event = ev,celltype = ct,cor = calc_rep_cor_simple(psi, ct))
  })
})

cor_df <- cor_df |> mutate( 
  celltype = factor(celltype, levels = c("hESC", "NPC", "Neuron")),
  event = factor(event, levels = c("SE", "RI", "A3SS", "A5SS", "MXE"))
)

pearsonBarplot <- ggplot(cor_df, aes(x = celltype, y = cor, fill = celltype)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.4) +
    geom_jitter(width = 0.12, size = 1, alpha = 0.6) +
    scale_fill_manual(values = cond_cols) +
    facet_wrap(~ event, nrow = 1) +
    coord_cartesian(ylim = c(0.8, 1.0)) +
    #geom_hline(yintercept = 0.8,  linetype="dashed", color = "#c43c3c", linewidth = 0.5) +
    labs(y = "Replicate PSI correlation (Pearson r)",x = NULL) +
    theme_classic() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.25),
        legend.position = "none",
        plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0, 0, 0, 0),
        axis.text.x  = element_text(size = 6, color = "black"),
        axis.text.y  = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 7),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 7, face = "bold")
)

save(pearsonBarplot, file = file.path(rmats_plots_dir, "pearsonBarplot.Rdata"))

#-------------------------------------------------------------------------------------------------
# Anova test across the samples to find differentially spliced events across the 3 groups
#-------------------------------------------------------------------------------------------------
load(file.path(rmats_data_dir, "hesc_neuron_jc.Rdata"))
load(file.path(rmats_data_dir, "hesc_npc_jc.Rdata"))
load(file.path(rmats_data_dir, "npc_neuron_jc.Rdata"))

groups <- factor(rep(c("hESC", "NPC", "Neuron"), each = 3), levels = c("hESC", "NPC", "Neuron"))

design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
fit <- lmFit(psi_mat, design)

contrast.matrix <- makeContrasts(
  hESC_vs_NPC = NPC - hESC,
  NPC_vs_Neuron = Neuron - NPC,
  hESC_vs_Neuron = Neuron - hESC,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

res_ebays <- topTable(fit2, number = Inf, adjust.method = "BH")
res_ebays$event_key <- rownames(res_ebays)
res_all_limma <- res_ebays |> left_join(all_matrix_psi_noNA |> 
                    select("event_key","EventType","Chr","GeneID","geneSymbol"),by = "event_key")

fdr_thresh <- 0.05
deltaPSI_thresh <- 0.1

save(res_ebays,res_all_limma, file = file.path(rmats_data_dir, "limma_anova_results.Rdata"))

load(file.path(rmats_data_dir, "limma_anova_results.Rdata")) # res_ebays,res_all_limma

# Classify patterns based on contrasts
res_classified <- res_all_limma |> mutate(
        # Determine significance and direction for each contrast
        sig_hESC_NPC = adj.P.Val < fdr_thresh & abs(hESC_vs_NPC) > deltaPSI_thresh,
        sig_NPC_Neuron = adj.P.Val < fdr_thresh & abs(NPC_vs_Neuron) > deltaPSI_thresh,
        sig_hESC_Neuron = adj.P.Val < fdr_thresh & abs(hESC_vs_Neuron) > deltaPSI_thresh,
        
        dir_hESC_NPC = sign(hESC_vs_NPC),
        dir_NPC_Neuron = sign(NPC_vs_Neuron),
        dir_hESC_Neuron = sign(hESC_vs_Neuron),
        
        # Define patterns
        # increasing: hESC < NPC < Neuron
        pattern = case_when( sig_hESC_Neuron & dir_hESC_Neuron > 0 & dir_hESC_NPC > 0 & dir_NPC_Neuron > 0 ~ "PSI_increasing",
            # decreasing: hESC > NPC > Neuron
            sig_hESC_Neuron & dir_hESC_Neuron < 0 & dir_hESC_NPC < 0 & dir_NPC_Neuron < 0 ~ "PSI_decreasing",
            # NPC-specific high: hESC < NPC & NPC > Neuron
            sig_hESC_NPC & sig_NPC_Neuron & dir_hESC_NPC > 0 & dir_NPC_Neuron < 0 ~ "NPC_peak",
            # NPC-specific low: hESC > NPC & NPC < Neuron
            sig_hESC_NPC & sig_NPC_Neuron & dir_hESC_NPC < 0 & dir_NPC_Neuron > 0 ~ "NPC_dip",
            # Early change (hESC to NPC)
            sig_hESC_NPC & !sig_NPC_Neuron ~ "hESC_NPC_specific",
            # Late change (NPC to Neuron)
            !sig_hESC_NPC & sig_NPC_Neuron ~ "NPC_Neuron_specific",
            TRUE ~ "Other"
        )
)

res_classified_grouped <- res_classified |> mutate(pattern_grouped = case_when(pattern %in% c("NPC_dip", "NPC_peak") ~ "NPC_specific",TRUE ~ pattern))
sig_all_limma <- res_classified_grouped |> filter(adj.P.Val < 0.05) |> filter(rowSums(across(c(sig_hESC_NPC, sig_NPC_Neuron, sig_hESC_Neuron))) > 0)

# Get PSI data for all 9 samples
heatmap_data <- all_matrix_psi_noNA |> filter(event_key %in% sig_all_limma$event_key) |>
    select(event_key, hESC_1, hESC_2, hESC_3,NPC_1, NPC_2, NPC_3,Neuron_1, Neuron_2, Neuron_3)

heatmap_data_fi <- sig_all_limma |> select(event_key, geneSymbol, pattern_grouped) |> 
    left_join(heatmap_data, by = "event_key") 
heatmap_data_fi <- heatmap_data_fi |> mutate(pattern_grouped = factor(pattern_grouped, 
                                                levels = c("PSI_increasing", "PSI_decreasing", 
                                                           "hESC_NPC_specific", "NPC_specific", 
                                                           "NPC_Neuron_specific"))) |> arrange(pattern_grouped)

heatmap_num_mat <- heatmap_data_fi |> 
                    select(hESC_1, hESC_2, hESC_3, NPC_1, NPC_2, NPC_3, Neuron_1, Neuron_2, Neuron_3) |> 
                    as.matrix()
rownames(heatmap_num_mat) <- heatmap_data_fi$event_key

# Scale by row 
heatmap_mat_scaled <- t(scale(t(heatmap_num_mat)))


sample_meta <- data.frame(Stage = factor(rep(c("hESC", "NPC", "Neuron"), each = 3), 
                  levels = c("hESC", "NPC", "Neuron")))
rownames(sample_meta) <- colnames(heatmap_mat_scaled)

colors_meta <- list(Stage = c("hESC" = "#419164", "NPC" = "#F9C555", "Neuron" = "#876EC4"))

colors_pattern <- c(
    "PSI_increasing" = "#c43d4d",
    "PSI_decreasing" = "#0d8fb3",
    "hESC_NPC_specific" = "#e07653",
    "NPC_specific" = "#808080",
    "NPC_Neuron_specific" = "#b9850c"
)

# Heatmap colors
heatmapColors <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)

col_anno <- HeatmapAnnotation(
    Stage = sample_meta$Stage,
    col = colors_meta,
    gap = unit(0.5, 'mm'),
    annotation_name_gp = gpar(fontsize = 7),
    simple_anno_size = unit(3, "mm"),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
        Stage = list(
            title = "CellType",
            title_gp = gpar(fontsize = 7),
            labels_gp = gpar(fontsize = 6),
            title_position = "leftcenter",
            grid_width = unit(2, "mm"),
            grid_height = unit(1, "mm"),
            at = c("hESC", "NPC", "Neuron"),
            labels = c("hESC", "NPC", "Neuron"),
            ncol = 1
        )
    )
)

heatmap_data_sorted <- heatmap_data_fi |>
    group_by(pattern_grouped) |>
    arrange(
        pattern_grouped,
        case_when(
            pattern_grouped == "PSI_increasing" ~ hESC_1 + hESC_2 + hESC_3,
            pattern_grouped == "PSI_decreasing" ~ -(hESC_1 + hESC_2 + hESC_3),
            pattern_grouped == "hESC_NPC_specific" ~ hESC_1 + hESC_2 + hESC_3,
            pattern_grouped == "NPC_specific" ~ NPC_1 + NPC_2 + NPC_3,
            pattern_grouped == "NPC_Neuron_specific" ~ NPC_1 + NPC_2 + NPC_3,
            TRUE ~ 0
        )) |> ungroup()

heatmap_num_mat <- heatmap_data_sorted |> 
    select(hESC_1, hESC_2, hESC_3, NPC_1, NPC_2, NPC_3, Neuron_1, Neuron_2, Neuron_3) |> 
    as.matrix()

rownames(heatmap_num_mat) <- heatmap_data_sorted$event_key

# Scale by row 
heatmap_mat_scaled <- t(scale(t(heatmap_num_mat)))

row_anno <- rowAnnotation(
    Pattern = heatmap_data_sorted$pattern_grouped,
    col = list(Pattern = colors_pattern),
    simple_anno_size = unit(0.4, "cm"),
    show_annotation_name = FALSE,
    show_legend = FALSE
)

pattern_levels <- levels(heatmap_data_sorted$pattern_grouped)
pattern_colors_for_labels <- colors_pattern[pattern_levels]

hmap <- Heatmap(
    heatmap_mat_scaled,
    col = colorRamp2(seq(-2, 2, length.out = 100), heatmapColors),
    
    # Row settings
    show_row_names = FALSE,
    cluster_rows = FALSE,
    show_row_dend = FALSE,
    row_split = heatmap_data_sorted$pattern_grouped,
    row_gap = unit(2, "mm"),
    row_title_rot = 0,
    row_title_side = "right",
    row_title_gp = gpar(
        fontsize = 8,
        col = pattern_colors_for_labels 
    ),
    border = TRUE,
    
    # Column settings
    show_column_names = FALSE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    row_title = NULL,
    
    # Annotations
    top_annotation = col_anno,
    right_annotation = row_anno
)

heatmapLegend <- Legend(
    at = c(-2, 2),
    col_fun = colorRamp2(breaks = seq(-2, 2, length.out = 100), colors = heatmapColors),
    border = NA,
    title_gp = gpar(fontsize = 0),
    labels_gp = gpar(fontfamily = "Helvetica", fontsize = 8),
    legend_width = unit(2.67, "in"),
    grid_height = unit(0.11, "in"),
    direction = "horizontal"
)

heatmapGrob <- grid.grabExpr(draw(hmap,
                                  show_annotation_legend = FALSE,
                                  show_heatmap_legend = FALSE,
                                  background = "transparent"))

heatmapLegendGrob <- grid.grabExpr(draw(heatmapLegend))

# Save grobs
save(heatmapGrob, file = file.path(rmats_plots_dir, "heatmapGrob.rda"))
save(heatmapLegendGrob, file = file.path(rmats_plots_dir, "heatmapLegendGrob.rda"))

#draw(hmap, merge_legend = TRUE)


# Summary table
pattern_summary <- sig_all_limma |> group_by(pattern_grouped) |>
    summarise(
        n_events = n(),
        n_genes = n_distinct(geneSymbol),
        .groups = "drop") |>
    arrange(desc(n_events))


save(res_classified_grouped, heatmap_data, pattern_summary,sig_all_limma,
     file = file.path(rmats_data_dir, "splicing_results_grouped_hmap.Rdata"))
load(file.path(rmats_data_dir, "splicing_results_grouped_hmap.Rdata")) # res_classified_grouped, heatmap_data, pattern_summary,sig_all_limma

# -----------------------------------------------------------------
# Pathway analysis -------------------------------------------------
pattern_groups <- c("PSI_increasing", "PSI_decreasing", "hESC_NPC_specific", "NPC_specific", "NPC_Neuron_specific")

for (pattern_name in pattern_groups) {  
  pattern_genes <- sig_all_limma |>
    filter(pattern_grouped == pattern_name) |>
    pull(GeneID) |>
    unique()
    pattern_genes_clean <- gsub("\\..*$", "", pattern_genes)
    file_prefix <- paste0("limma_cluster_",tolower(gsub(" ", "_", pattern_name))) 
    output_gene_file <- file.path(rmats_data_dir, paste0(file_prefix, "_genes.txt"))
    
    write.table(pattern_genes_clean, file = output_gene_file,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
    homer_cmd <- paste(
        homer_script,
        output_gene_file,
        file.path(rmats_data_dir, "all_rmats_background_genes.txt"),
        file.path(rmats_homer, file_prefix))

    system(homer_cmd)
}

# Kegg pathway

homer_comparisons <- list(
  "PSI_increasing" = file.path(rmats_homer, "limma_cluster_psi_increasing"),
  "PSI_decreasing" = file.path(rmats_homer, "limma_cluster_psi_decreasing"),
  "hESC_NPC_specific" = file.path(rmats_homer, "limma_cluster_hesc_npc_specific"),
  "NPC_specific" = file.path(rmats_homer, "limma_cluster_npc_specific"),
  "NPC_Neuron_specific" = file.path(rmats_homer, "limma_cluster_npc_neuron_specific")
)

pathway_types <- c("kegg.txt")
homer_kegg_list <- list()

# Finding significant Pathway
for (comp_name in names(homer_comparisons)) {
  homer_dir <- homer_comparisons[[comp_name]]
  input_path <- file.path(homer_dir, "kegg.txt")
  
  if (file.exists(input_path)) {
    pathway_table <- read_delim(input_path, delim = "\t", show_col_types = FALSE) |>
      mutate(pval = exp(1)^logP) |>
      filter(pval < 0.05) |> 
      distinct(Term, .keep_all = TRUE) |>
      mutate(`-log10pval` = -log10(pval)) |>
      dplyr::select(Term, Enrichment, `-log10pval`) |>
      arrange(desc(`-log10pval`))
    
    homer_kegg_list[[comp_name]] <- pathway_table
    
    output_file <- file.path(rmats_homer, "tables", paste0("limma_cluster_kegg_", comp_name, ".csv"))
    if (!dir.exists(dirname(output_file))) { dir.create(dirname(output_file), recursive = TRUE)}
    write_csv(pathway_table, output_file)
  }
}

# Then finding the redundant pathways 
combined_kegg_raw <- data.frame()

for (comp_name in names(homer_comparisons)) {
  input_path <- file.path(homer_comparisons[[comp_name]], "kegg.txt")
  
  if (file.exists(input_path)) {
    pathway_data <- read_delim(input_path, delim = "\t", show_col_types = FALSE)|>
      mutate(pval = exp(1)^logP) |>
      mutate(`-log10pval` = -log10(pval)) |>
      distinct(Term, .keep_all = TRUE) |>
      select(Term, `-log10pval`) |>   
      mutate(cluster = comp_name)
    
    combined_kegg_raw <- bind_rows(combined_kegg_raw, pathway_data)
  }
}

# Finding top 6 pathways per cluster
top_pathways_per_cluster <- combined_kegg_raw |>
  group_by(cluster) |>
  arrange(desc(`-log10pval`)) |>
  slice_head(n = 6) |>
  ungroup()

# Unique pathway list
selected_pathways <- top_pathways_per_cluster |> 
  pull(Term) |> 
  unique()

kegg_matrix_selected <- combined_kegg_raw |>
  filter(Term %in% selected_pathways) |>
  pivot_wider(names_from = cluster, 
              values_from = `-log10pval`,
              values_fill = 0) |>
  column_to_rownames("Term")

# order of cluster
cluster_order <- c("PSI_increasing", "PSI_decreasing", "hESC_NPC_specific", "NPC_specific", "NPC_Neuron_specific")
kegg_matrix_selected <- kegg_matrix_selected[, cluster_order]

pathway_max_cluster <- kegg_matrix_selected |>
  as.data.frame() |>
  rownames_to_column("Term") |>
  pivot_longer(cols = -Term, names_to = "cluster", values_to = "value") |>
  group_by(Term) |>
  slice_max(value, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(Term, max_cluster = cluster)

# ordering pathway
pathway_order <- pathway_max_cluster |>
  mutate(max_cluster = factor(max_cluster, levels = cluster_order)) |>
  arrange(max_cluster) |>
  pull(Term)

# reordering
kegg_matrix_ordered <- kegg_matrix_selected[pathway_order, ]

# annotate row
max_cluster_ordered <- pathway_max_cluster |>
  mutate(max_cluster = factor(max_cluster, levels = cluster_order)) |>
  arrange(max_cluster) |>
  pull(max_cluster)

cluster_colors <- c(
  "PSI_increasing" = "#c43d4d",
  "PSI_decreasing" = "#0d8fb3",
  "hESC_NPC_specific" = "#e07653",
  "NPC_specific" = "#808080",
  "NPC_Neuron_specific" = "#b9850c"
)

row_anno <- rowAnnotation(
  Top_in = max_cluster_ordered,
  col = list(Top_in = cluster_colors),
  width = unit(3, "mm"),
  show_annotation_name = FALSE
)

# Pathway heatmap with grob saving
pathway_hmap <- Heatmap(
  as.matrix(kegg_matrix_ordered),
  # name = "-log10 p",
  
  # Color
  col = colorRamp2(c(0, 2,4,6), 
                   c("#f4f1de", "#fcae61", "#f18179", "#f32c6b")),
  
  # Row 
   row_names_gp = gpar(fontsize = 7),
   row_names_side = "right",
   cluster_rows = FALSE,
  show_row_dend = FALSE,
  row_title = NULL,
  # show_row_names = FALSE,
  # Column
  # column_names_gp = gpar(fontsize = 8),
  # column_names_rot = 45,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  # Border
  border = TRUE,
  
  # Annotation
  right_annotation = NULL
)

pathway_legend <- Legend(
  at = c(0, 2, 4, 6),
  col_fun = colorRamp2(c(0, 2, 4, 6), 
                       c("#f4f1de", "#fcae61", "#f18179", "#f32c6b")),
  border = NA,
  title_gp = gpar(fontsize = 0),
  labels_gp = gpar(fontfamily = "Helvetica", fontsize = 8),
  legend_width = unit(2.9, "in"),
  grid_height = unit(0.11, "in"),
  direction = "horizontal"
)

pathway_hmap_grob <- grid.grabExpr(draw(pathway_hmap,
                                         show_annotation_legend = FALSE,
                                         show_heatmap_legend = FALSE,
                                         background = "transparent"))

pathway_legend_grob <- grid.grabExpr(draw(pathway_legend))

# Save grobs
save(pathway_hmap_grob, file = file.path(rmats_plots_dir, "pathway_hmap_grob.rda"))
save(pathway_legend_grob, file = file.path(rmats_plots_dir, "pathway_legend_grob.rda"))

#draw(pathway_hmap, merge_legend = TRUE)

# -------------------------------------------------------------------------------------------
# Group them as one figure

# Save grobs
load(file = file.path(rmats_plots_dir, "heatmapGrob.rda")) # heatmapGrob
load(file = file.path(rmats_plots_dir, "heatmapLegendGrob.rda")) # heatmapLegendGrob
pdf(file=file.path(rmats_plots_dir, "limma_ebays.pdf"),width=11,height=6.5,bg="transparent")

pageCreate(width = 11, height = 6.5, showGuides = FALSE)

# Event Level
plotText("a", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plotGG(plot=heatmapGrob, x=0.5, y=0.75, height=5.2, width=3)
plotGG(plot=heatmapLegendGrob, x=0.57, y=6, height=0.11, width=2.67)

x_start <- 0.59
y_start <- 0.63

# Create color groups for cell types
color_groups <- list(hESC = "#419164",NPC = "#F9C555",Neuron = "#876EC4")

cell_types <- names(color_groups)
rep_text <- rep(c("rep1", "rep2", "rep3"), times = 3)

for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    color <- color_groups[[cell_type]]
    x_pos <- x_start+0.425  + (i - 1) * 0.85
    
    plotText(cell_type, x = x_pos, y = y_start-0.15,
                     just = c("center", "top"), fontfamily = "Helvetica", fontsize = 8,
                     fontcolor = color, fontface = "bold")
}

for(i in seq_along(rep_text)){
    plotText(rep_text[i], x = x_start + (i-1)*0.3, y = y_start,
             just = c("left", "top"), fontfamily = "Helvetica",fontsize = 8)
}

colors_pattern <- c(
    "PSI_increasing" = "#c43d4d",
    "PSI_decreasing" = "#0d8fb3",
    "hESC_NPC_specific" = "#e07653",
    "NPC_specific" = "#808080",
    "NPC_Neuron_specific" = "#b9850c"
)


colors_pattern2 <- c(
    "Group1" = "#c43d4d",
    "Group2" = "#0d8fb3",
    "Group3" = "#e07653",
    "Group4" = "#808080",
    "Group5" = "#b9850c"
)

n_events <- c(1549, 1519, 877, 279, 178)
total_events <- sum(n_events)

gene_interval <- 0.1
current_y <- 1.05

for(i in 1:1){
    id_p <- names(colors_pattern)[i]
    color <- colors_pattern[[id_p]]
    pattern2 <- names(colors_pattern2)[i]
    top_genes <- sig_all_limma |> 
        filter(pattern_grouped == id_p) |>
        arrange(desc(F), adj.P.Val) |>
        distinct(geneSymbol, .keep_all = TRUE) |>
        head(16) |>
        pull(geneSymbol)
    
    title_y <- current_y + 0.8
    
    plotText(label = paste0(pattern2), 
             x = 3.33, y = title_y,
             just = "center", rot =90,
             fontsize = 8, fontcolor = color, fontface = "bold")
    
    for(gene in top_genes){
        plotText(label = gene, 
                 x = 3.5, y = current_y, 
                 just = c("left", "center"),
                 fontsize = 7, 
                 fontcolor = color, 
                 fontface = "italic", 
                 fontfamily = "Helvetica")
        
        current_y <- current_y + gene_interval
    }
    
    current_y <- current_y + 0.08 
}

for(i in 2){
    id_p <- names(colors_pattern)[i]
    color <- colors_pattern[[id_p]]
    pattern2 <- names(colors_pattern2)[i]
    top_genes <- sig_all_limma |> 
        filter(pattern_grouped == id_p) |>
        arrange(desc(F), adj.P.Val) |>
        distinct(geneSymbol, .keep_all = TRUE) |>
        head(16) |>
        pull(geneSymbol)
    
    title_y <- current_y + 0.8
    
    plotText(label = paste0(pattern2), 
             x = 3.33, y = title_y,
             just = "center", rot =90,
             fontsize = 8, fontcolor = color, fontface = "bold")
    
    for(gene in top_genes){
        plotText(label = gene, 
                 x = 3.5, y = current_y, 
                 just = c("left", "center"),
                 fontsize = 7, 
                 fontcolor = color, 
                 fontface = "italic", 
                 fontfamily = "Helvetica")
        
        current_y <- current_y + gene_interval
    }
    
}

current_y  <- 4.4
for(i in 3){
    id_p <- names(colors_pattern)[i]
    color <- colors_pattern[[id_p]]
    pattern2 <- names(colors_pattern2)[i]
    top_genes <- sig_all_limma |> 
        filter(pattern_grouped == id_p) |>
        arrange(desc(F), adj.P.Val) |>
        distinct(geneSymbol, .keep_all = TRUE) |>
        head(9) |>
        pull(geneSymbol)
    
    title_y <- current_y + 0.4
    
    plotText(label = paste0(pattern2), 
             x = 3.33, y = title_y,
             just = "center", rot =90,
             fontsize = 8, fontcolor = color, fontface = "bold")
    
    for(gene in top_genes){
        plotText(label = gene, 
                 x = 3.5, y = current_y, 
                 just = c("left", "center"),
                 fontsize = 7, 
                 fontcolor = color, 
                 fontface = "italic", 
                 fontfamily = "Helvetica")
        
        current_y <- current_y + gene_interval
    }
    
    current_y <- current_y+0.07
}

for(i in 4){
    id_p <- names(colors_pattern)[i]
    color <- colors_pattern[[id_p]]
    pattern2 <- names(colors_pattern2)[i]
    top_genes <- sig_all_limma |> 
        filter(pattern_grouped == id_p) |>
        arrange(desc(F), adj.P.Val) |>
        distinct(geneSymbol, .keep_all = TRUE) |>
        head(3) |>
        pull(geneSymbol)
    
    title_y <- current_y
    
    plotText(label = paste0(pattern2), 
             x = 3.33, y = title_y,
             just = "center", rot =90,
             fontsize = 8, fontcolor = color, fontface = "bold")
    
    for(gene in top_genes){
        plotText(label = gene, 
                 x = 3.5, y = current_y, 
                 just = c("left", "center"),
                 fontsize = 7, 
                 fontcolor = color, 
                 fontface = "italic", 
                 fontfamily = "Helvetica")
        
        current_y <- current_y + gene_interval
    }
    
    current_y <- current_y+0.07
}

for(i in 5){
    id_p <- names(colors_pattern)[i]
    color <- colors_pattern[[id_p]]
    pattern2 <- names(colors_pattern2)[i]
    top_genes <- sig_all_limma |> 
        filter(pattern_grouped == id_p) |>
        arrange(desc(F), adj.P.Val) |>
        distinct(geneSymbol, .keep_all = TRUE) |>
        head(2) |>
        pull(geneSymbol)
    
    title_y <- current_y + 0.18
    
    plotText(label = paste0(pattern2), 
             x = 3.33, y = title_y,
             just = "center", rot =90,
             fontsize = 8, fontcolor = color, fontface = "bold")
    
    for(gene in top_genes){
        plotText(label = gene, 
                 x = 3.5, y = current_y, 
                 just = c("left", "center"),
                 fontsize = 7, 
                 fontcolor = color, 
                 fontface = "italic", 
                 fontfamily = "Helvetica")
        
        current_y <- current_y + gene_interval
    }
}

plotText("b", x = 4.5, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file.path(rmats_plots_dir, "pathway_hmap_grob.rda")) # pathway_hmap_grob
load(file.path(rmats_plots_dir, "pathway_legend_grob.rda")) # pathway_legend_grob
plotGG(plot=pathway_hmap_grob, x=4.85, y=0.91, height=5, width=5.6)
plotGG(plot=pathway_legend_grob, x=4.85+0.1, y=6.0, height=0.11, width=2.9)

for(i in seq_along(colors_pattern2)){
  pattern <- names(colors_pattern2)[i]
  color <- colors_pattern2[[pattern]]

  plotText(pattern, x = 5.1 +(i-1)*0.59, y = 0.875,
           just = c("left", "center"), fontfamily = "Helvetica", fontsize = 7,
           fontcolor = color, fontface = "bold")
}



dev.off()



# Splicing pattern to see 

#violin_data <- sig_all_limma |>
#    select(event_key, EventType, geneSymbol,
#           hESC_vs_NPC, NPC_vs_Neuron, hESC_vs_Neuron,
#           sig_hESC_NPC, sig_NPC_Neuron, sig_hESC_Neuron) |>
#    filter(sig_hESC_NPC | sig_NPC_Neuron | sig_hESC_Neuron) |>
#    pivot_longer(
#        cols = c(hESC_vs_NPC, NPC_vs_Neuron, hESC_vs_Neuron),
#        names_to = "Comparison",
#        values_to = "deltaPSI"
#    ) |>
#    mutate(
#        Comparison = factor(Comparison, 
#                           levels = c("hESC_vs_NPC", "NPC_vs_Neuron", "hESC_vs_Neuron")),
#        EventType = factor(EventType,
#                          levels = c("SE", "RI", "MXE", "A3SS", "A5SS"))
#    ) |>
#    filter(
#        (Comparison == "hESC_vs_NPC" & sig_hESC_NPC) |
#        (Comparison == "NPC_vs_Neuron" & sig_NPC_Neuron) |
#        (Comparison == "hESC_vs_Neuron" & sig_hESC_Neuron))
#
#
#
#
##-----------------------------------------------------------------------------------------------------------------
## Create figure for results (rmats only first)--------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------------
#











# anova_three_way_pval <- apply(psi_mat, 1, function(x){
#     summary(aov(x~groups))[[1]][["Pr(>F)"]][1]
# })

# anova_adj_pval <- p.adjust(anova_three_way_pval, method = "BH")
# anova_res_df <- data.frame(event_key = rownames(psi_mat),
#                            pvalue = anova_three_way_pval,
#                            adj_pval = anova_adj_pval)

# # calculate each group PSI means
# anova_res_mean_df <- anova_res_df |> mutate(hESC_mean = rowMeans(psi_mat[,1:3]),
#                        NPC_mean = rowMeans(psi_mat[,4:6]),
#                        Neuron_mean = rowMeans(psi_mat[,7:9]),
#                        max_deltaPSI = pmax(abs(hESC_mean - NPC_mean), abs(NPC_mean - Neuron_mean), abs(hESC_mean - Neuron_mean))
#                        )

# sig_anova_df <- anova_res_mean_df |> filter(adj_pval < 0.05 & max_deltaPSI >= 0.1) # n = 3308
                       
# sig_anova_cases <- sig_anova_df |> mutate(pattern= case_when(hESC_mean < NPC_mean & NPC_mean < Neuron_mean ~ "increasing",
#                                 hESC_mean > NPC_mean & NPC_mean > Neuron_mean ~ "decreasing",
#                                 NPC_mean > hESC_mean & NPC_mean > Neuron_mean ~ "NPC_high",
#                                 NPC_mean < hESC_mean & NPC_mean < Neuron_mean ~ "NPC_low",
#                                 abs(hESC_mean - NPC_mean) >= 0.1 & abs(NPC_mean - Neuron_mean) < 0.05 ~ "early",
#                                 abs(hESC_mean - NPC_mean) < 0.05 & abs(NPC_mean - Neuron_mean) >= 0.1 ~ "late",
#                                 TRUE ~ "other"))




