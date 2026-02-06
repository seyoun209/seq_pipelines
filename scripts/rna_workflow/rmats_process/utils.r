library(rrvgo)
library(plotgardener)
library(grid)
library(RColorBrewer)
library(ggtext)

filterByChromosome <- function(maser_obj, chromosomes = paste0("chr", c(1:22, "X"))) {
    
    as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
    events <- as(maser_obj, "list")
    events_new <- list()
    
    for (type in as_types) {
        # Get GRanges
        gr <- events[[paste0(type, "_gr")]]
        
        if (length(gr) > 0) {
            # Filter by chromosome
            valid_idx <- lapply(gr, function(g) {
                as.character(seqnames(g)) %in% chromosomes
            })
            
            # All exons must be in valid chromosomes
            keep_idx <- which(Reduce("&", valid_idx))
            
            # Get event IDs to keep
            if (length(keep_idx) > 0) {
                res_id <- events[[paste0(type, "_events")]]$ID[keep_idx]
            } else {
                res_id <- integer(0)
            }
            
            # Filter all slots
            slots <- maser:::filterByIds(type, events, res_id)
            
            for (attrib in names(slots)) {
                events_new[[attrib]] <- slots[[attrib]]
            }
        }
    }
    
    # Copy metadata
    events_new[["n_cond1"]] <- events$n_cond1
    events_new[["n_cond2"]] <- events$n_cond2
    events_new[["conditions"]] <- events$conditions
    
    return(as(events_new, "Maser"))
}



#------------------
reduceGO <- function(goterms, category, ont = "BP", threshold = 0.8){

  # Calculate similarity matrix for GO terms based on ontology
  simMatrix <- calculateSimMatrix(goterms$TermID,
                                  orgdb = "org.Hs.eg.db",
                                  ont = ont,
                                  method = "Rel")

  # Create named vector of scores
  # Higher is better -> -log10 transforming p-values
  scores <- setNames(-log10(goterms$pval), goterms$TermID)

  # Group GO terms based on similarity threshold
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold = threshold,
                                  orgdb = "org.Hs.eg.db") |>
    dplyr::rename(`-log10pval` = score)

  # Join grouped terms with original data
  joined_reducedTerms <- left_join(goterms, reducedTerms,
                                   by = join_by("TermID" == "go", "Term" == "term")) |>
    dplyr::select(-cluster) |>
    relocate(parent, .after = Term) |>
    relocate(parentTerm, .after = parent) |>
    dplyr::rename(parentTermID = parent) |>
    mutate("category" = category)

  return(joined_reducedTerms)
}


# Making TABLE for GO results to save it as CSV
write_GO_table <- function(go_df) {
  go_df |> dplyr::select(-`Entrez Gene IDs`, -pval, -logP, -size, -termUniqueness, -termUniquenessWithinCluster, -termDispensability, -category) |>
    relocate(`-log10pval`, .after = Enrichment) |>
    arrange(desc(`-log10pval`))
}

# GO barplot function

gg_GO_barplot <- function(df,fill_colors,title = "GO Terms") {
  ggplot(df, aes(x = `-log10pval`, y = Term_unique, fill = category)) +
    geom_vline(xintercept = 2, color = "grey95", alpha = 0.4) +
    geom_vline(xintercept = 4, color = "grey95", alpha = 0.4) +
    geom_vline(xintercept = 6, color = "grey95", alpha = 0.4) +
    geom_vline(xintercept = 8, color = "grey95", alpha = 0.4) +
    geom_vline(xintercept = 10, color = "grey95", alpha = 0.4) +
    geom_vline(xintercept = 12, color = "grey95", alpha = 0.4) +
    geom_bar(stat = "identity") +
    scale_x_continuous(limits = c(0, 13), expand = c(0, 0),
                       name = "-log~10~pval",
                       breaks = seq(0, 12, 2)) +
    scale_fill_manual(values = fill_colors) +
facet_wrap(~category, ncol = 1, strip.position = "left",
               scales = "free_y") +
    geom_text(aes(x = 0, label = Term), hjust = 0, family = "Helvetica",
              size = 2) +  # Still display original Term name
    theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
          plot.background = element_rect(fill = 'transparent', color = "transparent"),
          text = element_text(family = "Helvetica"),
          legend.position = "None",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_markdown(size = 6),
          axis.text.x = element_text(color = "black", size = 4),
          strip.background = element_blank(),
          axis.ticks = element_blank(),
          axis.line.x = element_line(linewidth = 0.25),
          strip.text = element_text(size = 8, color = "black"),
          panel.spacing = unit(0, "mm"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 8)) +
    ggtitle(title)
}

kegg_barplot <- function(df, fill_colors, title = "Pathways") {
  ggplot(df, aes(x = `-log10pval`, y = Term_unique, fill = category)) +
    geom_vline(xintercept = 1, color = "grey95", alpha = 0.4) +
    geom_vline(xintercept = 2, color = "grey95", alpha = 0.4) +
    geom_vline(xintercept = 3, color = "grey95", alpha = 0.4) +
    geom_vline(xintercept = 4, color = "grey95", alpha = 0.4) +
    #geom_vline(xintercept = 10, color = "grey95", alpha = 0.4) +
    #geom_vline(xintercept = 12, color = "grey95", alpha = 0.4) +
    geom_bar(stat = "identity") +
    scale_x_continuous(limits = c(0, 5), expand = c(0, 0),
                       name = "-log~10~pval",
                       breaks = seq(0, 4, 1)) +
    scale_fill_manual(values = fill_colors) +
facet_wrap(~category, ncol = 1, strip.position = "left",
               scales = "free_y") +
    geom_text(aes(x = 0, label = Term), hjust = 0, family = "Helvetica",
              size = 2) +  # Still display original Term name
    theme(panel.background = element_rect(fill = 'transparent', color = "transparent"),
          plot.background = element_rect(fill = 'transparent', color = "transparent"),
          text = element_text(family = "Helvetica"),
          legend.position = "None",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_markdown(size = 7),
          axis.text.x = element_text(color = "black", size = 7),
          strip.background = element_blank(),
          axis.ticks = element_blank(),
          axis.line.x = element_line(linewidth = 0.25),
          strip.text = element_text(size = 7, color = "black"),
          panel.spacing = unit(0, "mm"),
          plot.title = element_text(hjust = 0.5, size = 8)) +
    ggtitle(title)
}
