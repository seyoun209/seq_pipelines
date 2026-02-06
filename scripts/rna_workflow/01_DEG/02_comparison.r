library(eulerr)
library(ggVennDiagram)

load(file= file.path(diff_data_dir,"shrink_normoxia_hgnc_tximeta.Rdata")) # res_Shrink_norm_df_hgnc
load(file= file.path(diff_data_dir,"shrink_hypoxia_hgnc_tximeta.Rdata")) # res_Shrink_hypo_df_hgnc
load(file= file.path(diff_data_dir,"shrink_interaction_hgnc_tximeta.Rdata")) # res_Shrink_inter_df_hgnc

diff_norm <- res_Shrink_norm_df_hgnc |> filter(class != "static")  # 1281
diff_hyp <- res_Shrink_hypo_df_hgnc |> filter(class != "static") # 593
diff_inter <- res_Shrink_inter_df_hgnc |> filter(class != "static") # 294

diff_norm_ensg <- diff_norm$gene_id |> unique()
diff_hyp_ensg <- diff_hyp$gene_id |> unique()
diff_inter_ensg <- diff_inter$gene_id |> unique()
onlyInteraction  <- setdiff(diff_inter_ensg, c(diff_norm_ensg, diff_hyp_ensg)) |> as.matrix()
norm_sub <- res_Shrink_norm_df_hgnc |> filter(gene_id  %in% onlyInteraction) 


venn_list <- list(
  "Normoxia" = diff_norm_ensg,
  "Hypoxia" = diff_hyp_ensg,
  "Interaction" = diff_inter_ensg
)

# Fit euler diagram (area-proportional)
venn_fit <- euler(venn_list)

# Plot with custom styling
deg_venn_grob <- plot(
  venn_fit,
  fills = list(
    fill = c("#A7D3D4", "#df91a3", "#f5d7a3"),
    alpha = 0.4
  ),
  edges = list(col = NA),
  quantities = list(cex = 0.8, col = "black", font = 2),
  labels = list(
    col = "black",
    font = 2,
    cex = 0.9
  )
)
deg_venn_grob

# Convert to ggplot for saving/combining
deg_venn_plot <- ggplot() +
  annotation_custom(
    grob = grid::grobTree(deg_venn_grob)
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

save(deg_venn_plot, file = file.path(diff_plot_dir, "deg_venn_plot.rda"))



#--------------Finding genes 

norm_hyp <- intersect(diff_norm_ensg, diff_hyp_ensg)
norm_inter <- intersect(diff_norm_ensg, diff_inter_ensg)
hyp_inter <- intersect(diff_hyp_ensg, diff_inter_ensg)
all_three <- Reduce(intersect, list(diff_norm_ensg, diff_hyp_ensg, diff_inter_ensg))

# Get only interaction genes
onlyInteraction <- setdiff(diff_inter_ensg, union(diff_norm_ensg, diff_hyp_ensg))

# Function to convert ENSG to SYMBOL (or keep ENSG if no SYMBOL)
ensg_to_symbol <- function(ensg_ids, df_with_symbols) {
  symbols <- sapply(ensg_ids, function(id) {
    symbol <- df_with_symbols$SYMBOL[df_with_symbols$gene_id == id][1]
    if (is.na(symbol) || symbol == "") {
      return(id)  # Return ENSG if no SYMBOL
    } else {
      return(symbol)
    }
  })
  return(symbols)
}

all_three_symbols <- ensg_to_symbol(all_three, res_Shrink_inter_df_hgnc)
onlyInteraction_symbols <- ensg_to_symbol(onlyInteraction, res_Shrink_inter_df_hgnc)

# Separate ENSG and non-ENSG for all_three
ensg_all_three <- all_three_symbols[grep("^ENSG", all_three_symbols)]
non_ensg_all_three <- all_three_symbols[!grepl("^ENSG", all_three_symbols)]

# Sort and combine
sorted_all_three <- c(sort(non_ensg_all_three), ensg_all_three)

# Separate ENSG and non-ENSG for onlyInteraction
ensg_only_inter <- onlyInteraction_symbols[grep("^ENSG", onlyInteraction_symbols)]
non_ensg_only_inter <- onlyInteraction_symbols[!grepl("^ENSG", onlyInteraction_symbols)]

# Sort and combine
sorted_only_inter <- c(sort(non_ensg_only_inter), ensg_only_inter)