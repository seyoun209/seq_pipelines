library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplotify)

count_summary <- violin_data |>
    group_by(Comparison, EventType) |>
    summarise(n_events = n(), 
              n_genes = n_distinct(geneSymbol),
              .groups = "drop")

# Total counts per comparison
total_counts <- violin_data |>
    group_by(Comparison) |>
    summarise(total_events = n(),
              total_genes = n_distinct(geneSymbol),
              .groups = "drop")



#==============================================================================================
# Create horizontal plots with dots
#==============================================================================================

# Color scheme for comparisons
comparison_colors <- c(
    "hESC_vs_NPC" = "#b8dcc8",
    "NPC_vs_Neuron" = "#fdeab3",
    "hESC_vs_Neuron" = "#d4c3e3"
)

# Title colors
title_colors <- c(
    "hESC_vs_NPC" = "#419164",
    "NPC_vs_Neuron" = "#F9C555",
    "hESC_vs_Neuron" = "#876EC4"
)

# Comparison labels
comparison_labels <- c(
    "hESC_vs_NPC" = "hESC vs NPC",
    "NPC_vs_Neuron" = "NPC vs Neuron",
    "hESC_vs_Neuron" = "hESC vs Neuron"
)

# Reorder Comparison and EventType
violin_data <- violin_data |>
    mutate(
        Comparison = factor(Comparison, 
                           levels = c("hESC_vs_NPC", "NPC_vs_Neuron", "hESC_vs_Neuron")),
        EventType = factor(EventType, 
                          levels = c("MXE","A5SS", "A3SS", "RI", "SE")))

count_summary <- count_summary |>
    mutate(
        Comparison = factor(Comparison, 
                           levels = c("hESC_vs_NPC", "NPC_vs_Neuron", "hESC_vs_Neuron")),
        EventType = factor(EventType, 
                          levels = c("MXE","A5SS", "A3SS", "RI", "SE")),
        label = paste0("Events (n=", n_events,")", "\nGenes (n=", n_genes,")"))

# Create summary labels for each comparison (bold total)
total_counts <- total_counts |>
    mutate(
        Comparison = factor(names(comparison_labels), 
                           levels = c("hESC_vs_NPC", "NPC_vs_Neuron", "hESC_vs_Neuron")),
        label = paste0("Events: ", total_events, "\nGenes: ", total_genes))

# Create horizontal plot with dots
p_violin <- ggplot(violin_data, aes(y = EventType, x = deltaPSI, fill = Comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_jitter(color = "grey30", height = 0.2, alpha = 0.4, size = 0.6) +
    geom_boxplot(width = 0.4, outlier.shape = NA, 
                color = "black", alpha = 0.7) +
    facet_wrap(~Comparison, ncol = 1, 
              labeller = labeller(Comparison = comparison_labels),
              strip.position = "left") +
    scale_fill_manual(values = comparison_colors) +
    scale_x_continuous(breaks = seq(-1, 1, 0.25)) +
    labs(
        y = "Event Type",
        x = expression(Delta*"PSI")
    ) +
    theme_classic() +
    theme(
        strip.background = element_blank(),
        strip.text.y.left = element_text(
            size = 10, 
            face = "bold", 
            angle = 90,
            hjust = 0.5
        ),
        strip.placement = "outside",
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.title.y = element_blank(),
        legend.position = "none"
        #panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        #panel.spacing.y = unit(0.5, "lines")
    )

# Manually color strip text by comparison
g <- ggplot_gtable(ggplot_build(p_violin))
strip_indices <- which(grepl('strip-l', g$layout$name))
colors_vector <- title_colors[levels(violin_data$Comparison)]

for(i in seq_along(strip_indices)) {
    j <- strip_indices[i]
    g$grobs[[j]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- colors_vector[i]
}

# Convert back to ggplot for annotation
p_violin_colored <- ggplotify::as.ggplot(g)

# Add event and gene counts for each event type
p_violin_annotated <- p_violin +
    geom_text(
        data = count_summary,
        aes(y = EventType, x = -1.05, label = label),
        size = 2.5,
        hjust = 1,
        inherit.aes = FALSE
    ) 
        # geom_text(
        #     data = total_counts,
        #     aes(x = -1.05, y = 5.5, label = label),
        #     hjust = 1,
        #     vjust = 1,
        #     size = 3,
        #     fontface = "bold",
        #     inherit.aes = FALSE
        # )

# Apply strip colors manually
g_final <- ggplot_gtable(ggplot_build(p_violin_annotated))
strip_indices <- which(grepl('strip-l', g_final$layout$name))

for(i in seq_along(strip_indices)) {
    j <- strip_indices[i]
    g_final$grobs[[j]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- colors_vector[i]
}
grid::grid.draw(g_final)

ggsave(file.path(rmats_plots_dir, "deltaPSI_horizontal_by_eventtype.pdf"), 
       g_final, width = 8, height = 10)
