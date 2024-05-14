#---- PLOT SCRIPTS -------------------------------------------------------------
#>                                                                            <#
#>    Erik Sandertun Roeed                                                    <#
#>                                                                            <#
#>    Plotting functions                                                      <#
#>                                                                            <#
#>    Contents                                                                <#
#>    function: custom_dapc_plot()                                            <#
#>    function: custom_snapclust_plot()                                       <#
#>    function: custom_snapclust_assignment_dot_plot()                        <#
#>                                                                            <#
#------------------------------------------------------------------------------#

custom_dapc_plot <- function(
    dapc_output,
    fill_colours = RColorBrewer::brewer.pal(n = 3, "Set1"),
    outline_colours = rep("gray10", 3)
)
{
  ind_coordinates <- dapc_output[[1]]
  population_centroids <- dapc_output[[2]]
  pca_axis_percentages <- dapc_output[[3]]
  
  dapc_plot <- ind_coordinates |> 
    dplyr::arrange(`Group`) |> 
    ggplot(aes(x = PC1, y = PC2, fill = `Group`)) +
    xlab(paste0("PC1 (", pca_axis_percentages[1], ") %")) +
    ylab(paste0("PC2 (", pca_axis_percentages[2], ") %")) +
    geom_hline(yintercept = 0, alpha = 0.15) +
    geom_vline(xintercept = 0, alpha = 0.15) +
    geom_point(
      aes(col = `Group`),
      pch = 21,
      size = 2
    ) +
    geom_label(
      data = population_centroids,
      mapping = aes(label = `Population`),
      show.legend = FALSE
    ) +
    guides(
      col = guide_legend(nrow = 1),
      fill = guide_legend(nrow = 1)
    ) +
    scale_colour_manual(values = outline_colours) +
    scale_fill_manual(values = fill_colours) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.background = element_rect(fill = "#ebebeb", colour = "black"),
      panel.grid = element_blank()
    )
  
  return(dapc_plot)
}



custom_snapclust_plot <- function(
    snapclust_output,
    true_group_labels = c(
      "H_americanus" = "H americanus",
      "Hybrid" = "Hybrid",
      "H_gammarus" = "H gammarus"
    ),
    suppress_ticks = FALSE,
    actual_group_label = "Actual group"
)
{
  snapclust_plot <- snapclust_output |> 
    ggplot(
      aes(
        y = Ind,
        x = Probability,
        fill = AssignedGroup
      )
    ) +
    geom_col(show.legend = FALSE, width = 1, ) +
    scale_y_discrete(
      name = actual_group_label,
      position = "right",
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      name = "Prob.",
      expand = c(0, 0),
      n.breaks = 2
    ) +
    scale_fill_viridis_d() +
    facet_grid(
      rows = vars(TrueGroup),
      drop = TRUE,
      scales = "free_y",
      labeller = labeller(
        TrueGroup = true_group_labels
      )
    ) +
    theme(
      axis.text.y = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black"),
      panel.spacing.y = unit(0, "lines"),
      strip.background = element_blank()
    )
  
  if (suppress_ticks)
  {
    snapclust_plot <- snapclust_plot +
      theme(axis.ticks.y = element_blank())
  }
  
  return(snapclust_plot)
}



custom_snapclust_assignment_dot_plot <- function(
    snapclust_output,
    assigned_group_relabels = list(
      `B` = "1",
      `0.125_A-0.875_B` = "0.875",
      `0.25_A-0.75_B` = "0.75",
      `0.5_A-0.5_B` = "0.5",
      `0.75_A-0.25_B` = "0.25",
      `0.875_A-0.125_B` = "0.125",
      `A` = "0"
    ),
    actual_group_title = "Actual group",
    assigned_group_title = "Assigned group"
)
{
  per_ind_data <- snapclust_output |> 
    mutate(
      AssignedGroup = vapply(
        X = AssignedGroup,
        FUN = function(assigned_group_label) {
          do.call(
            switch,
            args = append(
              as.character(assigned_group_label),
              assigned_group_relabels
            )
          )
        },
        FUN.VALUE = character(1)
      )
    ) |> 
    mutate(
      AssignedGroup = factor(
        AssignedGroup, levels = unlist(assigned_group_relabels)
      )
    ) |> 
    group_by(Ind) |> 
    filter(Probability == max(Probability)) |> 
    group_by(TrueGroup)
  
  per_assigned_group_data <- per_ind_data |> 
    count(AssignedGroup) |> 
    mutate(Fraction = n / sum(n))
  
  assignment_plot <- per_assigned_group_data |> 
    ggplot(
      aes(
        x = TrueGroup,
        y = AssignedGroup,
        size = Fraction
      )
    ) +
    xlab(actual_group_title) +
    ylab(assigned_group_title) +
    scale_x_discrete(position = "top") +
    geom_point(show.legend = FALSE) +
    geom_text(
      aes(label = paste(round(Fraction, 3) * 100, "%"), size = 0.25),
      show.legend = FALSE,
      vjust = 2.25,
      hjust = 0.35
    ) +
    theme(
      panel.border = element_rect(fill = NA, colour = "black")
    )
  
  assignment_prob_plot <- per_ind_data |> 
    mutate(
      AssignedGroup = factor(
        # Reverse levels of assigned group to plot more intuitively
        AssignedGroup, levels = rev(unlist(assigned_group_relabels))
      )
    ) |> 
    ggplot(
      aes(
        x = TrueGroup,
        y = Probability,
        fill = TrueGroup
      )
    ) +
    ggtitle(assigned_group_title) +
    xlab(actual_group_title) +
    scale_y_continuous(n.breaks = 3) +
    geom_hline(yintercept = 0.95, lty = 2, alpha = 0.5) +
    geom_violin(show.legend = FALSE, linewidth = 0.25) +
    # Should be scale = "free_x"?
    facet_grid(cols = vars(AssignedGroup), scale = "free", space = "free") +
    theme(
      panel.border = element_rect(fill = NA, colour = "black"),
      panel.spacing.x = unit(0, "lines"),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = -90, vjust = 0.2, hjust = 0),
      plot.title = element_text(size = 11, hjust = 0.5, vjust = -2.5)
    )
  
  plots_combined <- (assignment_plot + assignment_prob_plot) +
    plot_layout(heights = c(4, 1))
  
  return(plots_combined)
}
