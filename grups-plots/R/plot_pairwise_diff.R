#' @export
#' @importFrom dplyr %>%
#' @import plotly
#' @param plot_data data frame
#' @return plotly barplot
plot_pairwise_diff <- function(
  data,
  hide_self_comparisons = FALSE,
  norm_method = NULL,
  norm_metric = NULL,
  norm_values = NULL
) {
  plot_data   <- data$data
  norm_values <- data$norm_values
  # Add horizontal Ms lines as annotations if avg_self is defined
  if (!is.nan(norm_values$Self)) {
    plot_annotations <- list(
      grups.plots::ms_annotation("Ms", norm_values$Self),
      grups.plots::ms_annotation("(3/2)Ms", norm_values$First),
      grups.plots::ms_annotation("(7/4)Ms", norm_values$Second),
      grups.plots::ms_annotation("(15/8)Ms", norm_values$Third),
      grups.plots::ms_annotation("2Ms", norm_values$Unrelated)
    )
    plot_shapes <- list(
      grups.plots::hline(norm_values$Self),
      grups.plots::hline(norm_values$First),
      grups.plots::hline(norm_values$Second),
      grups.plots::hline(norm_values$Third),
      grups.plots::hline(norm_values$Unrelated)
    )
  } else {
    plot_annotations <- plot_shapes <- list()
  }

  # Hide self-comparisons from the plot if the user requested it.
  if (hide_self_comparisons) {
    plot_data <- plot_data[which(plot_data$Self == FALSE), ]
  }

  # Plot
  fig <- plotly::plot_ly(type    = "bar",
                  data    = plot_data,
                  color   = ~Rel
  )

  # We use a for loop to add individual traces because plotly's error bar
  # display is broken when using data factoring...
  plotted_rels <- c()
  for (i in seq_along(plot_data$Pair_name)) {

    # Adding traces individually has the effect of duplicating legend traces.
    # -> We need to keep track of what relationship was already plotted to only
    #    Display them once within the legend.
    already_plotted <- plot_data$Rel[i] %in% plotted_rels
    plotted_rels  <- unique(c(plotted_rels, as.character(plot_data$Rel[i])))

    # Add a trace for this pair.
    fig <- fig %>% plotly::add_trace(
      x           = plot_data$Pair_name[i],
      y           = plot_data$Norm.Avg[i],
      color       = plot_data$Rel[i],
      hovertext   = paste("Overlap:", plot_data$Raw.Overlap[i]),
      legendgroup = as.character(plot_data$Rel[i]),
      showlegend  = !already_plotted,
      error_y     = list(array = plot_data$Norm.CI[i], color = "#000000")
    )
  }

  fig %>% plotly::layout(
    title       = list(text = "Raw average genetic distances.",
                       y    = 0.99,
                       yref = "paper"
                      ),
    yaxis       = list(title = "Mean genetic distance"),
    xaxis       = list(title = "Pairs"),
    legend      = list(title       = list(text = "<b>Relationship</b>"),
                       orientation = "h",
                       y           = -0.2,
                       yref        = "paper"
                      ),
    shapes      = plot_shapes,
    annotations = plot_annotations,
    margin      = list(r = 60)
  ) %>%
  plotly::config(
    editable             = TRUE,
    displaylogo          = FALSE,
    scrollZoom           = TRUE,
    toImageButtonOptions = list(
      format   = "svg",
      filename = "raw-genetic-distances"
    )
  )
}
