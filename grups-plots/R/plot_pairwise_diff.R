#' @export
#' @importFrom dplyr %>%
#' @import plotly
#' @param plot_data data frame
#' @return plotly barplot
plot_pairwise_diff <- function(data, hide_self_comparisons = FALSE, norm_method = NULL, norm_metric = NULL, norm_values = NULL) {
  plot_data <- data$data
  norm_values <- data$norm_values

  # Add horizontal Ms lines as annotations if avg_self is not undefined
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
    plot_data <- plot_data[which(plot_data$self == FALSE),]
  }

  # Plot
  fig <- plotly::plot_ly(type    = "bar",
                  data    = plot_data,
                  #x       = ~pairs,
                  #y       = ~norm_avg,
                  color   = ~rel
                  #error_y = ~list(array = norm_ci,
                  #                color = "#000000"
                 #)
  )

  plotted_rels <- c()
  for (i in seq_along(plot_data$pairs)) {

    already_plotted <- plot_data$rel[i] %in% plotted_rels

    plotted_rels  <- unique(c(plotted_rels, as.character(plot_data$rel[i])))
    fig <- fig %>% plotly::add_trace(
      x           = plot_data$pairs[i],
      y           = plot_data$norm_avg[i],
      color       = plot_data$rel[i],
      hovertext   = paste("Overlap:", plot_data$overlap[i]),
      legendgroup = as.character(plot_data$rel[i]),
      showlegend  = !already_plotted,
      error_y     = list(array = plot_data$norm_ci[i], color = "#000000")
    )
  }

  fig <- fig %>% plotly::layout(title = list(text = "Raw average genetic distances.",
                              y    = 0.99,
                              yref = "paper"
                             ),
                 yaxis = list(title = "Mean genetic distance"
                              #range = c(0, 0.300)
                             ),
                 xaxis       = list(title = "Pairs"),
                 legend = list(title = list(text = "<b>Relationship</b>"),
                               orientation = "h",
                              y            = -0.2,
                              yref         = "paper"
                              ),
                 shapes      = plot_shapes,
                 annotations = plot_annotations,
                 margin      = list(r = 60)
  ) %>%
  plotly::config(editable = TRUE, displaylogo = FALSE, scrollZoom = TRUE)

  fig
}
