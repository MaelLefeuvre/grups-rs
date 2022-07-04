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
      grups.plots::ms_annotation("2Ms", norm_values$Unrelated)
    )
    plot_shapes <- list(
      grups.plots::hline(norm_values$Self),
      grups.plots::hline(norm_values$First),
      grups.plots::hline(norm_values$Unrelated)
    )
  } else {
    plot_annotations <- plot_shapes <- list()
  }

  # Hide self-comparisons from the plot if the user requested it.
  if (hide_self_comparisons) {
    plot_data <- plot_data[which(plot_data$self == FALSE),]
  }

  print(plot_data)
  # Plot
  plotly::plot_ly(type    = "bar",
                  data    = plot_data,
                  x       = ~pairs,
                  y       = ~norm_avg,
                  color   = ~rel,
                  error_y = ~list(array = norm_ci,
                                  color = "#000000"
                 )
  ) %>%
  plotly::layout(title = list(text = "Raw average genetic distances.",
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
}
