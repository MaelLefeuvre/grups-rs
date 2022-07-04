#' @export
#' @param norm_method User input normalization method
#' @param norm_metric User input normalization metric
#' @param plot_data `.pwd` results file dataframe.
#' @return summary statistic function.
get_norm_values <- function(norm_method, norm_metric, norm_value, plot_data) {
  norm_metric_function <- grups.plots::get_norm_metric(norm_metric)


  norm_value <- switch(norm_method,
    "Raw"              = mean(plot_data[which(plot_data$self == TRUE), ]$avg),
    "Value"            = norm_value,
    "Self-comparisons" = norm_metric_function(
      plot_data[which(plot_data$self == TRUE),]$norm_avg
    ),
    "All"              = norm_metric_function(
      plot_data[which(plot_data$self == FALSE),]$norm_avg
    )
  )

  switch(norm_method,
    "Raw"              = list(Self = norm_value, First = (3 / 2) * norm_value, Unrelated = 2 * norm_value),
    "Self-comparisons" = list(Self = norm_value, First = (3 / 2) * norm_value, Unrelated = 2 * norm_value),
    "All"              = list(Self = norm_value / 2, First = (norm_value / 2) *(3 / 2), Unrelated = norm_value),
    "Value"            = list(Self = 0.5, First = 0.75, Unrelated = 1)
  )
}