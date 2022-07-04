#' @export
#' @param norm_metric User input normalization metric
#' @return summary statistic function.
get_norm_metric <- function(norm_metric) {
  switch(norm_metric,
    "Median" = median,
    "Mean"   = mean,
    "Min"    = min,
    "Max"    = max,
  )
}