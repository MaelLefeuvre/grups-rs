
#' @export
#' @param norm_metric User input normalization metric
#' @param norm_value User input optional set normalization value
#' @return summary statistic function.
get_norm_method <- function(norm_method, norm_metric, norm_value = NULL) {

  norm_metric_function <- grups.plots::get_norm_metric(norm_metric)
  switch(norm_method,
    "Raw"              = function(x) {1},
    "Value"            = function(x) {norm_value},
    "Self-comparisons" = function(x) {norm_metric_function(x[which(x$self == TRUE),]$avg)},
    "All"              = function(x) {norm_metric_function(x[which(x$self == FALSE),]$avg)}
  )
}
