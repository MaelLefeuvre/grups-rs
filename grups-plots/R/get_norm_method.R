
#' @export
#' @param norm_metric User input normalization metric
#' @param norm_value User input optional set normalization value
#' @return summary statistic function.
get_norm_method <- function(
  norm_request,
  norm_method,
  norm_metric,
  norm_value = NULL
) {

  norm_metric_function <- grups.plots::get_norm_metric(norm_metric)
  if (norm_request) {
    switch(norm_method,
      "Pairwise" = function(x) {
        norm_metric_function(x[which(x$Self == FALSE), ]$Raw.Avg.PWD)
      },
      "Self"     = function(x) {
        norm_metric_function(x[which(x$Self == TRUE), ]$Raw.Avg.PWD)
      },
      "Value"                = function(x) {
        norm_value
      }
    )
  } else {
    function(x) 1
  }
}
