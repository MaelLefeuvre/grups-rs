#' @export
#' @param norm_method User input normalization method
#' @param norm_metric User input normalization metric
#' @param pwd_data `.pwd` results file dataframe.
#' @return summary statistic function.
get_norm_values <- function(
  norm_method,
  norm_request = FALSE,
  norm_metric,
  norm_value,
  pwd_data
) {
  norm_metric_function <- grups.plots::get_norm_metric(norm_metric)

  norm_value <- switch(norm_method,
    "Pairwise" = norm_metric_function(
      pwd_data[which(pwd_data$Self == FALSE), ]$Norm.Avg, na.rm = TRUE
    ),
    "Self"     = norm_metric_function(
      pwd_data[which(pwd_data$Self == TRUE), ]$Norm.Avg, na.rm = TRUE
    ),
    "Value"    = ifelse(norm_request, 1, norm_value)
  )

  switch(norm_method,
    "Pairwise" = list(Self      = norm_value / 2,
                      First     = (norm_value / 2) * (3  / 2),
                      Second    = (norm_value / 2) * (7  / 4),
                      Third     = (norm_value / 2) * (15 / 8),
                      Unrelated = norm_value
                     ),
    "Self"     = list(Self      = norm_value,
                      First     = (3  / 2) * norm_value,
                      Second    = (7  / 4) * norm_value,
                      Third     = (15 / 8) * norm_value,
                      Unrelated = 2 * norm_value
                     ),
    "Value"   = list(Self      = norm_value / 2,
                     First     = (norm_value / 2) * (3  / 2),
                     Second    = (norm_value / 2) * (7  / 4),
                     Third     = (norm_value / 2) * (15 / 8),
                     Unrelated = norm_value
                    )
  )
}