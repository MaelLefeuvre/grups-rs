#' @export
#' @param pwd_df a dataframe of raw PWD results (.pwd    file extension)
#' @param res_df a dataframe of pedsims results (.result file extension)
#' @return either pwd_df or a merged dataframe of both if res_df is defined
#'         (merge across pair names).
merge_pwd_results <- function(pwd_df = NULL, res_df) {
  # Merge only if res_df is defined and exists. If not, simply return the
  # original raw PWD dataframe.

  # Select which columns of pwd_df should be kept
  subset <- c("Pair_name", "Raw.Overlap")
  if (hasArg(res_df)) {
    merge(
      res_df,
      pwd_df$data[, subset],
      by.x = "Pair_name",
      by.y = "Pair_name"
    )
  } else {
    pwd_df
  }
}