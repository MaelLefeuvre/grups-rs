#' @export
#' @importFrom utils read.table
#' @param path path leading to a GRUPS .result summary file
#' @return dataframe containing simulation results for all pairs.
load_res_file <- function(path) {
  read.table(path, sep = "\t", header = TRUE)
}