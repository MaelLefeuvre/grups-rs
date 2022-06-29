#' @export
#' @importFrom utils read.table
#' @param path path leading to a GRUPS .blk results file
#' @return block dataframe with columns "chr" "start" "end" "overlap" "pwd"
load_blockfile <- function(path) {
  read.table(
    path,
    sep = "\t",
    header = FALSE,
    col.names = c("chr", "start", "end", "overlap", "pwd")
  )
}