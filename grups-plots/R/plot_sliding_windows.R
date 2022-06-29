#' @export
#' @import plotly
#' @param block_dataframe dataframe of block jackknife information.
#' @param subset string or vector of chromosomes to subset.
#' @return plotly plot 
plot_sliding_window <- function(block_dataframe, subset = NA) {

    plotly::plot_ly(block_dataframe, x = ~start, y = ~pwd/overlap, color = ~as.factor(chr), mode="lines")
}