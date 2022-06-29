#' Simple horizontal line plotly
#' @export
#' @param y vertical coordinate of the horizontal line
#' @param color color of the horizontal line.
#' @return plotly list type line
hline <- function(y, color = "red") {
  list(
    type = "line",
    y0   = y,
    y1   = y,
    x0   = 0,
    x1   = 1,
    xref = "paper",
    line = list(color = color, dash = "dot")
  )
}
