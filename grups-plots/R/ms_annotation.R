#' Simple annotation for MS values of PWD plot.
#' @export
#' @param label label of the annotation
#' @param y coordinate of the text annotation
#' @return plotly list annotation
ms_annotation <- function(label, y) {
  list(
    x         = 1,
    y         = y,
    xref      = "paper",
    yref      = y,
    text      = label,
    angle     = 45,
    showarrow = FALSE,
    xanchor   = "left",
    yanchor   = "center"
  )
}
