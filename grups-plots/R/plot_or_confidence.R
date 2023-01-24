#' @export
#' @import plotly
#' @importFrom dplyr %>%
#' @importFrom reshape2 melt
#' @param matrix matrix of bhattacharya coefficients between distributions.
plot_or_confidence <- function(or_matrix, predictor = NULL) {

  odds <- or_matrix[, predictor]
  odds <- odds[!(names(odds) %in% predictor)] # remove auto-comparison
  odds <- odds[order(abs(odds))]  # sort according to increasing values of ORs
  odds <- reshape2::melt(odds)

  fig  <- plotly::plot_ly(
    type   = "scatter",
    mode   = "markers",
    marker = list(
      size = 14,
      line = list(color = "black", width = 2)
    )
  )

  sapply(rownames(odds), FUN = function(rel) {
    this_log_odds  <- odds[[rel, "value"]]
    raw_odds_ratio <- exp(odds[-!(row.names(odds) %in% rel), ])
    odds_conf      <- exp(this_log_odds + 1.96 * sqrt(1))
    fig <<- fig %>% plotly::add_trace(
      x    = this_log_odds,
      y    = rel,
      name = rel
    )
  })
  # add invidible marker (ensure x-axis zero line is displayed)
  fig_title <- sprintf(
    "Odds ratio of obs. within most likely relationship (%s) vs. others",
    predictor
  )
  fig <- fig %>% plotly::layout(
    title  = list(text = fig_title),
    xaxis  = list(title = list(text = "<b>Log Odds Ratio</b>")),
    legend = list(
      orientation = "h",
      title       = list(text = "<b>Relationship</b>"),
      y           = -0.4,
      x           = 0,
      xref        = "paper",
      yref        = "paper",
      xanchor     = "left",
      yanchor     = "bottom"
    ),
    shapes = list(
      list(type = "line",
        y0   = 0,
        y1   = 1,
        yref = "paper",
        x0   = 0,
        x1   = 0,
        line = list(color = "black", dash = "dot")
      )
    )
  ) %>%
  plotly::config(
    editable    = TRUE,
    displaylogo = FALSE,
    scrollZoom  = TRUE,
    toImageButtonOptions = list(
      format   = "svg",
      filename = paste0("Odds-Ratio-Matrix")
    )
  )
}