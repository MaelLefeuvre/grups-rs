#' @export
#' @import plotly
#' @importFrom dplyr %>%
#' @importFrom reshape2 melt
#' @param matrix matrix of bhattacharya coefficients between simulation distributions.
plot_or_confidence <- function(or_matrix, predictor = NULL) {
    #print(or_matrix[[predictor,"Unrelated"]])

    odds <- or_matrix[, predictor]
    odds <- odds[!(names(odds) %in% predictor)] # remove predictor auto-comparison

    odds <- sort(abs(odds))  # sort according to increasing values of ORs
    odds <- reshape2::melt(odds)

    fig  <- plotly::plot_ly(type = "scatter", mode = "markers")

    sapply(rownames(odds), FUN=function(rel) {
        this_log_odds = odds[[rel,"value"]]
        raw_odds_ratio = exp(odds[-!(row.names(odds) %in% rel),])

        odds_conf <- exp(this_log_odds + 1.96 * sqrt(1))

        fig <<- fig %>% plotly::add_trace(x=this_log_odds, y = rel, name=rel)
    })
    # add invidible marker (ensure x-axis zero line is displayed)
    fig <- fig %>% plotly::layout(shapes = list(
           list(type = "line",
               y0   = 0,
               y1   = 1,
               yref = "paper",
               x0   = 0,
               x1   = 0,
               line = list(color = "black", dash = "dot")
              )
    ))


    fig
}