#' @export
#' @import plotly
#' @import RColorBrewer
#' @importFrom magrittr %>%
#' @param block_dataframe dataframe of block jackknife information.
#' @param pair Pairwise comparison identifier, or label.
#' @param subset string or vector of chromosomes to subset.
#' @return plotly scatterplot
plot_sliding_window <- function(block_dataframe, pair) {

    # Set color palette, and suppress warnings:
    colorpalette <- suppressWarnings(
      RColorBrewer::brewer.pal(
        length(levels(as.factor(block_dataframe$chr))),
        "Set2"
      )
    )

    # Filter out chromosome which were not requested by the user.
    plot_title <- paste("Average mismatch rate for", pair)
    plotly::plot_ly(
        type  = "scatter",
        data  = block_dataframe,
        x     = ~start,
        y     = ~avg_pwd,
        color = ~as.factor(chr),
        colors = colorpalette,
        mode  = "lines+markers"
    ) %>%

    plotly::layout(
        title  = list(text = plot_title,
                      x    = 0.1,
                      y    = 0.99,
                      xref = "paper",
                      yref = "paper"
                     ),
        xaxis  = list(title = "Position (Mb)"),
        yaxis  = list(title = "Average mismatch rate", range = c(0, 1)),
        legend = list(title       = list(text = "<b> Chromosome </b>"),
                      orientation = "h",
                      y            = -0.2,
                      yref         = "paper"
                     )
    ) %>%
    plotly::config(editable    = TRUE,
                   displaylogo = FALSE,
                   scrollZoom  = TRUE,
                   toImageButtonOptions = list(
                    format = "svg",
                    filename = paste0("sliding_window-", pair)
                    )
                  )
}