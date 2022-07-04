#' @export
#' @import plotly
#' @param block_dataframe dataframe of block jackknife information.
#' @param pair Pairwise comparison identifier, or label.
#' @param subset string or vector of chromosomes to subset.
#' @return plotly scatterplot
plot_sliding_window <- function(block_dataframe, pair, chromosomes = NA) {

    # Filter out chromosome which were not requested by the user.
    block_dataframe <- block_dataframe[which(block_dataframe$chr %in% chromosomes),]

    plot_title <- paste("Average mismatch rate for", pair)

    plotly::plot_ly(
        data  = block_dataframe,
        x     = ~start,
        y     = ~avg_pwd,
        color = ~as.factor(chr),
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
        yaxis  = list(title="Average mismatch rate"),
        legend = list(title       = list(text = '<b> Chromosome </b>'),
                      orientation = 'h',
                      y            = -0.2,
                      yref         = "paper"
                     )
    ) %>%
    plotly::config(editable    = FALSE,
                   displaylogo = FALSE,
                   scrollZoom  = TRUE
                  )

}