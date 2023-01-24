#' @export
#' @import shiny
#' @import plotly
#' @import catmaply
#' @import viridis
#' @importFrom dplyr %>%
#' @importFrom reshape2 melt
#' @param kinship_matrix Kinship matrix
#' @param dimensions Shiny input dimensions (for dynamic resizing)
#' @param order Order in which kinship degrees should be arranged.
#' @return a Plotted kinship matrix
plot_kinship_matrix <- function(kinship_matrix, dimensions, order) {
  plot_title   <- "<b>Kinship matrix</b>"
  legend_title <- "<b>Relationship</b>"

  kinship_matrix$rel <- factor(
    x      = kinship_matrix$rel,
    levels = order
  )

  color_palette <- viridis::viridis_pal(
    option    = "D",
    direction = -1
  )(length(unique(kinship_matrix$rel)))

  tickfont        <- list(size = 16)
  legendfont      <- list(size = 16)
  legendtitlefont <- list(size = 20)


  catmaply::catmaply(
    kinship_matrix,
    x              = Left_ind,
    x_tickangle    = 60,
    y              = Right_ind,
    y_tickangle    = 20,
    z              = rel,
    rangeslider    = FALSE,
    color_palette  = color_palette,
    hover_template = paste(
      "<b>Pair</b>:", Left_ind, "-", Right_ind,
      "<br><b>Relationship:</b>:", rel,
      "<extra></extra>"
    ),
  ) %>% plotly::layout(
    title    = list(text = plot_title),
    autosize = TRUE,
    width    = 1.3 * as.numeric(dimensions[2]),
    height   = 1.3 * as.numeric(dimensions[2]),
    xaxis    = list(
      fixedrange = FALSE,
      autorange  = TRUE,
      automargin = TRUE,
      tickfont   = tickfont
    ),
    yaxis    = list(
      fixedrange = FALSE,
      autorange  = TRUE,
      automargin = TRUE,
      tickfont   = tickfont
    ),
    legend   = list(
        font    = legendfont,
        title   = list(text = legend_title, font = legendtitlefont),
        y       = 0,
        x       = 1,
        yanchor = "bottom",
        xanchor = "right",
        yref    = "paper",
        xref    = "paper"
    )
  ) %>%
  plotly::config(
    editable             = TRUE,
    displaylogo          = FALSE,
    scrollZoom           = TRUE,
    toImageButtonOptions = list(
      format   = "svg",
      filename = "kinship-matrix"
    )
  )
}
