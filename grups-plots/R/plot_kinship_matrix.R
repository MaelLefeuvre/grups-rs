#' @export
#' @import shiny
#' @import plotly
#' @import catmaply
#' @import viridis
#' @importFrom dplyr %>%
#' @importFrom reshape2 melt
#' @param kinship_matrix Kinship matrix
#' @param dimensions Shiny input dimensions (helps us retrieve dimensions of the window)
#' @param order Order in which kinship degrees should be arranged.
#' @return a Plotted kinship matrix
plot_kinship_matrix <- function(kinship_matrix, dimensions, order) {
  plot_title <- "Kinship matrix"

  kinship_matrix$rel <- factor(
    kinship_matrix$rel,
    levels = order
  )

  catmaply::catmaply(
    kinship_matrix,
    x           = Left_ind,
    x_tickangle = 60,
    y           = Right_ind,
    y_tickangle = 20,
    z           = rel,
    rangeslider = FALSE,
    categorical_color_range = TRUE,
    categorical_col = rel,
    color_palette = viridis::viridis,
    hover_template = paste(
      "<b>Pair</b>:", Left_ind, "-", Right_ind,
      "<br><b>Relationship:</b>:", rel,
      "<extra></extra>"
  ),
  ) %>% plotly::layout(
    title = list(text = "Kinship Matrix"),
    autosize=TRUE,
    width =  1.3 * as.numeric(dimensions[2]),
    height = 1.3 * as.numeric(dimensions[2]),
    xaxis = list(
      fixedrange=FALSE,
      autorange=TRUE,
      automargin=TRUE
    ),
    yaxis = list(
      #scaleanchor="x",
      #scaleratio=0.5,
      fixedrange=FALSE,
      autorange=TRUE,
      automargin=TRUE
    ),
    legend = list(
        title=list(text='<b> Relationship </b>'),
        orientation = "h",
        #xanchor = "center",
        #yanchor = "center",
        #x = 1,
        y = 1.1,
        yanchor = "bottom",
        yref = "paper"
    )
  ) %>% plotly::config(editable    = TRUE,
      displaylogo = FALSE,
      scrollZoom  = TRUE,
      toImageButtonOptions = list(
        format   = "svg",
        filename = "kinship-matrix"
      )
    )
}