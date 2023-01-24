#' @export
#' @import plotly
#' @importFrom dplyr %>%
#' @importFrom reshape2 melt
#' @param matrix matrix of bhattacharya coefs between sim distributions.
plot_bc_matrix <- function(
  bc_matrix,
  pair,
  plot_title      = NULL,
  marker_text     = NULL,
  cutoff_values,
  cutoff_labels,
  cutoff_colors   =  c("#44AA99", "#DDCC77", "#CC6677"),
  right_align     = FALSE,
  absolute_values = FALSE
) {

  bc_matrix <- reshape2::melt(bc_matrix)

  bc_matrix$signif <- cut(
    if (absolute_values) abs(bc_matrix$value) else bc_matrix$value,
    breaks         = cutoff_values,
    include.lowest = TRUE,
    right          = right_align,
    ordered_result = TRUE,
    labels         = cutoff_labels
  )

  xrange <- c(0.5, length(unique(bc_matrix$Val2)) + 0.5)
  yrange <- c(0.5, length(unique(bc_matrix$Val1)) + 0.5)

  x_axis1 <- list(
    autoaxis    = FALSE,
    showgrid    = FALSE,
    showline    = FALSE,
    zeroline    = FALSE,
    tickangle   = 45,
    title       = "",
    range       = xrange,
    scaleanchor = "y",
    side        = "bottom"
  )

  y_axis1 <- list(
    autoaxis    = FALSE,
    showgrid    = FALSE,
    showline    = FALSE,
    zeroline    = FALSE,
    tickangle   = 0,
    title       = "",
    rangemode   = "tozero",
    range       = yrange,
    scaleanchor = "x"
  )

  marker_size <- (180 / length(unique(bc_matrix$Var2))) - 1

  plotly::plot_ly(
    data   = as.data.frame(bc_matrix),
    x      = ~Var2,
    y      = ~Var1,
    text   = ~paste(marker_text, value),
    type   = "scatter",
    mode   = "markers",
    marker = list(
      size         = marker_size,
      opacity      = 1,
      reversescale = TRUE,
      line         = list(color = "#111111")
    ),
    symbol = I("square"),
    color  = ~signif,
    colors = cutoff_colors
  ) %>%

  plotly::layout(
    xaxis      = x_axis1,
    yaxis      = y_axis1,
    title      = list(text = plot_title),
    hovermode = "closest",
    legend     = list(
      title       = list(text = "<b>Value</b>", size = 11),
      orientation = "h",
      font        = list(size = 11),
      x           = -0.2,
      y           = -0.4,
      xref        = "paper",
      yref        = "paper",
      xanchor     = "left",
      yanchor     = "top"
    )
  ) %>%
  plotly::config(
    editable    = TRUE,
    displaylogo = FALSE,
    scrollZoom  = TRUE,
    toImageButtonOptions = list(
      format = "svg",
      filename = paste0("BC-matrix")
    )
  )
}