#' @export
#' @importFrom magrittr %>%
#' @importFrom utils read.table
#' @import plotly
#' @import RColorBrewer
#' @param block_dataframe dataframe of raw simulation results
#' @param pair Pairwise comparison identifier, or label.
#' @return plotly violin plot
plot_pedigree_sims <- function(
  sims_dataframe,
  results_data,
  labels_to_plot,
  pair
) {

  # remove labels which were not requested by the user:
  sims_dataframe <- sims_dataframe[
    which(sims_dataframe$label %in% labels_to_plot),
  ]

  # Compute mean_pwd according to label:
  simulation_avgs <- with(sims_dataframe, tapply(avg, label, mean))

  # order dataframe levels according to average.
  label_levels <- names(simulation_avgs[order(simulation_avgs)])
  sims_dataframe$label <- factor(sims_dataframe$label, levels = label_levels)



  # Load observed pairwise_difference data.
  pwd_data <- results_data[which(results_data$Pair_name == pair), ]
  # extract mean and 95%CI.
  plotdist_mean <- pwd_data$Corr.Avg.PWD[1]
  plotdist_std  <- pwd_data$Corr.CI.95[1]
  # ---- Prepare observed pwd lines.
  fin_rel <- list(type = "line",
                 line = list(color = "black"),
                 xref = "paper",
                 yref = "y",
                 x0   = 0,
                 x1   = 1,
                 size = 2
  )
  sd_plus <- list(type = "line",
                  line = list(color = "black",
                              dash  = "dash"
                  ),
                  xref = "paper",
                  yref = "y",
                  x0   = 0,
                  x1   = 1,
                  size = 0.3
  )
  sd_min  <- list(type = "line",
                  line = list(color = "black",
                              dash  = "dash"
                  ),
                  xref = "paper",
                  yref = "y",
                  x0   = 0,
                  x1   = 1,
                  size = 0.3
  )
  uncert  <- list(type      = "rect",
                  fillcolor = "#AAAAAA",
                  line      = "#AAAAAA",
                  opacity   = 0.3,
                  x0        = 0,
                  x1        = 1,
                  xref      = "paper",
                  y0        = plotdist_mean - plotdist_std,
                  y1        = plotdist_mean + plotdist_std,
                  yref      = "y"
  )

  fin_rel[c("y0", "y1")] <- plotdist_mean
  sd_plus[c("y0", "y1")] <- plotdist_mean + plotdist_std
  sd_min[c("y0",  "y1")] <- plotdist_mean - plotdist_std
  lines <- list(fin_rel, sd_plus, sd_min, uncert)



colorpalette <- viridis::viridis_pal(
    option    = "D",
    direction = -1
  )(length(labels_to_plot))


  # ---- A. Main plot.
  plotly::plot_ly(type     = "violin",
                  data     = sims_dataframe,
                  x        = ~label,
                  y        = ~avg,
                  color    = ~label,
                  colors   = colorpalette,
                  points   = "outliers",
                  box      = list(visible = TRUE, notched = TRUE, notchwidth = 0.5),
                  meanline = list(visible = TRUE)
                 ) %>%

  # ---- B. Add observed avg and uncertainty annotations
  plotly::layout(title  = list(text = paste("Simulation results of", pair),
                               y    = 0.99,
                               yref = "paper"
                              ),
                xaxis  = list(
                  title          = "Pairwise comparisons",
                  showticklabels = FALSE
                ),
                yaxis  = list(title = "Mean genetic distance"),
                legend = list(
                  title       = list(text = "<b>Label</b>"),
                  orientation = "h",
                  x = 0,
                  y = -0.2,
                  xref = "paper",
                  yref = "paper"
                ),
                shapes = lines
                ) %>%

  plotly::config(editable    = TRUE,
                 displaylogo = FALSE,
                 scrollZoom  = TRUE,
                 toImageButtonOptions = list(
                   format = "svg",
                   filename = paste0("simulation_plots-", pair)
                 )
                )
}
