#' @export
#' @importFrom dplyr %>%
#' @importFrom utils read.table
#' @import plotly
#' @param path path leading to a GRUPS `.pwd` results file.
#' @return plotly barplot
plot_pairwise_diff <- function(path) {
  # Load dataset
  pwd_data <- read.table(path, sep = "\t", header = TRUE)

  # Companion dataset, ordered according to avg.pwd
  plot_data <- data.frame(
    avg   = sort(pwd_data[, 4]),
    ci    = pwd_data[, 5][order(pwd_data[, 4])],
    pairs = pwd_data[, 1][order(pwd_data[, 4])]
  )

  # Comparisons are considered as self-comparisons if both individuals share
  # the same name.
  plot_data$self <- lapply(
    strsplit(as.character(plot_data$pairs), "-"),
    FUN = function(x) x[1] == x[2]
  )

  # Order pairs according to their avg pwd.
  plot_data$pairs <- factor(
    plot_data$pairs,
    levels = plot_data$pairs[order(plot_data$avg, decreasing = FALSE)]
  )

  # Compute mean of avg pwd for self-comparisons.
  # - [WARN]: NaN if there are no self-comparisons.
  avg_self <- mean(plot_data$avg[which(plot_data$self == TRUE)])

  # Add horizontal Ms lines as annotations if avg_self is not undefined
  if (!is.nan(avg_self)) {
    plot_annotations <- list(
      grups.plots::ms_annotation("Ms", avg_self),
      grups.plots::ms_annotation("(3/2)Ms", (3 / 2) * avg_self),
      grups.plots::ms_annotation("2Ms", 2 * avg_self)
    )
    plot_shapes <- list(
      grups.plots::hline(avg_self),
      grups.plots::hline((3 / 2) * avg_self),
      grups.plots::hline(2 * avg_self)
    )
  } else {
    plot_annotations <- plot_shapes <- list()
  }

  # Plot
  plotly::plot_ly(type    = "bar",
                  data    = plot_data,
                  x       = ~pairs,
                  y       = ~avg,
                  error_y = ~list(array = ci,
                                  color = "#000000"
                 )
  ) %>%
  plotly::layout(title = list(text = "Raw average genetic distances.",
                              y    = 0.99,
                              yref = "paper"
                             ),
                 yaxis = list(title = "Mean genetic distance",
                              range = c(0, 0.300)
                             ),
                 xaxis       = list(title = "Pairs"),
                 shapes      = plot_shapes,
                 annotations = plot_annotations,
                 margin      = list(r = 60)
  ) %>%
  plotly::config(editable = TRUE, displaylogo = FALSE, scrollZoom = TRUE)
}
