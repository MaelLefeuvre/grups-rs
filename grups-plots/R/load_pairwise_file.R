#' @export
#' @importFrom utils read.table
#' @param path path leading to a GRUPS `.pwd` results file.
#' @return dataframe
load_pairwise_file <- function(path, min_overlap = 0, norm_method, norm_metric, norm_value) {
  # Load dataset
  pwd_data <- read.table(path, sep = "\t", header = TRUE)

  # Filter out individuals with an overlap lower than the req. treshold
  pwd_data <- pwd_data[which(pwd_data[, 2] >= min_overlap),]

  pwd_data <- pwd_data[order(pwd_data[, 4]), ]

  # Companion dataset, ordered according to avg.pwd
  plot_data <- data.frame(
    avg     = pwd_data[, 4],
    ci      = pwd_data[, 5],
    pairs   = pwd_data[, 1],
    overlap = pwd_data[, 2]
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

  # Format for raw display
  plot_data$self <- as.character(plot_data$self)
  plot_data[, c(3, 4, 1, 2)]

  # Normalize values
  norm_method_function <- grups.plots::get_norm_method(norm_method, norm_metric, norm_value)

  plot_data$norm_avg <- plot_data$avg / norm_method_function(plot_data)
  plot_data$norm_ci  <- plot_data$ci  / norm_method_function(plot_data)

  norm_values <- grups.plots::get_norm_values(norm_method, norm_metric, norm_value, plot_data)

  # Assign least z-score + putative relationship
  z_scores <- lapply(plot_data$norm_avg, FUN = function(x) {
    z_temps <- abs(x - unlist(norm_values))
    z_temps <- z_temps[which.min(z_temps)]
  })
  plot_data$z_score <- unlist(z_scores)

  plot_data$rel <- factor(names(unlist(z_scores)), levels = c("Unrelated", "First", "Second", "Third", "Fourth", "Fifth", "Self"))

  list(data = plot_data, norm_values = norm_values)
}