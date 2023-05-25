#' @export
#' @importFrom utils read.table
#' @param path path leading to a GRUPS .sims results file
#' @return dataframe containing simulation results for a given pair.
load_simfile <- function(path) {
  data           <- read.table(path, sep = "\t", header = FALSE)
  colnames(data) <- c(
    "replicate",
    "label",
    "parent0",
    "parent1",
    "parent0.id",
    "parent1.id",
    "pwd",
    "overlap",
    "avg"
  )

  # Reorder labels according to distribution average.
  data$label <- factor(data$label)

  tryCatch(
    expr = {
      data$label <- factor(
        data$label,
        levels = levels(data$label)[
          order(aggregate(avg ~ label, data, FUN = mean)$avg)
        ]
      )
    },
    error = function(cond) {
      msg <- paste(
        "Failed to aggregate relatedness labels according to average PWD",
        "values in simfile '", path, "'. This might be caused by a complete",
        "lack of overlap between the two individuals."
      )
    }
  )
  data
}