#' @export
get_bc_matrix <- function(sims_data, labels_to_keep) {
  labels_relationships <- levels(sims_data$label)

  bc_matrix <- matrix(
    data = 0,
    nrow = length(labels_relationships),
    ncol = length(labels_relationships)
  )
  colnames(bc_matrix) <- rownames(bc_matrix) <- labels_relationships

  lapply(labels_relationships,
    FUN = function(rel1) {
      lapply(labels_relationships,
        FUN = function(rel2) {
          rel1_rows <- which(sims_data$label == rel1)
          rel2_rows <- which(sims_data$label == rel2)
          merged_data <- c(
            sims_data[rel1_rows, ]$avg,
            sims_data[rel2_rows, ]$avg
          )

          merged_data_hist  <- hist(
            x      = merged_data,
            breaks = length(merged_data) / 10,
            plot   = FALSE
          )
          x_hist <- hist(
            x      = sims_data[rel1_rows, ]$avg,
            breaks = merged_data_hist$breaks,
            plot   = FALSE
          )
          y_hist <- hist(
            x      = sims_data[rel2_rows, ]$avg,
            breaks = merged_data_hist$breaks, plot = FALSE
          )
          x_hist_prob <- x_hist$counts / sum(x_hist$counts)
          y_hist_prob <- y_hist$counts / sum(y_hist$counts)
          n_bc_bins <- length(x_hist_prob)
          # estimate Bhattacharyya co-efficient
          for (i in seq(1, n_bc_bins)){
            bc_matrix[rel1, rel2] <<- bc_matrix[rel1, rel2] +
                                      sqrt(x_hist_prob[i] * y_hist_prob[i])
          }
        }
      )
    }
  )

  # Filter-out unwanted relationships
  bc_matrix <- bc_matrix[which(rownames(bc_matrix) %in% labels_to_keep), ]
  bc_matrix <- bc_matrix[, which(colnames(bc_matrix) %in% labels_to_keep)]

  bc_matrix
}