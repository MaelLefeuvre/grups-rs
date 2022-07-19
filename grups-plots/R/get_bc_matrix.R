#' @export
get_bc_matrix <- function(sims_data, labels_to_keep) {
  labels_relationships <- levels(sims_data$label)

  BC_matrix = matrix(0, length(labels_relationships), length(labels_relationships))
  colnames(BC_matrix) <- rownames(BC_matrix) <- labels_relationships

  lapply(labels_relationships,
    FUN=function(rel1) {
      lapply(labels_relationships,
        FUN = function(rel2) {
          rel1_rows <- which(sims_data$label == rel1)
          rel2_rows <- which(sims_data$label == rel2)
          merged_data <- c(sims_data[rel1_rows,]$avg, sims_data[rel2_rows,]$avg)

          merged_data_hist  <- hist(merged_data, breaks=length(merged_data)/10, plot = FALSE)
          histX <- hist(sims_data[rel1_rows, ]$avg, breaks=merged_data_hist$breaks, plot = FALSE)
          histY <- hist(sims_data[rel2_rows, ]$avg, breaks=merged_data_hist$breaks, plot = FALSE)
          histXPr <- histX$counts / sum(histX$counts)
          histYPr <- histY$counts / sum(histY$counts)
          numBCbins <- length(histXPr)
          # estimate Bhattacharyya co-efficient
          for (i in seq(1, numBCbins)){
            BC_matrix[rel1, rel2] <<- BC_matrix[rel1, rel2] + sqrt(histXPr[i]*histYPr[i])
          }
        }
      )
    }
  )

  # Filter-out unwanted relationships
  BC_matrix <- BC_matrix[which(rownames(BC_matrix) %in% labels_to_keep),]
  BC_matrix <- BC_matrix[, which(colnames(BC_matrix) %in% labels_to_keep)]
  BC_matrix
}