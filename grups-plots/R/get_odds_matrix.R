

#' @export
get_odds_matrix <- function(sims_data, observed_results, pair, labels_to_keep) {
  labels_relationships <- levels(sims_data$label)

  obs_pwd <- observed_results[
    which(observed_results$Pair_name == pair),
  ]$Corr.Avg.PWD

  # Compute a matrix of z-scores. (@TODO this computation is all over the place)
  obs_dist_z <- rep(NA, length(labels_relationships))
  obs_prob   <- rep(NA, length(labels_relationships))

  #MODIFIED: Median Average Deviation. More robust than Z-score
  mad        <- rep(NA, length(labels_relationships))

  names(obs_dist_z)  <- labels_relationships
  names(obs_prob)    <- labels_relationships
  names(mad)         <- labels_relationships     #MODIFIED:

  for (rel in labels_relationships) {
    relrows           <- which(sims_data$label == rel)
    rel_avg           <- sims_data[relrows, ]$avg
    obs_dist_z[[rel]] <- abs(obs_pwd - mean(rel_avg))   / sd(rel_avg)
    mad[[rel]]        <- abs(obs_pwd - median(rel_avg)) / mad(rel_avg) #MODIFIED
    obs_prob[[rel]]   <- pnorm(-obs_dist_z[[rel]])
  }

  # Compute the Odds ratio matrix for each relationship.
  or_matrix <- matrix(
    data = NA,
    nrow = length(labels_relationships),
    ncol = length(labels_relationships)
  )
  colnames(or_matrix) <- rownames(or_matrix) <- labels_relationships

  lapply(labels_relationships, FUN = function(rel1) {
    lapply(labels_relationships, FUN = function(rel2) {
      or_matrix[[rel1, rel2]] <<- (obs_prob[[rel1]] / (1 - obs_prob[[rel1]])) /
                                   (obs_prob[[rel2]] / (1 - obs_prob[[rel2]]))
    })
  })

  #//best_odds <- obs_prob[which(min(obs_dist_z) == obs_dist_z)] /
  #//            (1 - obs_prob[which(min(obs_dist_z) == obs_dist_z)])
  #//for (t in seq_along(labels_relationships)){
  #//  if (t != which(min(obs_dist_z) == obs_dist_z)) {
  #//    these_odds <- obs_prob[t] / (1 - obs_prob[t])
  #//    cat(paste(
  #//      labels_relationships[t],
  #//      best_odds / these_odds,
  #//      "\n",
  #//      sep = "\t"
  #//    ))
  #//  }
  #//}

  # Filter-out unwanted relationships
  or_matrix <- or_matrix[which(rownames(or_matrix) %in% labels_to_keep), ]
  or_matrix <- or_matrix[, which(colnames(or_matrix) %in% labels_to_keep)]

  # Odds of belonging to the same relationship are always equal to 1
  or_matrix[col(or_matrix) == row(or_matrix)] <- 1

  log(or_matrix)
}
