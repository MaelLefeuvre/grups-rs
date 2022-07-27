

#' @export
get_odds_matrix <- function(sims_data, observed_results, pair, labels_to_keep) {
  labels_relationships <- levels(sims_data$label)

  obs_pwd = observed_results[which(observed_results$Pair_name == pair),]$Corr.Avg.PWD

  # Compute a matrix of z-scores. (this computation is all over the place)
  obsDistZ = rep(NA, length(labels_relationships))
  names(obsDistZ) <- labels_relationships
  obsProb = rep(NA, length(labels_relationships))
  names(obsProb) <- labels_relationships
  for (rel in labels_relationships) {
    relrows         <- which(sims_data$label == rel)
    rel_avg         <- sims_data[relrows,]$avg
    obsDistZ[[rel]] <- abs(obs_pwd - mean(rel_avg)) / sd(rel_avg)
    obsProb[[rel]]  <- pnorm(-obsDistZ[[rel]])

    print(paste("Rel: ", rel, "| Rel avg:", mean(rel_avg), sd(rel_avg), "| obsDistZ:", obsDistZ[[rel]], "| obsProb:", obsProb[[rel]]))
  }

  print("Z-scores")
  print(obsDistZ)
  cat("\n")
  print("Probability of observation within simulated relationship distribution")
  print(obsProb)

  # Compute the Odds ratio matrix for each relationship.
  ORs_matrix = matrix(NA, length(labels_relationships), length(labels_relationships))
  colnames(ORs_matrix) <- rownames(ORs_matrix) <- labels_relationships

  lapply(labels_relationships, FUN = function(rel1) {
    lapply(labels_relationships, FUN=function(rel2) {
      ORs_matrix[[rel1, rel2]] <<- (obsProb[[rel1]]/(1-obsProb[[rel1]])) / (obsProb[[rel2]]/(1-obsProb[[rel2]]))
    })
  })

  #cat("Odds ratios matrix\n")
  #print(ORs_matrix)

  print(paste("Most likely relationship Z-score:", obsDistZ[which(min(obsDistZ) == obsDistZ)], labels_relationships[which(min(obsDistZ) == obsDistZ)], sep=" "))
  print(paste("Odds ratio of obs within most likely relationship (", labels_relationships[which(min(obsDistZ) == obsDistZ)], ") vs. other relationship", sep=""))
  best_odds = obsProb[which(min(obsDistZ) == obsDistZ)]/(1-obsProb[which(min(obsDistZ) == obsDistZ)])
  for (t in 1:length(labels_relationships)){
      if (t != which(min(obsDistZ) == obsDistZ)){
          these_odds = obsProb[t]/(1-obsProb[t])
          cat(paste(labels_relationships[t], best_odds/these_odds, "\n", sep="\t"))
      }
  }

  # Filter-out unwanted relationships
  ORs_matrix <- ORs_matrix[which(rownames(ORs_matrix) %in% labels_to_keep),]
  ORs_matrix <- ORs_matrix[, which(colnames(ORs_matrix) %in% labels_to_keep)]

  log(ORs_matrix)
}
