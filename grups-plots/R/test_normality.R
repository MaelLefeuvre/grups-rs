#' @export
test_normality <- function(sims_file, alpha = 0.05) {
  labels_relationships <- levels(sims_file$label)

  kstest <- rep(NA, length(labels_relationships))
  names(kstest) <- labels_relationships

  data <- sims_file$avg
  for (rel in labels_relationships) {
    relrows <- which(sims_file$label == rel)
    kstest[[rel]] <- ks.test(
      x = data[relrows],
      y = rnorm(
        length(data[relrows]),
        mean = mean(data[relrows]),
        sd = sd(data[relrows])
      )
    )$p.val
  }

  t(data.frame(p.val = kstest, reject = as.character(kstest < alpha)))
}