#' @import e1071
#' @import future.apply
#' @import future
#' @param sim_files a list of '.sim' files for each
#' @param results_file the main '.result' file of GRUPS
#' @param progressor optional progressr::progressor to report progress.
#' @export
get_svmop_probs <- function(
  results_file,
  sim_files,
  threads    = 1,
  progressor = NULL
) {
  future::plan(future::multisession, workers = threads)
  probs <- data.frame(
    Pair_name    = results_file$Pair_name,
    Corr.Avg.PWD = results_file$Corr.Avg.PWD
  )

  i <- 0
  svm_probs <- future.apply::future_sapply(
    X           = results_file$Pair_name,
    future.seed = TRUE,
    FUN         = function(pair_name) {
      if (!is.null(progressor)) {
        i <<- i + 1
        msg <- sprintf("Processing [%d/%d]: %s", i, NROW(sim_files), pair_name)
        progressor(msg)
      }
      svms     <- fit_svms(sim_files[pair_name, ])
      pair_row <- which(probs$Pair_name == pair_name)
      get_within_class_probs(svms, probs[pair_row, ])
    }
  )
  probs <- merge(probs, t(svm_probs), by.x = "Pair_name", by.y = "row.names")

  probs
}

fit_svms <- function(sim_file) {
  train <- grups.plots::load_simfile(sim_file)

  ## Separate in test/train
  #all_aboard_the_train   <- sample(nrow(data), prop * nrow(data))
  #train <- data[all_aboard_the_train,]
  #test  <- data[-all_aboard_the_train,]
  #
  ## Ensure we have the correct number of data points
  #stopifnot(nrow(data) == nrow(train)+nrow(test))

  svmop_fits <- list()
  for (rel_index in 1:(nlevels(train$label) - 1)) {
    colname     <- levels(train$label)[rel_index]
    svmop_train <- data.frame(
      label          = as.numeric(train$label) > rel_index,
      Corr.Avg.PWD   = train$avg
    )

    # We can assign a different cost according to the relatedness
    # Let's consider it twice as worse, the more we're checking for fine grain
    # relatedness
    #cost <- 1000 * 2^(nlevels(train$label) - rel_index - 1)
    #print(paste("Rel:", colname, "cost:", cost))
    svmop_fits[[colname]] <- e1071::svm(label ~ Corr.Avg.PWD,
      data        = svmop_train,
      probability = TRUE,
      kernel      = "linear",
      type        = "C-classification"
    )
  }
  list(fits = svmop_fits, rel_order =  levels(train$label))
}

# svmop: a list of Ordinally partitionned SVM classifiers
# return the list of being within an ordinal class.
get_within_class_probs <- function(svmops, obs) {
  # Get the probability of being greater than a given class
  # class: character string
  # fits: a list of SVM fits. One for each class
  # obs: a dataframe with an 'avg' column.
  prob_greater_than <- function(class, fits, obs) {
    pred <- predict(fits[[class]], obs, probability = TRUE)
    attr(pred, "probabilities")[, c("TRUE")]
  }

  probs <- list()

  # First class: P(Y=1) = 1 - P(Y>1)
  first_rel          <- svmops$rel_order[1]
  probs[[first_rel]] <- 1 - prob_greater_than(first_rel, svmops$fits, obs)

  # Intermediate classes: P(Y=n) = P(Y>n-1) - P(Y>n)
  for (i in 2:(length(svmops$fits))) {
    rel_n          <- svmops$rel_order[i]
    rel_prev       <- svmops$rel_order[i - 1]
    probs[[rel_n]] <- prob_greater_than(rel_prev, svmops$fits, obs) -
                      prob_greater_than(rel_n,    svmops$fits, obs)
  }

  # Final class: P(Y=N) = P(Y > N-1)
  final_rel          <- rev(svmops$rel_order)[1]
  previous_rel       <- rev(svmops$rel_order)[2]
  probs[[final_rel]] <- prob_greater_than(previous_rel, svmops$fits, obs)

  probs
}