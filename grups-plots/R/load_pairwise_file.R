#' @export
#' @importFrom utils read.table
#' @param path path leading to a GRUPS `.pwd` results file.
#' @return dataframe
load_pairwise_file <- function(
  path,
  res_data,
  norm_avg_type,
  sample_regex,
  min_overlap  = 0,
  norm_request = FALSE,
  norm_method  = "Raw",
  norm_metric  = "Pairwise-comparisons",
  norm_value   = NULL
) {
  # Load dataset
  pwd_data <- read.table(path, sep = "\t", header = TRUE)

    subset_cols <- c(
      "Pair_name", "Corr.Overlap", "Corr.Avg.PWD", "Corr.Sum.PWD",
      "Corr.Avg.PWD", "Corr.CI.95", "Corr.Avg.Phred"
    )
    pwd_data <- merge(
      pwd_data,
      res_data[, subset_cols],
      by.x = "Pair_name",
      by.y = "Pair_name"
    )

  norm_avg_col     <- paste(norm_avg_type, "Avg.PWD", sep = ".")
  norm_overlap_col <- paste(norm_avg_type, "Overlap", sep = ".")
  norm_ci_col      <- paste(norm_avg_type, "CI.95", sep = ".")

  # Filter out individuals with an overlap lower than the req. treshold
  pwd_data <- pwd_data[which(pwd_data[[norm_overlap_col]] >= min_overlap), ]



  pwd_data <- pwd_data[order(pwd_data[[norm_avg_col]]), ]

  # Comparisons are considered as self-comparisons if both individuals share
  # the same name. Determined by searching for a name pattern & ensuring:
  # - the name pattern is found twice
  # - both matches are only separated by a dash
  # - the match spans the entire string.
  # i.e. : ^([A-Za-z0-9]+([-0-9]+){0,1})-(\1)$
  twin_regex <- paste0("^(", sample_regex, ")-(\\1)$")
  pwd_data$Self <- stringr::str_detect(
    pwd_data$Pair_name,
    twin_regex
  )

  # Order pairs according to their avg pwd.
  pwd_data$Pair_name <- factor(
    pwd_data$Pair_name,
    levels = pwd_data$Pair_name[
      order(pwd_data[[norm_avg_col]], decreasing = FALSE)
    ]
  )

  # Format for raw display
  pwd_data$Self <- as.character(pwd_data$Self)

  # Normalize values
  norm_method_function <- grups.plots::get_norm_method(
    norm_request,
    norm_method,
    norm_metric,
    norm_value
  )
  pwd_data$Norm.Avg <- pwd_data[[norm_avg_col]] / norm_method_function(pwd_data)
  pwd_data$Norm.CI  <- pwd_data[[norm_ci_col]]  / norm_method_function(pwd_data)

  norm_values <- grups.plots::get_norm_values(
    norm_method,
    norm_request,
    norm_metric,
    norm_value,
    pwd_data
  )


  # Assign least z-score + putative relationship
  z_scores  <- lapply(pwd_data$Norm.Avg, FUN = function(x) {
    if (is.na(x)) return(NA) # Deal with NaN (when overlap == 0)
    z_temps <- abs(x - unlist(norm_values))
    z_temps <- z_temps[which.min(z_temps)]
    z_temps
  })

  pwd_data$Z_score <- unlist(z_scores)

  pwd_data$Rel <- factor(
    names(unlist(z_scores)),
    levels = c("Unrelated", "First", "Second", "Third", "Fourth", "Fifth",
               "Self"
              )
  )

  list(data = pwd_data, norm_values = norm_values)
}
