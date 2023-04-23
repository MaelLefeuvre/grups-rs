#' @export
#' @import tidyr
#' @import reshape2
#' @param data GRUPS Results file, w/ pairs of individuals & their assigned rel
#' @return A Kinship Matrix
get_kinship_matrix <- function(data, sample_regex) {
  data.frame(
    pair = data$Pair_name,
    rel  = data$Most_Likely_rel
  ) %>% tidyr::separate_wider_regex(
    col   = pair,
    patterns = c(Left_ind = sample_regex, Right_ind = sample_regex),
  )
}