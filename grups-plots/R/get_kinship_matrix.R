#' @export
#' @import tidyr
#' @import reshape2
#' @param data GRUPS Results file, w/ pairs of individuals & their assigned rel
#' @return A Kinship Matrix
get_kinship_matrix <- function(data, order) {
  data.frame(
    pair = data$Pair_name,
    rel  = data$Most_Likely_rel
  ) %>% tidyr::separate(
    col   = pair,
    into  = c("Left_ind", "Right_ind"),
    sep   = "-",
  )
}