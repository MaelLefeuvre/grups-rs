#' @export
#' @import tidyr
#' @import reshape2
#' @param data A GRUPS Results file, containing pair of individuals + their most_likely relationship
#' @return A Kinship Matrix
get_kinship_matrix <- function(data, order) {
  data.frame(pair = data$Pair_name, rel= data$Most_Likely_rel) %>% tidyr::separate(
    col   = pair,
    into  = c("Left_ind", "Right_ind"),
    sep   = "-",

  ) #%>% reshape2::dcast(formula = Left_ind~Right_ind, fun.aggregate = NULL, value.var = "rel")
}