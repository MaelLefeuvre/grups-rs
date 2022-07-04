#' @export 
#' @import shiny
#' @param values vector of checkbox elements.
#' @param title  title of the div / section.
#' @param input_id html/css id
#' @param ncol number of display columns.
#' @return list containing shiny css style, header, and checkboxgroupinput.
shiny_reactive_checkbox <- function(values, title, input_id, ncol = 2) {

  css_style <- paste0(
    "column-count: ", ncol, ";
    -webkit-column-count: ", ncol, "; /* Chrome, Safari, Opera */ 
    -moz-column-count: ", ncol, ";    /* Firefox */ 
    -moz-column-fill: auto;
    -column-fill: balance;
    "
  )

  list(
    shiny::h5(title),
    shiny::tags$div(
      align = "justify",
      style = css_style,
      shiny::checkboxGroupInput(
        inputId  = input_id,
        label    = NULL,
        selected = values,
        inline   = FALSE,
        values
      )
    )
  )
}