#' Main shiny dashboard application
#' @import shiny
#' @import plotly
#' @import stringr
#' @export
app <- function(ui, server, data_dir = "./grups_output", ...) {
  #pwd_path <- "../Esperstedt-output/Esperstedt.q30Q30.g1k-phase3-variant_sites.chr1-22.only_SNPs-only-tv.pwd"
  #data_dir <- "../Esperstedt-output/"

  # Search for pwd_file(s)
  pwd_files <- list.files(
    path = data_dir,
    full.names = TRUE,
    pattern = "\\.pwd$"
  )

  # Search for blk_file(s)
  blk_files <- list.files(
    path = paste(data_dir, "blocks", sep = "/"),
    full.names = TRUE,
    pattern = "\\.blk$"
  )

  ui <- shiny::fluidPage(
    theme = bslib::bs_theme(bootswatch = "darkly"),
    shiny::titlePanel("GRUPS-plots"),
    shiny::navbarPage("GRUPS-plots",
      # ---- Render pwd barplots tab.
      shiny::tabPanel("PWD-barplot",
        shiny::uiOutput("pwd_barplots")
      ),

      # ---- Render block scatterplots tab
      shiny::tabPanel("Blocks",
        shiny::uiOutput("blk_scatterplots")
      )
    )
  )
  server <- function(input, output, session) {

    # ---- 1. Render all pwd-barplots
    # See: https://gist.github.com/wch/5436415/
    #      https://stackoverflow.com/questions/22823843/shiny-r-print-multiple-plots-using-a-loop
    output$pwd_barplots <- shiny::renderUI({
      plot_output_list <- lapply(seq_along(pwd_files), function(i) {
        plotname <- paste("pwd_barplot", i, sep = "_")
        plotly::plotlyOutput(plotname)
      })

      base::do.call(shiny::tagList, plot_output_list)
    })

    for (i in seq_along(pwd_files)) {
      local({
        my_i <- i
        plotname <- paste("pwd_barplot", my_i, sep = "_")

        output[[plotname]] <- plotly::renderPlotly({
          grups.plots::plot_pairwise_diff(path = pwd_files[my_i])
        })
      })
    }

    # ---- 2. Render all blk-scatterplots
    blk_pairs <- stringr::str_extract(
      blk_files,
      "(?<=-)([^-]+)-([^-]+)(?=.blk$)"
    )
    output$blk_scatterplots <- shiny::renderUI({
      plot_output_list <- lapply(seq_along(blk_files), function(i) {
        plotname  <- paste("blk_scatterplot", blk_pairs[i], sep = "_")
        print(plotname)
        shiny::tabPanel(
          blk_pairs[i],
          plotly::plotlyOutput(plotname)
        )
      })
      base::do.call(shiny::tagList, plot_output_list)
    })

    for (i in seq_along(blk_files)) {
      local({
        my_i <- i
        plotname  <- paste("blk_scatterplot", blk_pairs[my_i], sep = "_")
        output[[plotname]] <- plotly::renderPlotly({
          grups.plots::plot_sliding_window(
            grups.plots::load_blockfile(blk_files[my_i])
          )
        })
      })
    }


  }

  shiny::shinyApp(ui, server, ...)
}