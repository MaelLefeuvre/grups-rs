#' Main shiny dashboard application
#' @import shiny
#' @import plotly
#' @import stringr
#' @export
app <- function(ui, server, data_dir = "./grups_output", ...) {
  #pwd_path <- "../Esperstedt-output/Esperstedt.q30Q30.g1k-phase3-variant_sites.chr1-22.only_SNPs-only-tv.pwd"
  #data_dir <- "../Esperstedt-output/"

  # ---- 1. Search for pwd_file(s)
  pwd_files <- list.files(
    path = data_dir,
    full.names = TRUE,
    pattern = "\\.pwd$"
  )

  # ---- Sanity checks
  shiny::validate(
    shiny::need(
      length(pwd_files) == 1,
      "[ERROR]: Exactly one `.pwd` file must exist within `data_dir`. Exiting."
    )
  )

  # ---- 2a. Search for blk_file(s)                                          [A FUNCTION]
  blk_files <- list.files(
    path = paste(data_dir, "blocks", sep = "/"),
    full.names = TRUE,
    pattern = "\\.blk$"
  )

  # ---- 2b. Extract block pair names, parse all that data into a df.        [B FUNCTION]
  blk_files <- data.frame(
    path      = blk_files,
    row.names = stringr::str_extract(  # Extract pair names
      blk_files,
      "(?<=-)([^-]+)-([^-]+)(?=.blk$)"
    ),
    stringsAsFactors = FALSE
  )
  # ---- 3a. Search for simulation files                                     [A FUNCTION]
  sim_files <- list.files(
    path = paste(data_dir, "simulations", sep = "/"),
    full.names = TRUE,
    pattern = "\\.sims$"
  )
  
  # ---- 3b. Extract simulations. pair names, parse everything into a df.    [B FUNCTION]
  sim_files <- data.frame(
    path = sim_files,
    row.names = stringr::str_extract( # Extract pair names
      sim_files,
      "(?<=-)([^-]+)-([^-]+)(?=.sims$)"
    ),
    stringsAsFactors = FALSE
  )

  ui <- shiny::fluidPage(
    theme = bslib::bs_theme(bootswatch = "darkly"),
    ##shiny::titlePanel("GRUPS-plots"),
    shiny::navbarPage("GRUPS-plots",
      # ---- 1. Render pwd barplots tab.
      shiny::tabPanel("PWD-barplot",
        shiny::uiOutput("pwd_barplots")
      ),

      # ---- 2. Render block scatterplots tab
      shiny::tabPanel("Sliding mismatch rate",
        shiny::fluidPage(
          shiny::sidebarLayout(
            shiny::sidebarPanel(
              shiny::selectInput("block_pair",
                "Pair name",
                rownames(blk_files)
              ),
              shiny::sliderInput("block_width",
                "Block-width",
                value = 20,
                min   = 2,
                max   = 100
              ),
              shiny::sliderInput("block_step",
                "Sliding window step",
                value = 1,
                min = 1,
                max = 10
              )
            ),
            shiny::mainPanel(
              shiny::tabsetPanel(
                shiny::tabPanel("Plot",
                  plotly::plotlyOutput("block_scatterplot")
                ),
                shiny::tabPanel("Raw data",
                  shiny::tableOutput("block_dataframe")
                )
              )
            )
          )
        )
      ),

      # ---- 3. Render violin plots
      shiny::tabPanel("Simulations plot",
        shiny::fluidPage(
          shiny::sidebarLayout(
            shiny::sidebarPanel(
              shiny::selectInput("sim_pair",
                "Pair name",
                rownames(sim_files)
              )
            ),
            shiny::mainPanel(
              shiny::tabsetPanel(
                shiny::tabPanel("Plot",
                  plotly::plotlyOutput("sims_violinplot")
                ),
                shiny::tabPanel("Raw data",
                  shiny::tableOutput("sims_dataframe")
                )
              )
            )
          )

        )
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

    # ---- 2a. Load / Update block dataframe
    load_block_dataframe <- shiny::reactive(
      grups.plots::load_blockfile(
        path  = blk_files[input$block_pair, ],
        width = input$block_width,
        step    = input$block_step
      )
    )

    # ---- 2b. Render block-scatterplot
    output$block_scatterplot <- plotly::renderPlotly({
      grups.plots::plot_sliding_window(
        load_block_dataframe(),
        input$block_pair
      )
    })

    # ---- 2c. Render block dataframe
    output$block_dataframe <- shiny::renderTable(
      load_block_dataframe()
    )

    # ---- 3a. Load / Update simulations dataframe
    load_sims_dataframe <- shiny::reactive(
      grups.plots::load_simfile(
        path = sim_files[input$sim_pair, ]
      )
    )

    # ---- 3b. Render simulations violin plot
    output$sims_violinplot <- plotly::renderPlotly({
      grups.plots::plot_pedigree_sims(
        sims_dataframe = load_sims_dataframe(),
        pwd_path       = pwd_files[1],
        pair           = input$sim_pair
      )
    })

    # ---- 2c. Render simulations dataframe
    output$sims_dataframe <- shiny::renderTable(
      load_sims_dataframe()
    )

  }

  shiny::shinyApp(ui, server, ...)
}