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
    shiny::navbarPage("GRUPS-plots",

      # ---- 1. Render pwd barplots tab.
      shiny::tabPanel("Raw Pairwise differences.",
        shiny::fluidPage(
          shiny::sidebarLayout(
            shiny::sidebarPanel(
              shiny::fluidRow(
                shiny::column(5,
                  shiny::radioButtons("norm_method",
                    "Normalization method",
                    c("Raw", "Self-comparisons", "All", "Value"),
                    selected = "Raw",
                    width = "100%"
                  )
                ),
                shiny::column(5,
                  shiny::conditionalPanel(
                    condition = "input.norm_method != 'Raw' && input.norm_method != 'Value'",
                    shiny::radioButtons("norm_metric",
                      "Normalization metric",
                      c("Median", "Mean", "Min", "Max")
                    )
                  ),
                  shiny::conditionalPanel(
                    condition = "input.norm_method == 'Value'",
                    shiny::numericInput("norm_value",
                      "Normalization value",
                      0.25,
                      min = 0.01,
                      max = 2,
                      step = 0.01
                    )
                  )
                )
              ),
              shiny::numericInput("min_overlap",
                "Minimum SNP overlap treshold",
                0,
                min = 0,
                max = .Machine$integer.max,
                step = 10000
              ),
              shiny::checkboxInput("hide_self_comparisons",
                "Hide self-comparisons",
              )
            ),
            shiny::mainPanel(
              shiny::tabsetPanel(
                shiny::tabPanel("Plot",
                  plotly::plotlyOutput("pwd_barplot")
                ),
                shiny::tabPanel("Raw data",
                  shiny::tableOutput("pwd_dataframe")
                )
              )
            )
          )
        )
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
              ),
              shiny::uiOutput("chromosome_subset_checkboxgroup")
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
              ),
              shiny::uiOutput("simulation_labels_checkboxgroup"),
              shiny::numericInput("ks_alpha", "Kolmogorov-Smirnov alpha treshold", min=0, max=1, step=0.01, value=0.05)
            ),
            shiny::mainPanel(
              shiny::tabsetPanel(
                shiny::tabPanel("Plot",
                  plotly::plotlyOutput("sims_violinplot"),
                  DT::dataTableOutput("ks_normality_test"),
                  shiny::tableOutput("bc_matrix")
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

    # ---- 1a. Load / Update pairwise dataframe.
    load_pairwise_dataframe <- shiny::reactive(
      grups.plots::load_pairwise_file(
        path = pwd_files[1],
        min_overlap = input$min_overlap,
        norm_method = input$norm_method,
        norm_metric = input$norm_metric,
        norm_value = input$norm_value
      )
    )

    # ---- 1b. Render pairwise difference barplot.
    output$pwd_barplot <- plotly::renderPlotly(
      grups.plots::plot_pairwise_diff(
        data = load_pairwise_dataframe(),
        hide_self_comparisons = input$hide_self_comparisons,
        norm_method = input$norm_method,
        norm_metric = input$norm_metric,
      )
    )

    # ---- 1c. Render pairwise dataframe
    output$pwd_dataframe <- shiny::renderTable(
      load_pairwise_dataframe()$data
    )

    # ---- 2a. Load / Update block dataframe
    load_block_dataframe <- shiny::reactive(
      grups.plots::load_blockfile(
        path  = blk_files[input$block_pair, ],
        width = input$block_width,
        step  = input$block_step
      )
    )

    # ---- 2b. Output chromosome subset checkbox group
    output$chromosome_subset_checkboxgroup <- shiny::renderUI({
      chromosomes <- levels(as.factor(load_block_dataframe()$chr))
      grups.plots::shiny_reactive_checkbox(
        values  = chromosomes,
        title   = "Display chromosome(s):",
        input_id = "chromosome_labels",
        ncol    = 6
      )
    })

    # ---- 2c. Render block-scatterplot
    output$block_scatterplot <- plotly::renderPlotly({
      grups.plots::plot_sliding_window(
        load_block_dataframe(),
        input$block_pair,
        chromosomes = input$chromosome_labels
      )
    })

    # ---- 2d. Render block dataframe
    output$block_dataframe <- shiny::renderTable(
      load_block_dataframe()
    )

    # ---- 3a. Load / Update simulations dataframe
    load_sims_dataframe <- shiny::reactive(
      grups.plots::load_simfile(
        path = sim_files[input$sim_pair, ]
      )
    )

    # ---- 3b. Output simulation labels checkbox group.
    output$simulation_labels_checkboxgroup <- shiny::renderUI({
      sim_labels <- levels(load_sims_dataframe()$label)
      grups.plots::shiny_reactive_checkbox(
        values   = sim_labels,
        title    = "Display Relationship Label(s):",
        input_id = "violin_labels",
        ncol     = 2
      )
    })

    # ---- 3c. Render simulations violin plot
    output$sims_violinplot <- plotly::renderPlotly({
      grups.plots::plot_pedigree_sims(
        sims_dataframe = load_sims_dataframe(),
        pwd_path       = pwd_files[1],
        pair           = input$sim_pair,
        labels_to_plot = input$violin_labels
      )
    })

    # ---- 3d. Render simulations dataframe
    output$sims_dataframe <- shiny::renderTable(
      load_sims_dataframe()
    )

    # ---- 3e. Render KS normality test
    output$ks_normality_test <- DT::renderDataTable(
      grups.plots::test_normality(load_sims_dataframe(), alpha = input$ks_alpha) %>% 
      DT::datatable(style = "bootstrap", options = list(dom = 't', ordering=F, scrollX = T, width="100%")) %>% 
      formatStyle(1:999, rows="p.val", color = JS(paste("value  < ", input$ks_alpha," ? 'red' : ''"))),
      rownames = TRUE
    )

    # ---- 3f. Render BC matrix
    output$bc_matrix <- shiny::renderTable(
      grups.plots::get_bc_matrix(load_sims_dataframe()),
      rownames = TRUE
    )

  }

  shiny::shinyApp(ui, server, ...)
}