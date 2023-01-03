#' Main shiny dashboard application
#' @import shiny
#' @import plotly
#' @import stringr
#' @import shinycssloaders
#' @export
app <- function(ui, 
		server, 
		data_dir = "./grups_output", 
		sample_regex = "[A-Za-z0-9]+([-0-9]+){0,1}", 
		...
) {

  # ----- Format a file pair regular expression
  pair_regex = paste0("(?<=-)(",sample_regex,"-",sample_regex,")")

  # ---- 0a. Configure loading spinner animation
  options(spinner.type = 8, spinner.color = "#0dc5c1")

  # ---- 1. Search for .result file(s)
  res_files <- list.files(
    path = data_dir,
    full.names = TRUE,
    pattern = "\\.result$"
  )

  # ---- 2. Search for .pwd file(s)
  pwd_files <- list.files(
    path = data_dir,
    full.names = TRUE,
    pattern = "\\.pwd$"
  )

  # ---- 3a. Search for blk_file(s)                                   [A FUNC]
  blk_files <- list.files(
    path = paste(data_dir, "blocks", sep = "/"),
    full.names = TRUE,
    pattern = "\\.blk$"
  )

  # ---- 3b. Extract block pair names, parse all that data into a df. [B FUNC]
  blk_files <- data.frame(
    path      = blk_files,
    row.names = stringr::str_extract(  # Extract pair names
      blk_files,
      paste0(pair_regex, "(?=.blk$)")
    ),
    stringsAsFactors = FALSE
  )

  # ---- 4a. Search for simulation files                              [A FUNC]
  sim_files <- list.files(
    path = paste(data_dir, "simulations", sep = "/"),
    full.names = TRUE,
    pattern = "\\.sims$"
  )
  # ---- 4b. Extract simulations, pair names, parse them into a df. [B FUNC]
  sim_files <- data.frame(
    path = sim_files,
    row.names = stringr::str_extract( # Extract pair names
      sim_files,
      paste0(pair_regex, "(?=.sims$)")
    ),
    stringsAsFactors = FALSE
  )

  # ---- 5a. Search for .yaml config files
  config_files <- list.files(
    path = data_dir,
    full.names = TRUE,
    pattern = "\\.yaml$"
  )

  # ---- 5b. order them in decreasing order. -> last yaml becomes the first.
  config_files <- config_files[order(config_files, decreasing = TRUE)]


  # ---- Sanity checks
  shiny::validate(
    shiny::need(
      length(pwd_files) == 1,
      "[ERROR]: Exactly one `.pwd` file must exist within `data_dir`. \
       Exiting."
    ),

    shiny::need(
      length(res_files) == 1,
      "[ERROR]: Exactly one `.result` file must exist within `data_dir`. \
       Exiting."
    ),

    shiny::need(
      length(config_files) >= 1,
      "[ERROR]: At least one `.yaml` configuration file must exist within \
       `data_dir`. Exiting"
    )
  )

  ui <- shiny::fluidPage(
    theme = bslib::bs_theme(bootswatch = "darkly", version = 4),
    shiny::navbarPage("GRUPS-plots",
      # ---- 1. Render summary table.
      shiny::tabPanel("Summary",
        DT::dataTableOutput("simulation_results_df") %>%
          shinycssloaders::withSpinner()
      ),
      # ---- 1. Render pwd barplots tab.
      shiny::tabPanel("Raw Pairwise differences.",
        shiny::fluidPage(
          shiny::sidebarLayout(
            shiny::sidebarPanel(
              shiny::fluidRow(
                shiny::column(6,
                  shiny::checkboxInput("norm_request",
                    label = "Normalize"
                  )
                ),
                shiny::column(6,
                  shiny::checkboxInput("hide_self_comparisons",
                    label = "Hide self-comparisons"
                  )
                )
              ),

              shiny::fluidRow(
                shiny::column(6,
                  shiny::radioButtons("norm_method",
                    "Normalization subset",
                    choiceNames  = c("Pairwise-comparisons (2Ms)",
                                     "Self-comparisons (Ms)",
                                     "Value (2Ms)"
                                    ),
                    choiceValues = c("Pairwise", "Self", "Value"),
                    selected = "Pairwise",
                    width = "100%"
                  )
                ),
                shiny::column(6,
                  shiny::conditionalPanel(
                    condition = "input.norm_method != 'Value'",
                    shiny::radioButtons("norm_metric",
                      label    = "Normalization metric",
                      choices  = c("Median", "Mean", "Min", "Max"),
                      selected = "Median"
                    )
                  ),
                  shiny::conditionalPanel(
                    condition = "input.norm_method == 'Value'",
                    shiny::numericInput("norm_value",
                      label = "Normalization value",
                      value = 0.25,
                      min   = 0.01,
                      max   = 2,
                      step  = 0.01
                    )
                  )
                )
              ),
              shiny::fluidRow(
                shiny::column(6,
                  shiny::numericInput("min_overlap",
                    label = "Minimum SNP overlap treshold",
                    value = 0,
                    min   = 0,
                    max   = .Machine$integer.max,
                    step  = 10000
                  ),
                ),
                shiny::column(6,
                  shiny::radioButtons("norm_avg_type",
                    label = "Type",
                    choiceNames = c("Raw Average PWD", "Corrected Average PWD"),
                    choiceValues = c("Raw", "Corr"),
                  )
                )
              )
            ),
            shiny::mainPanel(
              shiny::tabsetPanel(
                shiny::tabPanel("Plot",
                  plotly::plotlyOutput("pwd_barplot") %>%
                    shinycssloaders::withSpinner()
                ),
                shiny::tabPanel("Raw data",
                   DT::dataTableOutput("pwd_dataframe") %>%
                    shinycssloaders::withSpinner()
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
                max   = 50
              ),
              shiny::sliderInput("block_step",
                "Sliding window step",
                value = 1,
                min = 1,
                max = 50
              ),
              shiny::uiOutput("chromosome_subset_checkboxgroup") %>%
                shinycssloaders::withSpinner()
            ),
            shiny::mainPanel(
              shiny::tabsetPanel(
                shiny::tabPanel("Plot",
                  plotly::plotlyOutput("block_scatterplot") %>%
                    shinycssloaders::withSpinner()
                ),
                shiny::tabPanel("Raw data",
                  shiny::tableOutput("block_dataframe") %>%
                    shinycssloaders::withSpinner()
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
              shiny::uiOutput("simulation_labels_checkboxgroup") %>%
                shinycssloaders::withSpinner(),
              shiny::numericInput("ks_alpha",
                "Kolmogorov-Smirnov alpha treshold",
                min   = 0,
                max   = 1,
                step  = 0.01,
                value = 0.05
              )
            ),
            shiny::mainPanel(
              shiny::tabsetPanel(
                shiny::tabPanel("Plot",
                  plotly::plotlyOutput("sims_violinplot") %>%
                    shinycssloaders::withSpinner(),
                  DT::dataTableOutput("ks_normality_test") %>%
                    shinycssloaders::withSpinner(),
                  shiny::hr(),
                  shiny::fluidRow(
                    shiny::column(6,
                      plotly::plotlyOutput("bc_matrix_plot") %>%
                        shinycssloaders::withSpinner()
                    ),
                    shiny::column(6,
                      plotly::plotlyOutput("plot_or_matrix") %>%
                        shinycssloaders::withSpinner()
                    )
                  ),
                  shiny::hr(),
                  plotly::plotlyOutput("OR_confidence")
                ),
                shiny::tabPanel("Raw data",
                  shiny::tableOutput("sims_dataframe") %>%
                    shinycssloaders::withSpinner()
                )
              )
            )
          )
        )
      ),

      # ---- 4. Render yaml configuration file:
      shiny::tabPanel("Configuration",
        shiny::verbatimTextOutput("config_file") %>%
          shinycssloaders::withSpinner()
      )
    )
  )
  server <- function(input, output, session) {

  # 0 ---- Update block slider inputs
  shiny::observe({
    max_blockstep_value <- input$block_width - 1
    new_step_value <- ifelse(
      input$block_step < max_blockstep_value,
      input$block_step,
      max_blockstep_value
    )
    shiny::updateSliderInput(
      session,
      "block_step",
      value = new_step_value,
      min = 1,
      max = max_blockstep_value
    )
  })


    load_results_file <- shiny::reactive({
      grups.plots::load_res_file(res_files[1])
    })

    # ---- 1a. Load results  summary dataframe
    output$simulation_results_df <- DT::renderDataTable({
      DT::datatable(
        grups.plots::merge_pwd_results(
          pwd_df = load_pairwise_dataframe(),
          res_df = load_results_file()
        ),
        style = "bootstrap5",
        options = list(ordering = TRUE, scrollX = FALSE),
        class = "table-condensed",
        filter = "top",
      )
  })

    # ---- 1b. Load / Update pairwise dataframe.
    load_pairwise_dataframe <- shiny::reactive(
      grups.plots::load_pairwise_file(
        path          = pwd_files[1],
        res_data      = load_results_file(),
        norm_avg_type = input$norm_avg_type,
        sample_regex  = sample_regex,
        min_overlap   = input$min_overlap,
        norm_request  = input$norm_request,
        norm_method   = input$norm_method,
        norm_metric   = input$norm_metric,
        norm_value    = input$norm_value
      )
    )

    # ---- 2a. Render pairwise difference barplot.
    output$pwd_barplot <- plotly::renderPlotly(
      grups.plots::plot_pairwise_diff(
        data                  = load_pairwise_dataframe(),
        hide_self_comparisons = input$hide_self_comparisons,
        norm_method           = input$norm_method,
        norm_metric           = input$norm_metric,
      )
    )

    # ---- 2b. Render pairwise dataframe
    output$pwd_dataframe <- DT::renderDataTable(
      DT::datatable(
        load_pairwise_dataframe()$data,
        style = "bootstrap5",
        options = list(ordering = TRUE, scrollX = TRUE),
        class = "table-condensed",
        filter = "top",
      )

    )

    # ---- 3a. Load / Update block dataframe
    load_block_dataframe <- shiny::reactive(
      grups.plots::load_blockfile(
        path  = blk_files[input$block_pair, ],
        width = input$block_width,
        step  = input$block_step
      )
    )

    # ---- 3b. Output chromosome subset checkbox group
    output$chromosome_subset_checkboxgroup <- shiny::renderUI({
      chromosomes <- levels(as.factor(load_block_dataframe()$chr))
      grups.plots::shiny_reactive_checkbox(
        values  = chromosomes,
        title   = "Display chromosome(s):",
        input_id = "chromosome_labels",
        ncol    = 6
      )
    })

    # ---- 3c. Render block-scatterplot
    output$block_scatterplot <- plotly::renderPlotly({
      grups.plots::plot_sliding_window(
        load_block_dataframe(),
        pair            = input$block_pair
      )
    })

    ## WIP: Updates traces without loading up the entire plotly widget.
    ## Ugly, finicky, and most likely error prone code at the moment.
    #shiny::observeEvent(input$block_width, {
    #  plotly::plotlyProxy("block_scatterplot", session) %>%
    #  plotly::plotlyProxyInvoke("addTraces",
    #    lapply(1:22, FUN = function(x) {
    #      data=load_block_dataframe()
    #      list(x = data[which(data$chr == x),]$start,
    #           y = data[which(data$chr == x),]$avg_pwd,
    #           color = ~as.factor(x),
    #           #colors = RColorBrewer::brewer.pal(n=22, "Set2"),
    #           name = x,
    #           group = as.factor(x),
    #           type = 'scatter',
    #           mode = 'lines+markers'
    #      )
    #    })
    #  ) %>%
    #  plotly::plotlyProxyInvoke("deleteTraces", as.integer(0:21))
    #})
    # ---- 3d. Filter out chromosome which were not requested by the user.
    shiny::observeEvent(input$chromosome_labels, {
      chr_to_hide <- unique(
        load_block_dataframe()$chr[
          !(load_block_dataframe()$chr %in% input$chromosome_labels)
        ] - 1
      )
      chr_to_keep <- unique(
        load_block_dataframe()$chr[
          (load_block_dataframe()$chr %in% input$chromosome_labels)
        ] - 1
      )
      plotly::plotlyProxy("block_scatterplot",
        session
      ) %>%
      plotly::plotlyProxyInvoke("restyle",
        list(visible = FALSE),
        chr_to_hide
      ) %>%
      plotly::plotlyProxyInvoke("restyle",
        list(visible = TRUE),
        chr_to_keep
      )

    })

    # ---- 3e. Reset user-selected chromosome if he used select/deselect All.
    shiny::observeEvent(input$chromosome_labels_select, {
      shiny::updateCheckboxGroupInput(
        session  = session,
        inputId  = "chromosome_labels",
        choices  = levels(as.factor(load_block_dataframe()$chr)),
        selected = levels(as.factor(load_block_dataframe()$chr))
      )
    })
    shiny::observeEvent(input$chromosome_labels_deselect, {
      shiny::updateCheckboxGroupInput(
        session  = session,
        inputId  = "chromosome_labels",
        choices  = levels(as.factor(load_block_dataframe()$chr)),
        selected = c()
      )
    })

    # ---- 3f. Render block dataframe
    output$block_dataframe <- shiny::renderTable(
      load_block_dataframe()
    )

    # ---- 4a. Load / Update simulations dataframe
    load_sims_dataframe <- shiny::reactive(
      grups.plots::load_simfile(
        path = sim_files[input$sim_pair, ]
      )
    )

    # ---- 4b. Output simulation labels checkbox group.
    output$simulation_labels_checkboxgroup <- shiny::renderUI({
      sim_labels <- levels(load_sims_dataframe()$label)
      grups.plots::shiny_reactive_checkbox(
        values   = sim_labels,
        title    = "Display Relationship Label(s):",
        input_id = "violin_labels",
        ncol     = 2
      )
    })

    # ---- 4c. Render simulations violin plot
    output$sims_violinplot <- plotly::renderPlotly({
      grups.plots::plot_pedigree_sims(
        sims_dataframe = load_sims_dataframe(),
        results_data   = load_results_file(),
        pair           = input$sim_pair,
        labels_to_plot = input$violin_labels
      )
    })

    # ---- 4f. Reset User-selected chromosome if he used select/deselect All.
    shiny::observeEvent(input$violin_labels_select, {
      shiny::updateCheckboxGroupInput(
        session  = session,
        inputId  = "violin_labels",
        choices = levels(load_sims_dataframe()$label),
        selected = levels(load_sims_dataframe()$label)
      )
    })
    shiny::observeEvent(input$violin_labels_deselect, {
      shiny::updateCheckboxGroupInput(
        session  = session,
        inputId  = "violin_labels",
        choices  = levels(load_sims_dataframe()$label),
        selected = c()
      )
    })

    # ---- 4g. Render simulations dataframe
    output$sims_dataframe <- shiny::renderTable(
      load_sims_dataframe()
    )

    # ---- 4h. Render KS normality test
    output$ks_normality_test <- DT::renderDataTable(
      grups.plots::test_normality(
        sims_file = load_sims_dataframe(),
        alpha     = input$ks_alpha
      ) %>%
      DT::datatable(
        style   = "bootstrap5",
        options = list(dom      = "t",
                       ordering = FALSE,
                       scrollX  = TRUE,
                       width    = "100%"
                      )
      ) %>%
      DT::formatStyle(1:999,
                      rows  = "p.val",
                      color = DT::JS(paste("value  < ",
                                           input$ks_alpha,
                                           " ? 'red' : ''"
                                          )
                                    )
      ),
      rownames = TRUE
    )

    # ---- 4i. Render BC matrix plot
    output$bc_matrix_plot <- plotly::renderPlotly(
      grups.plots::plot_bc_matrix(
        bc_matrix = grups.plots::get_bc_matrix(
          sims_data = load_sims_dataframe(),
          labels_to_keep = input$violin_labels
        ),
        plot_title = "<b>Matrix of Bhattacharya coefficients</b>",
        cutoff_values = c(-Inf, 0.01, 0.05, Inf),
        cutoff_labels = c("< 0.01", "< 0.05", ">= 0.05"),
        absolute_values = FALSE,
        marker_text = "<b>BC coefficient:</b>"
      )
    )

    # ---- 4j. Render  OR matrix
    output$plot_or_matrix <- plotly::renderPlotly({
      fig <- grups.plots::plot_bc_matrix(
        bc_matrix = grups.plots::get_odds_matrix(
                      sims_data        = load_sims_dataframe(),
                      observed_results = load_results_file(),
                      pair             = input$sim_pair,
                      labels_to_keep   = input$violin_labels
                    ),
        plot_title = "<b>Matrix of log(Odds Ratio)</b>",
        marker_text = NULL,
        cutoff_values = c(-Inf, log(1), log(100), Inf),
        cutoff_labels = c("<= 1/1", "<= 100/1 ", "> 100/1"),
        cutoff_colors = c("#CC6677", "#DDCC77", "#44AA99"),
        right_align = TRUE,
        absolute_values = TRUE
      )

      # Dirty trick to obtain a proper marker text annotation.
      fig[["x"]][["attrs"]][[1]][["text"]] <- ~paste(
        "<b>Odds Ratio:</b>", exp(value), "<br>",
        "<b>Log(OR):</b>", value
      )
      fig
    })

    # ---- 4k. Render OR confidence
    output$OR_confidence <- plotly::renderPlotly(
      grups.plots::plot_or_confidence(
        or_matrix = grups.plots::get_odds_matrix(
                      sims_data        = load_sims_dataframe(),
                      observed_results = load_results_file(),
                      pair             = input$sim_pair,
                      labels_to_keep = input$violin_labels
        ),
        predictor = load_results_file()[
          which(load_results_file()[, 1] == input$sim_pair), 2
        ]
      )
    )

    # ---- 5a. Render yaml config file
    output$config_file <- shiny::renderText(
      readChar(config_files[1], file.info(config_files[1])$size)
    )
  }

  options(shiny.autoreload = TRUE)

  shiny::shinyApp(ui, server, ...)
}
