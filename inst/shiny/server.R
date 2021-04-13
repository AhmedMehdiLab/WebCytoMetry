library(magrittr)
library(shiny)
library(WebCytoMetry)

NONE <- c("None" = "__NULL__")
CLUSTER <- c("Cluster" = "cluster_id")
CONDITION <- c("Condition" = "condition")
PATIENT <- c("Patient" = "patient_id")
SAMPLE <- c("Sample" = "sample_id")

check <- function(x) if (x != 0 && x != NONE) x

server <- function(input, output, session) {
  return()
  # data storage
  cache <- reactiveValues()
  cache$items <- list()

  # populate cache
  isolate({
    item <- import_fcs_root(system.file("extdata", package = "WebCytoMetry"))
    for (name in names(item)) cache$items[[name]] <- item[[name]]
  })

  # load essential data
  item <- reactive({req(input$home_source); cache$items[[input$home_source]]})
  flow <- reactive({if (is.null(input$home_trans)) item()$data else flowCore::transform(item()$data, flowCore::transformList(input$home_trans, flowCore::arcsinhTransform(a = input$home_trans_a, b = input$home_trans_b, c = input$home_trans_c)))})
  meta <- reactive({req(flow()); flowCore::pData(flow())})
  panel <- reactive({item()$panel})

  ch_antigen <- reactive({setNames(panel()$antigen, panel()$external)})
  ch_colname <- reactive({setNames(panel()$fcs_colname, panel()$external)})

  # extract expression data
  flow_expr <- reactive({extract_expressions(flow(), panel())})
  flow_exps <- reactive({scale_expressions(flow_expr(), input$xp_range)})
  exp_label <- reactive({extract_labels(flow())})

  # update primary inputs
  observe({updateSelectInput(session, "home_source", choices = names(cache$items))})
  observe({updateSelectInput(session, "home_trans", choices = ch_colname())})

  observe({updateSelectInput(session, "som_channels", choices = ch_antigen())})
  observe({updateSelectInput(session, "som_meta", choices = stringr::str_c("meta", seq_len(input$som_max)))})

  observe({updateSelectInput(session, "heat_expr", choices = ch_antigen())})
  observe({updateSelectInput(session, "heat_freq", choices = c("Abundances" = "abundances", ch_antigen()))})

  observe({updateSelectInput(session, "dr_color", choices = c(CLUSTER, CONDITION, SAMPLE, PATIENT, ch_antigen()))})
  observe({updateSelectInput(session, "dr_facet", choices = c(NONE, CLUSTER, CONDITION, SAMPLE, PATIENT))})

  observe({updateSelectInput(session, "mem_channels", choices = ch_antigen())})
  observe({updateSelectInput(session, "mem_background", choices = ch_antigen())})
  observe({updateSelectInput(session, "sim_channels", choices = ch_colname())})
  observe({updateSelectInput(session, "ch_channels", choices = ch_antigen())})

  # primary calculations
  cluster <- reactive({
    validate(need(length(input$som_channels) > 1, "Select two or more channels for clustering"))
    do_cluster(flow(), panel(), input$som_channels, input$som_x, input$som_y, input$cl_clusters)
  })
  som <- reactive({
    validate(need(length(input$som_channels) > 1, "Select two or more channels for clustering"))
    do_som(flow(), input$som_channels, input$som_x, input$som_y, input$cl_clusters)
  })

  # extract annotated expression data
  expr_anno <- reactive({collate_expressions(flow_expr(), exp_label(), cluster(), input$files_min)})
  exps_anno <- reactive({collate_expressions(flow_exps(), exp_label(), cluster(), input$files_min)})

  # update secondary inputs
  observe({updateSelectInput(session, "ch_colour", choices = names(expr_anno()))})
  observe({updateSelectInput(session, "ch_facet", choices = names(expr_anno()))})

  # secondary calculations
  mem <- reactive({calc_mem(expr_anno())})
  sim_cosine <- reactive({calc_sim_cosine(exps_anno())})
  sim_hilbert <- reactive({
    validate(need(!is.null(input$sim_channels), "Select one or more channels"))
    calc_sim_hilbert(exps_anno(), input$sim_channels, input$sim_bins, input$sim_dims)
  })

  radviz_plot <- reactive({plot_radviz(exps_anno(), input$ch_channels)})
  glmm <- reactive({compute_glmm(expr_anno(), input$som_meta, input$glmm_min) %>% collate_results()})
  lmer <- reactive({compute_lmer(expr_anno(), input$som_meta, input$lmer_min) %>% collate_results()})

  # generate outputs
  observeEvent(input$home_source, {
    new_index <- seq_along(flow())
    new_data <- list()

    for (i in cache$home_view) removeTab("home_view", i)
    lapply(rev(new_index), function(i) {
      new_id <- stringr::str_c("view_", i)
      new_data[[new_id]] <- flow()[[i]]
      prependTab("home_view", tabPanel(i, verbatimTextOutput(new_id)), TRUE)
      output[[new_id]] <- renderPrint(flow()[[i]])
    })

    cache$home_view <- new_index
  })

  # display outputs
  output$home_main <- renderPrint({flow()})

  output$plot_abnd <- renderPlot({CATALYST::plotAbundances(cluster(), input$som_meta, input$abnd_mode, check(input$abnd_group), check(input$abnd_shape))})
  output$plot_code <- renderPlot({CATALYST::plotCodes(cluster(), input$som_meta)})
  output$plot_clxp <- renderPlot({CATALYST::plotClusterExprs(cluster(), input$som_meta, NULL)})
  output$plot_heat <- renderPlot({
    validate(need(!is.null(input$heat_expr), "Select one or more channels"))
    CATALYST::plotMultiHeatmap(cluster(), input$heat_expr, input$heat_freq, input$som_meta, bars = TRUE, perc = TRUE)
  })
  output$plot_dr <- renderPlot({
    cluster() %>%
      CATALYST::runDR(input$dr_method, check(input$dr_cells), features = input$som_channels) %>%
      CATALYST::plotDR(NULL, input$dr_color, check(input$dr_facet))
  })
  output$plot_star <- renderPlot({FlowSOM::PlotStars(som()$FlowSOM, backgroundValues = som()$metaclustering)})

  output$xp_matrix <- renderDataTable({flow_expr()})
  output$xp_scaled <- renderDataTable({flow_exps()})
  output$file_view <- renderText({stringr::str_c("Selected:\n", stringr::str_c(files_inc()$files, collapse = "\n"))})

  output$mem_view <- renderDataTable({mem()})
  output$mem_plot <- renderPlot({
    mem() %>%
      ggplot2::ggplot(ggplot2::aes(x = file_name, y = MEM, color = !!rlang::sym(input$mem_colour))) +
      ggplot2::geom_point() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
  })

  output$simh_view <- renderDataTable({sim_hilbert()$hilbert})
  output$simh_plot <- renderPlot({
    sim_hilbert()$hilbert %>%
      ggplot2::ggplot(ggplot2::aes(x = V1, y = V2, color = !!rlang::sym(input$sim_colour), size = Cells)) +
      ggplot2::geom_point()
  })
  output$simh_hier <- renderPlot({plot(hclust(sim_hilbert()$dists))})
  output$simc_view <- renderDataTable({sim_cosine()})
  output$simc_heat <- renderPlot({heatmap(sim_cosine(), symm = TRUE)})
  output$projection <- renderPlot({
    validate(need(length(input$ch_channels) > 1, "Select two or more channels"))
    plot_radviz(flow_exps(), metaclust(), sim_cosine(), input$ch_channels, input$ch_colour, input$ch_facet)
  })

  output$som_volcano <- renderPlot({calc_glmmadmb(metaclust(), files_inc(), input$som_meta) %>% ggplot2::ggplot(ggplot2::aes(Estimate, logPValue, color = Meta))})
  output$mem_volcano <- renderPlot({calc_contrast(mem(), input$mem_background)})

  output$debug <- renderText({browser(); TRUE})
}
