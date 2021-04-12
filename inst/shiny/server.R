library(magrittr)
library(shiny)
library(WebCytoMetry)

server <- function(input, output, session) {

  # data storage
  cache <- reactiveValues()
  cache$items <- list()

  # populate cache
  isolate({
    #item <- import_fcs_root(system.file("extdata", package = "WebCytoMetry"))
    #for (name in names(item)) cache$items[[name]] <- item[[name]]
    cache$items$KatJanin <- readRDS(system.file("extdata", "KatJanin.rds", package = "WebCytoMetry"))
  })

  # extract data
  item <- reactive({
    req(input$home_source)
    cache$items[[input$home_source]]
  })
  flow <- reactive({
    if (is.null(input$home_trans)) item()$data
    else flowCore::transform(item()$data, flowCore::transformList(input$home_trans, flowCore::arcsinhTransform(a = input$home_trans_a, b = input$home_trans_b, c = input$home_trans_c)))
  })
  meta <- reactive({
    req(flow())
    flowCore::pData(flow())
  })
  panel <- reactive({item()$panel})
  ch_antigen <- reactive({setNames(panel()$antigen, panel()$external)})
  ch_colname <- reactive({setNames(panel()$fcs_colname, panel()$external)})

  meta_labs <- reactive({extract_labels(flow())})
  metaclust <- reactive({
    validate(need(isTruthy(cluster()), "Cluster data before proceeding"))

    labels <- meta_labs()
    labels$cluster <- cluster()$cluster_id
    labels %>% dplyr::left_join(CATALYST::cluster_codes(cluster()), by = c("cluster" = "som100"))
  })
  code_name <- reactive({req(cluster()); cluster() %>% CATALYST::cluster_codes() %>% names() %>% stringr::str_subset("meta")})
  flow_expr <- reactive({extract_expressions(flow(), panel())})
  flow_exps <- reactive({scale_expressions(flow_expr(), input$xp_range)})
  files_inc <- reactive({
    req(input$files_min)
    extract_files(meta_labs(), input$files_min)
  })

  # update inputs
  observe({updateSelectInput(session, "home_source", choices = names(cache$items))})
  observe({updateSelectInput(session, "home_trans", choices = ch_colname())})
  observe({updateSelectInput(session, "cl_channels", choices = ch_antigen())})
  observe({req(cluster()); updateSelectInput(session, "cl_k", choices = code_name())})
  observe({updateSelectInput(session, "heat_channels", choices = ch_antigen())})
  observe({updateSelectInput(session, "heat_axis", choices = c("Abundances" = "abundances", "State" = "state", ch_antigen()))})
  observe({updateSelectInput(session, "dr_color", choices = c(intersect(names(meta()), c("condition", "sample_id", "patient_id")), code_name(), ch_antigen()))})
  observe({updateSelectInput(session, "dr_facet", choices = c("None" = "__none__", intersect(names(meta()), c("condition", "sample_id", "patient_id"))))})
  observe({updateSelectInput(session, "mem_channels", choices = ch_antigen())})
  observe({updateSelectInput(session, "mem_background", choices = ch_antigen())})
  observe({updateSelectInput(session, "sim_channels", choices = ch_antigen())})
  observe({updateSelectInput(session, "ch_channels", choices = ch_antigen())})
  observe({updateSelectInput(session, "ch_colour", choices = names(metaclust()))})
  observe({updateSelectInput(session, "ch_facet", choices = names(metaclust()))})

  # calculate
  cluster <- reactive({
    validate(need(length(input$cl_channels) > 1, "Select two or more channels"))
    calc_cluster(flow(), panel(), input$cl_channels, input$cl_clusters)
  })

  mem <- reactive({calc_mem(flow_expr(), meta_labs(), input$mem_channels, input$mem_sample)})
  som <- reactive({calc_som(flow(), panel()$fcs_colname[which(panel()$antigen %in% input$cl_channels)], input$cl_clusters)})

  som_mapping <- reactive({som()$FlowSOM$map$mapping[, 1]})
  sim_hilbert <- reactive({
    validate(need(!is.null(input$sim_channels), "Select one or more channels"))
    calc_sim_hilbert(flow_exps(), meta_labs(), files_inc()$files, input$sim_channels, input$sim_bins, input$sim_dims)
  })
  sim_cosine <- reactive({calc_sim_cosine(flow_exps(), meta_labs(), files_inc()$files)})

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

  output$cl_abundance <- renderPlot({CATALYST::plotAbundances(cluster(), input$cl_k, input$cl_abundance_by)})
  output$cl_code <- renderPlot({CATALYST::plotCodes(cluster(), input$cl_k)})
  output$cl_expr <- renderPlot({CATALYST::plotClusterExprs(cluster(), input$cl_k, NULL)})
  output$expr_heat <- renderPlot({
    validate(
      need(!is.null(input$heat_channels), "Select one or more channels"),
      need(input$heat_by != "both" || length(input$heat_channels) == 1, "Select one channel to aggregate by both Cluster ID and Sample ID")
    )
    CATALYST::plotExprHeatmap(cluster(), input$heat_channels, input$heat_by, input$cl_k, bars = TRUE, perc = TRUE)
  })
  output$freq_heat <- renderPlot({CATALYST::plotFreqHeatmap(cluster(), input$cl_k)})
  output$dr_plot <- renderPlot({
    cluster() %>%
      CATALYST::runDR(input$dr_method, if (input$dr_cells != 0) input$dr_cells, features = NULL) %>%
      CATALYST::plotDR(NULL, input$dr_color, if (input$dr_facet != "__none__") input$dr_facet)
  })

  output$xp_matrix <- renderDataTable({flow_expr()})
  output$xp_scaled <- renderDataTable({flow_exps()})
  output$file_view <- renderText({stringr::str_c("Selected:\n", stringr::str_c(files_inc()$files, collapse = "\n"))})
  output$file_plot <- renderPlot({files_inc()$plot})

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

  output$som_star <- renderPlot({FlowSOM::PlotStars(som()$FlowSOM, backgroundValues = som()$metaclustering)})
  output$som_volcano <- renderPlot({calc_glmmadmb(metaclust(), files_inc(), input$cl_k) %>% ggplot2::ggplot(ggplot2::aes(Estimate, logPValue, color = Meta))})
  output$mem_volcano <- renderPlot({calc_contrast(mem(), input$mem_background)})

  output$debug <- renderText({browser(); TRUE})
}
