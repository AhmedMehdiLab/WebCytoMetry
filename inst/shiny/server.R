library(magrittr)
library(shiny)
library(WebCytoMetry)

server <- function(input, output, session) {

  # store values
  cache <- reactiveValues()
  flow_data <- reactiveValues()
  flow <- reactive({req(input$home_source); flow_data[[input$home_source]]})

  # load data
  isolate({
    item <- import_fcs_root(system.file("extdata", package = "WebCytoMetry"))
    for (name in names(item)) flow_data[[name]] <- item[[name]]
  })

  # extract data
  data <- reactive({if (is.null(input$home_trans)) flow()$data else flowCore::transform(flow()$data, flowCore::transformList(input$home_trans, flowCore::arcsinhTransform(a = input$home_trans_a, b = input$home_trans_b, c = input$home_trans_c)))})
  meta <- reactive({req(data()); flowCore::pData(data())})
  panel <- reactive({flow()$panel})
  channels <- reactive({setNames(panel()$antigen, panel()$external)})

  # update inputs
  observe({updateSelectInput(session, "home_source", choices = names(flow_data))})
  observe({updateSelectInput(session, "home_trans", choices = setNames(panel()$fcs_colname, panel()$external))})
  observe({updateSelectInput(session, "cl_channels", choices = channels())})
  observe({updateSelectInput(session, "heat_channels", choices = channels())})
  observe({updateSelectInput(session, "heat_axis", choices = c("Abundances" = "abundances", "State" = "state", channels()))})
  observe({updateSelectInput(session, "dr_color", choices = c("meta20", intersect(names(meta()), c("condition", "sample_id", "patient_id")), channels()))})
  observe({updateSelectInput(session, "dr_facet", choices = c("None" = "__none__", intersect(names(meta()), c("condition", "sample_id", "patient_id"))))})

  # calculate
  cluster <- reactive({validate(need(length(input$cl_channels) > 1, message = "Select two or more channels")); data() %>% CATALYST::prepData(panel(), meta()) %>% CATALYST::cluster(input$cl_channels)})
  reduction <- reactive({CATALYST::runDR(cluster(), input$dr_method, if (input$dr_cells != 0) input$dr_cells, features = NULL)})

  expression <- reactive({get_expressions(data(), panel())})
  exp_scaled <- reactive({scale_expressions(expression()$matrix, input$xp_range)})
  files_incl <- reactive({req(input$filt_cell); filter_files(expression()$labels, input$filt_cell)})

  gate <- reactive({flowWorkspace::GatingSet(data())})

  # generate outputs
  observeEvent(input$home_source, {
    new_index <- seq_along(data())
    new_data <- list()

    for (i in cache$home_view) removeTab("home_view", i)
    lapply(new_index, function(i) {
      new_id <- stringr::str_c("view_", i)
      new_data[[new_id]] <- data()[[i]]
      new_panel <- tabPanel(i, verbatimTextOutput(new_id))
      output[[new_id]] <- renderPrint(data()[[i]])
      appendTab("home_view", new_panel)
    })

    cache$home_view <- new_index
  })

  # display outputs
  output$home_info <- renderPrint({data()})
  output$cl_abundance <- renderPlot({CATALYST::plotAbundances(cluster(), by = input$cl_abundance)})
  output$cl_code <- renderPlot(CATALYST::plotCodes(cluster()))
  output$cl_expr <- renderPlot({CATALYST::plotClusterExprs(cluster(), features = NULL)})
  output$heat_expr <- renderPlot({
    validate(
      need(length(input$heat_channels) != 0, "Select one or more channels to display"),
      need(input$heat_by != "both" || length(input$heat_channels) == 1, "Select one channel to aggregate by both Cluster ID and Sample ID")
    )
    CATALYST::plotExprHeatmap(cluster(), input$heat_channels, input$heat_by, bars = TRUE, perc = TRUE)
  })
  output$heat_freq <- renderPlot({CATALYST::plotFreqHeatmap(cluster())})
  output$xp_matrix <- renderDataTable({expression()$matrix})
  output$xp_scaled <- renderDataTable({exp_scaled()})
  output$filt_view <- renderPlot({files_incl()$plot})
  output$filt_file <- renderText({stringr::str_c("Selected:\n", stringr::str_c(files_incl()$files, collapse = "\n"))})
  output$dr_view <- renderPlot({withProgress(CATALYST::plotDR(reduction(), NULL, input$dr_color, if (input$dr_facet != "__none__") input$dr_facet), message = "Loading...")})

  output$debug <- renderText({browser(); TRUE})

  # load input data
  # raw_flow <- reactive({raw_flows[[input$source]]})
  # raw_data <- reactive({raw_flow()$data})
  # raw_pane <- reactive({raw_flow()$pane})
  # raw_meta <- reactive({flowCore::pData(raw_data())})
  #
  # # compute matrix
  # channels <- reactive({
  #   cols <- setNames(flowCore::colnames(raw_data()), raw_pane()$antigen)
  #   cols <- cols[cols %in% input$channels]
  #   names(cols)})
  # phenofunc <- reactive({
  #   allc <- union(input$phenotypic, input$functional)
  #   cols <- setNames(flowCore::colnames(raw_data()), raw_pane()$antigen)
  #   cols <- cols[cols %in% allc]
  #   names(cols)})
  # m_data <- reactive({flow_matrix(gate(), channels = channels())})
  # matrix <- reactive({m_data()$matrix})
  # scaled <- reactive({m_data()$scaled})
  #
  # # annotate channels
  # annos <- reactive({flow_annos(gate())})
  # files <- reactive({flow_filter(annos(), input$filter)})
  # mem <- reactive({flow_calc_mem(matrix(), files(), annos(), channels())})
  # sim <- reactive({flow_similarity(annos(), scaled(), files(), phenofunc())})
  # chan <- reactive({flow_order_channels(scaled(), sim()$sim, phenofunc())})
  # som <- reactive({flow_start_flowSOM(gate(), input$functional, input$phenotypic)})
  # proj <- reactive({flow_proj_plot(annos(), scaled(), chan()$S, clusts = som()$clust, metas = som()$meta)})
  # sx <- reactive({flow_FlowSOM(matrix(), files(), som()$som, annos(), phenofunc(), som()$meta, som()$clust)})

  # output$clust_exp <- renderPlot(CATALYST::plotClusterExprs(cluster(), features = input$features))
  #
  # output$strains <- renderPlot(flow_plot_strain(annos()))
  # output$mem <- renderPlot(flow_plot_mem(mem()$mem))
  # output$sim <- renderPlot(sim()$plot)
  # output$simc <- renderPlot({plot(hclust(as.dist(as.matrix(sim()$dist)[files(), files()])))})
  # output$simp <- renderPlot({heatmap(sim()$sim, symm = TRUE)})
  # output$chano <- renderPlot(chan()$plot)
  # output$chand <- renderPlot(proj()$density)
  # output$chans <- renderPlot(proj()$strain)
  # output$chanm <- renderPlot({proj()$meta})
  # output$proje <- renderPlot(flow_proj_explore(mem()$df, chan()$S, files()))
  # output$ps <- renderPlot(sx()$stars)
  # output$px <- renderPlot(sx()$meta)
  # output$pv <- renderPlot(sx()$volcano)
  # output$pc <- renderPlot(sx()$contrast)
}
