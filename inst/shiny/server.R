library(magrittr)
library(shiny)
library(shinyjs)
library(WebCytoMetry)

NONE <- c("None" = "__NULL__")
CLUSTER <- c("Cluster" = "cluster_id")
CONDITION <- c("Condition" = "condition")
PATIENT <- c("Patient" = "patient_id")
SAMPLE <- c("Sample" = "sample_id")
COMMON_OPTS <- c("Channel" = "channel", "File" = "file", "Strain" = "strain", "Time" = "time")

check <- function(x) if (x != 0 && x != NONE) x
check_facet <- function(x) if (x != NONE) ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(x)))

server <- function(input, output, session) {
  shinyjs::disable("glmm_form")
  shinyjs::disable("lmer_form")
  shinyjs::disable("hypothesis")

  # data storage
  cache <- reactiveValues()
  cache$items <- list()
  cache$sce <- list()

  # populate cache
  isolate({
    item <- import_fcs_root(system.file("extdata", "fcs", package = "WebCytoMetry"))
    for (name in names(item)) cache$items[[name]] <- item[[name]]

    sces <- list.files(system.file("extdata", "sce", package = "WebCytoMetry"), full.names = TRUE)
    for (name in sces) cache$sce[[tools::file_path_sans_ext(basename(name))]] <- readRDS(name)
    cache$items$`Use SCE data` <- ""
  })

  # load essential data
  flow_not_sce <- reactive({input$home_source != "Use SCE data"})

  flow_item <- reactive({req(input$home_source); if (flow_not_sce()) cache$items[[input$home_source]]})
  flow_proc <- reactive(flow_item()) # {rescale_item(flow_item(), db_home_trans(), input$trans_a, input$trans_b, input$trans_c)})

  observeEvent(input$home_source, {shinyjs::toggleState("home_sce", !flow_not_sce())})
  flow_sce <- reactive({cache$sce[[input$home_sce]]})
  flow_meta <- reactive({if (flow_not_sce()) flow_item()$panel else tibble::tibble(fcs_colname = names(flow_sce()), antigen = make.names(names(flow_sce())), external = names(flow_sce()))})

  ch_antigen <- reactive({setNames(flow_meta()$antigen, flow_meta()$external)})
  ch_colname <- reactive({setNames(flow_meta()$fcs_colname, flow_meta()$external)})
  metaclusts <- reactive({stringr::str_c("meta", seq(2, input$som_max))})

  # load data
  observeEvent(input$home_upload, cache$items[[input$home_upload$name]] <- readRDS(input$home_upload$datapath))
  observeEvent(input$home_up_sce, cache$sce[[basename(input$home_up_sce$name)]] <- readRDS(input$home_up_sce$datapath))

  # update primary inputs
  observe({updateSelectInput(session, "home_source", choices = names(cache$items))})
  observe({updateSelectInput(session, "home_sce", choices = names(cache$sce))})
  observe({updateSelectInput(session, "home_trans", choices = ch_colname())})
  observe({updateSelectInput(session, "som_channels", choices = if (flow_not_sce()) ch_antigen() else ch_colname())})
  observe({updateSelectInput(session, "som_meta", choices = metaclusts())})
  observe({updateSelectInput(session, "mlth_freq", choices = c("Abundances" = "abundances", ch_antigen()))})
  observe({updateSelectInput(session, "dr_color", choices = c(CLUSTER, CONDITION, SAMPLE, PATIENT, if (flow_not_sce()) ch_antigen() else ch_colname(), metaclusts()), selected = CLUSTER)})
  observe({updateSelectInput(session, "dr_facet", choices = c(NONE, CLUSTER, CONDITION, SAMPLE, PATIENT))})

  observe({updateSelectInput(session, "hil_channels", choices = ch_antigen())})
  observe({updateSelectInput(session, "rad_channels", choices = ch_antigen())})

  # debounce
  # db_home_trans <- reactive({input$home_trans}) %>% debounce(3000)
  db_som_channels <- reactive({input$som_channels}) %>% debounce(3000)
  # db_hil_channels <- reactive({input$hil_channels}) %>% debounce(3000)
  # db_rad_channels <- reactive({input$rad_channels}) %>% debounce(3000)

  # extract expression data
  # expr_info <- reactive({collate_expressions(flow_proc(), min_cells = 0)}) #input$files_min)})
  # exps_info <- reactive({rescale_expressions(expr_info(), input$xp_range)})

  # update secondary inputs
  observe({
    choices <- if (input$sub_by != "__NULL__") levels(cluster()$sce@colData[[input$sub_by]]) else NONE
    updateSelectInput(session, "sub_to", choices = choices)})
  #observe({updateSelectInput(session, "rad_facet", choices = c(NONE, expr_info()$labels, metaclusts()))})

  # analyze expression data
  # mem <- reactive({calc_mem(expr_info())})
  # sim_cosine <- reactive({calc_sim_cosine(exps_info())})
  # sim_hilbert <- reactive({
  #   validate(need(!is.null(db_hil_channels()), "Select one or more channels"))
  #   calc_sim_hilbert(exps_info(), db_hil_channels(), input$hil_bins, input$hil_dims)
  # })

  # cluster data
  cluster <- reactive({
    validate(need(length(db_som_channels()) > 1, "Select two or more channels in 'Cluster' tab for clustering"))
    if (flow_not_sce()) {
      result <- do_cluster(flow_proc(), db_som_channels(), input$som_x, input$som_y, input$som_max)
      sce <- result$sce
      SummarizedExperiment::colData(sce) <- dplyr::left_join(
        sce@colData %>% as.data.frame(),
        CATALYST::cluster_codes(sce) %>% as.data.frame(),
        c("cluster_id" = "som100")) %>% S4Vectors::DataFrame()
      result$sce <- sce
      return(result)
    } else {
      sce <- CATALYST::cluster(flow_sce(), db_som_channels(), input$som_x, input$som_y, input$som_max)
      sce@colData$sample <- sce$sample_id # block MultiHeatmap error
      SummarizedExperiment::colData(sce) <- dplyr::left_join(
        sce@colData %>% as.data.frame(),
        CATALYST::cluster_codes(sce) %>% as.data.frame(),
        c("cluster_id" = "som100")) %>% S4Vectors::DataFrame()
      list(sce = sce)
    }
  })
  sce <- reactive({if (input$sub_to == "__NULL__") cluster()$sce else cluster()$sce[, c(which(cluster()$sce@colData[[input$sub_by]] == input$sub_to))]})

  # combine expression data with cluster information
  # xprc_info <- reactive({collate_expressions(flow_proc(), sce(), input$files_min)})
  # xpsc_info <- reactive({rescale_expressions(xprc_info(), input$xp_range)})

  # analyze clustered expression data
  # glmm <- reactive({compute_glmm(xprc_info(), input$som_meta, input$glmm_form, input$glmm_min) %>% compute_results(input$hypothesis)})
  # lmer <- reactive({compute_lmer(xprc_info(), input$som_meta, input$lmer_form, input$lmer_min) %>% compute_results(input$hypothesis)})
  # radviz <- reactive({
  #   validate(need(length(db_rad_channels()) > 1, "Select two or more channels"))
  #   prep_radviz(xpsc_info(), db_rad_channels())
  # })

  # display outputs
  output$home_main <- renderPrint({flow_proc()$data})
  observeEvent(input$home_source, {
    new_index <- seq_along(flow_proc()$data)
    new_data <- list()

    for (i in cache$home_view) removeTab("home_view", as.character(i))
    lapply(rev(new_index), function(i) {
      new_id <- stringr::str_c("view_", i)
      new_data[[new_id]] <- flow_proc()$data[[i]]
      prependTab("home_view", tabPanel(i, verbatimTextOutput(new_id)), TRUE)
      output[[new_id]] <- renderPrint(flow_proc()$data[[i]])
    })

    cache$home_view <- new_index
  })

  output$plot_abnd <- renderPlot({CATALYST::plotAbundances(sce(), input$som_meta, input$abnd_mode, check(input$abnd_group), check(input$abnd_shape))})
  output$plot_code <- renderPlot({CATALYST::plotCodes(sce(), input$som_meta)})
  output$plot_count <- renderPlot({CATALYST::plotCounts(sce(), input$som_meta, "condition", input$count_prop)})
  output$plot_clxp <- renderPlot({CATALYST::plotClusterExprs(sce(), input$som_meta, NULL)})
  output$plot_expr <- renderPlot({
    validate(need(!is.null(db_som_channels()), "Select one or more channels"))
    CATALYST::plotExprHeatmap(sce(), db_som_channels(), input$exph_by, input$som_meta, bars = TRUE, perc = TRUE)
  })
  output$plot_heat <- renderPlot({
    validate(need(!is.null(db_som_channels()), "Select one or more channels"))
    CATALYST::plotMultiHeatmap(sce(), db_som_channels(), input$mlth_freq, input$som_meta, bars = TRUE, perc = TRUE)
  })
  # output$plot_star <- renderPlot({FlowSOM::PlotStars(cluster()$som$FlowSOM, backgroundValues = cluster()$som$metaclustering)})

  output$plot_dr <- renderPlot({
    set.seed(1)
    sce() %>%
      CATALYST::runDR(input$dr_method, check(input$dr_cells), features = db_som_channels()) %>%
      CATALYST::plotDR(NULL, input$dr_color, check(input$dr_facet), ncol = 5)
  })
  # output$plot_file <- renderPlot({
  #   expr_info()$expr_anno %>%
  #     dplyr::select(dplyr::all_of(expr_info()$labels)) %>%
  #     dplyr::group_by(strain, file) %>%
  #     dplyr::summarise(n = dplyr::n()) %>%
  #     ggplot2::ggplot(ggplot2::aes(file, n, fill = strain)) +
  #     ggplot2::geom_col() +
  #     ggplot2::geom_hline(yintercept = input$files_min) +
  #     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
  # })

  # output$plot_mem <- renderPlot({
  #   mem() %>%
  #     ggplot2::ggplot(ggplot2::aes_string(input$mem_x, "mem",
  #                                   color = check(input$mem_color),
  #                                   group = check(input$mem_group))) +
  #     ggplot2::geom_point() +
  #     ggplot2::geom_line() +
  #     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
  # })
  # output$plot_hil <- renderPlot({
  #   sim_hilbert()$hilbert %>%
  #     ggplot2::ggplot(ggplot2::aes_string("V1", "V2", color = check(input$hil_colour), size = "cells")) +
  #     ggplot2::geom_point()
  # })
  # output$plot_hil_hier <- renderPlot({plot(hclust(sim_hilbert()$dists))})
  # output$plot_cos <- renderPlot({heatmap(sim_cosine(), symm = TRUE)})

  # output$plot_radviz <- renderPlot({
  #   prep <- radviz()
  #
  #   if (input$rad_mode == "Standard") {
  #     main <- plot(prep) + ggplot2::geom_point()
  #   } else if (input$rad_mode == "Smooth") {
  #     main <- Radviz::smoothRadviz(prep)
  #   } else if (input$rad_mode == "Contour") {
  #     main <- contour(prep) + ggplot2::geom_blank()
  #   } else if (input$rad_mode == "Hex") {
  #     main <- Radviz::hexplot(prep)
  #   }
  #
  #   main + check_facet(input$rad_facet)
  # })
  # output$plot_glmm <- renderPlot({glmm() %>% ggplot2::ggplot(ggplot2::aes_string(input$vol_x, input$vol_y, color = check(input$vol_color))) + ggplot2::geom_point() + check_facet(input$vol_facet)})
  # output$plot_lmer <- renderPlot({lmer() %>% ggplot2::ggplot(ggplot2::aes_string(input$vol_x, input$vol_y, color = check(input$vol_color))) + ggplot2::geom_point() + check_facet(input$vol_facet)})

  output$debug <- renderText({browser(); TRUE})
}
