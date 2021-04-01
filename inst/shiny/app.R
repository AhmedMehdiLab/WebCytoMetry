library(shiny)
source("load.R")
source("backend.R")
source("hilbert_similarity.R")
source("projection.R")
source("marker_enrichment.R")

raw_flows <- flow_load(system.file("extdata", package = "WebCytoMetry"))

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("WebCytoMetry"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput("source", "Source", names(raw_flows)), br(),
      selectInput("transform", "Transform (asinh)", NULL, multiple = TRUE),
      numericInput("trans_push", "a", 1),
      numericInput("trans_scale", "b", 0.25),
      numericInput("trans_lift", "c", 0), br(),
      selectInput("type", "Type", NULL, multiple = TRUE),
      selectInput("state", "State", NULL, multiple = TRUE), br(),
      selectInput("reduce", "Reduction",
                  c("Uniform Manifold Approximation and Projection" = "UMAP",
                    "t-Distributed Stochastic Neighbor Embedding" = "TSNE",
                    "Principal Component Analysis" = "PCA",
                    "Multi-Dimensional Scaling" = "MDS",
                    "Diffusion Map" = "DiffusionMap")),
      numericInput("cells", "Cells", 500, 1, step = 1),
      selectInput("features", "Features", c("State" = "state", "Type" = "type")),
      selectInput("dr_clust", "DR Colour", NULL),
      selectInput("dr_facet", "DR Facet", NULL), br(),
      selectInput("exp_heat", "Expression Heatmap", c("sample_id", "cluster_id", "both")),
      selectInput("mult_hm2", "Multi Heatmap", NULL),
      selectInput("abund_by", "Abundance", c("sample_id", "cluster_id")), br(),
      selectInput("channels", "Channels of Interest", NULL, multiple = TRUE), br(),
      numericInput("filter", "Sample threshold", 5000, 0, step = 1), br(),
      selectInput("phenotypic", "Phenotypic Channels", NULL, multiple = TRUE),
      selectInput("functional", "Functional Channels", NULL, multiple = TRUE), br(),
    ),

    mainPanel(tabsetPanel(
      tabPanel("Info", verbatimTextOutput("info"), verbatimTextOutput("debug")),
      tabPanel("Dimensional Reduction", plotOutput("reduction")),
      tabPanel("Expression Heatmap", plotOutput("expr_heat")),
      tabPanel("Multi Heatmap", plotOutput("mult_heat")),
      tabPanel("Cluster Expressions", plotOutput("clust_exp")),
      tabPanel("Codes", plotOutput("code")),
      tabPanel("Abundances", plotOutput("abundance")),
      tabPanel("Strain Plot", plotOutput("strains")),
      tabPanel("MEM Plot", plotOutput("mem")),
      tabPanel("Similarity Plot", plotOutput("sim")),
      tabPanel("Similarity Cluster", plotOutput("simc")),
      tabPanel("Similarity Heatmap", plotOutput("simp")),
      tabPanel("Channel Order", plotOutput("chano")),
      tabPanel("Channel Density", plotOutput("chand")),
      tabPanel("Channel Strain", plotOutput("chans")),
      tabPanel("Channel Metacluster", plotOutput("chanm")),
      tabPanel("Projection Exploration", plotOutput("proje")),
      tabPanel("Projection Stars", plotOutput("ps")),
      tabPanel("Projection View", plotOutput("px")),
      tabPanel("Projection Volcano", plotOutput("pv")),
      tabPanel("Projection Contrast", plotOutput("pc"))
    ))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  # load input data
  raw_flow <- reactive({raw_flows[[input$source]]})
  raw_data <- reactive({raw_flow()$data})
  raw_pane <- reactive({raw_flow()$pane})
  raw_meta <- reactive({flowCore::pData(raw_data())})

  # update inputs
  observe({
    cols <- setNames(flowCore::colnames(raw_data()), raw_pane()$antigen)
    updateSelectInput(session, "transform", choices = cols, selected = cols)
    updateSelectInput(session, "type", choices = cols, selected = cols)
    updateSelectInput(session, "channels", choices = cols, selected = cols)
    updateSelectInput(session, "phenotypic", choices = cols, selected = cols)
    updateSelectInput(session, "funcitonal", choices = cols, selected = cols)
    updateSelectInput(session, "mult_hm2", choices = c("abundances", "state", raw_pane()$antigen))})
  observe({
    cols <- setNames(flowCore::colnames(raw_data()), raw_pane()$antigen)
    cols <- cols[!cols %in% input$type]
    updateSelectInput(session, "state", choices = cols, selected = cols)})
  observe({
    columns <- names(raw_meta())
    updateSelectInput(session, "dr_clust", choices = c("meta20", columns, raw_pane()$antigen))
    updateSelectInput(session, "dr_facet", choices = c("None" = "", columns))})

  # transform flow set
  flow <- reactive({
    flowCore::transform(raw_data(), flowCore::transformList(
      from = input$transform, tfun = flowCore::arcsinhTransform(
        a = input$trans_push,
        b = input$trans_scale,
        c = input$trans_lift)))})
  gate <- reactive({flowWorkspace::GatingSet(flow())})
  pane <- reactive({flow_panel(raw_pane(), input$type, input$state)})

  # cluster flow set
  cluster <- reactive({
    prep <- CATALYST::prepData(flow(), pane(), raw_meta())
    CATALYST::cluster(prep, input$features)})

  # compute matrix
  channels <- reactive({
    cols <- setNames(flowCore::colnames(raw_data()), raw_pane()$antigen)
    cols <- cols[cols %in% input$channels]
    names(cols)})
  phenofunc <- reactive({
    allc <- union(input$phenotypic, input$functional)
    cols <- setNames(flowCore::colnames(raw_data()), raw_pane()$antigen)
    cols <- cols[cols %in% allc]
    names(cols)})
  m_data <- reactive({flow_matrix(gate(), channels = channels())})
  matrix <- reactive({m_data()$matrix})
  scaled <- reactive({m_data()$scaled})

  # annotate channels
  annos <- reactive({flow_annos(gate())})
  files <- reactive({flow_filter(annos(), input$filter)})
  mem <- reactive({flow_calc_mem(matrix(), files(), annos(), channels())})
  sim <- reactive({flow_similarity(annos(), scaled(), files(), phenofunc())})
  chan <- reactive({flow_order_channels(scaled(), sim()$sim, phenofunc())})
  som <- reactive({flow_start_flowSOM(gate(), input$functional, input$phenotypic)})
  proj <- reactive({flow_proj_plot(annos(), scaled(), chan()$S, clusts = som()$clust, metas = som()$meta)})
  sx <- reactive({flow_FlowSOM(matrix(), files(), som()$som, annos(), phenofunc(), som()$meta, som()$clust)})

  # display output
  output$info <- renderPrint({raw_data()})
  output$debug <- renderPrint({})
  output$reduction <- renderPlot({
    reduction <- CATALYST::runDR(cluster(), input$reduce, input$cells, input$features)
    CATALYST::plotDR(reduction, input$reduce, input$dr_clust, if (input$dr_facet != "") input$dr_facet)
  })
  output$expr_heat <- renderPlot({CATALYST::plotExprHeatmap(cluster(), input$features, input$exp_heat, bars = TRUE, perc = TRUE)})
  output$mult_heat <- renderPlot({CATALYST::plotMultiHeatmap(cluster(), input$features, input$mult_hm2, bars = TRUE, perc = TRUE)})
  output$clust_exp <- renderPlot(CATALYST::plotClusterExprs(cluster(), features = input$features))
  output$code <- renderPlot(CATALYST::plotCodes(cluster()))
  output$abundance <- renderPlot(CATALYST::plotAbundances(cluster(), by = input$abund_by))
  output$strains <- renderPlot(flow_plot_strain(annos()))
  output$mem <- renderPlot(flow_plot_mem(mem()$mem))
  output$sim <- renderPlot(sim()$plot)
  output$simc <- renderPlot({plot(hclust(as.dist(as.matrix(sim()$dist)[files(), files()])))})
  output$simp <- renderPlot({heatmap(sim()$sim, symm = TRUE)})
  output$chano <- renderPlot(chan()$plot)
  output$chand <- renderPlot(proj()$density)
  output$chans <- renderPlot(proj()$strain)
  output$chanm <- renderPlot({proj()$meta})
  output$proje <- renderPlot(flow_proj_explore(mem()$df, chan()$S, files()))
  output$ps <- renderPlot(sx()$stars)
  output$px <- renderPlot(sx()$meta)
  output$pv <- renderPlot(sx()$volcano)
  output$pc <- renderPlot(sx()$contrast)}

shinyApp(ui = ui, server = server)
