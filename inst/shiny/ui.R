library(shiny)
library(shinyjs)

NONE <- c("None" = "__NULL__")
CLUSTER <- c("Cluster" = "cluster_id")
CONDITION <- c("Condition" = "condition")
PATIENT <- c("Patient" = "patient_id")
SAMPLE <- c("Sample" = "sample_id")

COMMON_OPTS <- c("Channel" = "channel", "File" = "file", "Strain" = "strain", "Time" = "time")
COMPUTE_VAL <- c("Estimate" = "estimate", "Standard Error" = "std.err", "P-value" = "p.value", "Adjusted P-value" = "adj.p.v", "Log P-value" = "log.p.v", "Log Adjusted P-value" = "loga.pv")
COMPUTE_EXT <- c("Cluster" = "cluster", "Metacluster" = "metacluster", "Channel" = "channel")

pane_main <- tabPanel(
  "Home",
  shinyjs::useShinyjs(),
  sidebarLayout(
    sidebarPanel(
      selectInput("home_source", "Flow Cytometry Project", NULL),
      selectInput("home_sce", "Single Cell Experiment", NULL), br(),

      # selectInput("home_trans", "Transform Channels", NULL, multiple = TRUE),
      # fluidRow(
      #   column(4, numericInput("trans_a", "a", 1)),
      #   column(4, numericInput("trans_b", "b", 0.25)),
      #   column(4, numericInput("trans_c", "c", 0))),
      # helpText("Channel values will be transformed using the asinh function"), br(),

      fileInput("home_upload", "Upload Flow Cytometry RDS", accept = ".rds"),
      fileInput("home_up_sce", "Upload Single Cell Experiment RDS", accept = ".rds")),
    mainPanel(
      verbatimTextOutput("home_main"), hr(),
      tabsetPanel(id = "home_view", type = "pills")) ))

pane_cluster <- tabPanel(
  "Cluster",
  sidebarLayout(
    sidebarPanel(
      selectInput("som_channels", "Channels", NULL, multiple = TRUE),
      fluidRow(column(4, numericInput("som_x", "SOM X", 10, 2, 20)),
               column(4, numericInput("som_y", "SOM Y", 10, 2, 20)),
               column(4, numericInput("som_max", "Metaclusters", 10, 2, 20))),
      selectInput("som_meta", "Metacluster", NULL), br(),

      selectInput("sub_by", "Subset by", c(NONE, CONDITION, PATIENT, SAMPLE)),
      selectInput("sub_to", "Subset as", NULL), br(),

      selectInput("abnd_mode", "Abundance", c(CLUSTER, SAMPLE)),
      fluidRow(column(6, selectInput("abnd_group", "Group", c(NONE, CONDITION, PATIENT))),
               column(6, selectInput("abnd_shape", "Shape", c(NONE, CONDITION, PATIENT)))),
      checkboxInput("count_prop", "Proportional Count"), br(),

      selectInput("exph_by", "Expression Group", c("sample_id", "cluster_id", "both")),
      selectInput("mlth_freq", "Multi Heatmap Frequency", NULL)),
    mainPanel(tabsetPanel(
      tabPanel("Abundances", plotOutput("plot_abnd")),
      tabPanel("Codes", plotOutput("plot_code")),
      tabPanel("Counts", plotOutput("plot_count")),
      tabPanel("Expressions", plotOutput("plot_clxp")),
      tabPanel("Expression Heatmap", plotOutput("plot_expr")),
      tabPanel("Multi Heatmap", plotOutput("plot_heat"))
      # tabPanel("Star Plot", plotOutput("plot_star"))
))))

pane_reduce <- tabPanel(
  "Reduce",
  sidebarLayout(
    sidebarPanel(
      selectInput("dr_method", "Reduction Method", c("Uniform Manifold Approximation and Projection" = "UMAP", "t-Distributed Stochastic Neighbor Embedding" = "TSNE", "Principal Component Analysis" = "PCA", "Multi-Dimensional Scaling" = "MDS", "Diffusion Map" = "DiffusionMap")),
      fluidRow(column(6, selectInput("dr_color", "Color by", NULL, multiple = TRUE)),
               column(6, selectInput("dr_facet", "Facet by", NULL))),
      numericInput("dr_cells", "Cells Per Sample", 50, 0),
      helpText("All cells will be used if set to 0")),
    mainPanel(plotOutput("plot_dr")) ))

pane_filter <- tabPanel(
  "Filter",
  sidebarLayout(
    sidebarPanel(
      numericInput("files_min", "Minimum Cells Per File", 0, 0), br(),
      sliderInput("xp_range", "Quantile Range", 0, 1, c(0, 1))),
    mainPanel(plotOutput("plot_file")) ))

pane_analyze <- tabPanel(
  "Analyze",
  sidebarLayout(
    sidebarPanel(
      selectInput("mem_x", "MEM Plot", COMMON_OPTS),
      fluidRow(column(6, selectInput("mem_color", "Color", c(NONE, COMMON_OPTS))),
               column(6, selectInput("mem_group", "Group", c(NONE, COMMON_OPTS)))), br(),

      selectInput("hil_channels", "Hilbert Similarity", NULL, multiple = TRUE),
      fluidRow(column(6, numericInput("hil_bins", "Bins", 5, 3, 20)),
               column(6, numericInput("hil_dims", "Dimensions", 4, 2, 10))),
      selectInput("hil_colour", "Color", c(NONE, "File" = "file", "Strain" = "strain", "Time" = "time")),
      helpText("Only the first two dimensions will be plotted")),
    mainPanel(tabsetPanel(
      tabPanel("MEM Plot", plotOutput("plot_mem")),
      tabPanel("Hilbert Similarity", plotOutput("plot_hil")),
      tabPanel("Hilbert Hierarchy", plotOutput("plot_hil_hier")),
      tabPanel("Cosine Similarity", plotOutput("plot_cos")) ))))

pane_project <- tabPanel(
  "Project",
  sidebarLayout(
    sidebarPanel(
      selectInput("rad_channels", "Radviz Channels", NULL, multiple = TRUE),
      fluidRow(column(6, selectInput("rad_mode", "Mode", c("Standard", "Smooth", "Contour", "Hex"))),
               column(6, selectInput("rad_facet", "Facet", NULL))), br(),

      fluidRow(column(6, selectInput("vol_x", "Volcano X", COMPUTE_VAL)),
               column(6, selectInput("vol_y", "Volcano Y", COMPUTE_VAL))),
      fluidRow(column(6, selectInput("vol_color", "Color", c(NONE, COMPUTE_EXT))),
               column(6, selectInput("vol_facet", "Facet", c(NONE, COMPUTE_EXT)))), br(),

      textInput("glmm_form", "GLMM Formula", "n ~ strain + offset(logTotal)"),
      textInput("lmer_form", "LMER Formula", "median ~ strain + (1|strain)"),
      textInput("hypothesis", "Hypothesis", "strainE7 = 0"), br(),

      numericInput("glmm_min", "GLMM Minimum Cells Per File", 50, 0),
      numericInput("lmer_min", "LMER Minimum Files Per Cluster", 5, 0)),
    mainPanel(tabsetPanel(
      tabPanel("Radviz Plot", plotOutput("plot_radviz")),
      tabPanel("GLMM Volcano", plotOutput("plot_glmm")),
      tabPanel("LMER Volcano", plotOutput("plot_lmer")) ))))

ui <- navbarPage("WebCytoMetry", pane_main, pane_cluster, pane_reduce
                 # pane_filter, pane_analyze, pane_project,
                 # , tabPanel("Debug", verbatimTextOutput("debug"))
)
