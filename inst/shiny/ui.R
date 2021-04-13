library(shiny)

NONE <- c("None" = "__NULL__")
CLUSTER <- c("Cluster" = "cluster_id")
CONDITION <- c("Condition" = "condition")
PATIENT <- c("Patient" = "patient_id")
SAMPLE <- c("Sample" = "sample_id")

pane_main <- tabPanel(
  "Home",
  sidebarLayout(
    sidebarPanel(
      selectInput("home_source", "Source", NULL), br(),
      selectInput("home_trans", "Transform Channels", NULL, multiple = TRUE),
      fluidRow(
        column(4, numericInput("home_trans_a", "a", 1)),
        column(4, numericInput("home_trans_b", "b", 0.25)),
        column(4, numericInput("home_trans_c", "c", 0))),
      helpText("Channel values will be transformed using the asinh function")),
    mainPanel(
      verbatimTextOutput("home_main"),
      tabsetPanel(id = "home_view", type = "pills")) ))

pane_cluster <- tabPanel(
  "Cluster",
  sidebarLayout(
    sidebarPanel(
      selectInput("som_channels", "Channels", NULL, multiple = TRUE),
      fluidRow(column(4, numericInput("som_x", "SOM X", 10, 2, 20)),
               column(4, numericInput("som_y", "SOM Y", 10, 2, 20)),
               column(4, numericInput("som_max", "Metaclusters", 10, 2, 20))),
      selectInput("som_meta", "Metacluster", NULL), hr(),

      selectInput("abnd_mode", "Abundance", c(CLUSTER, SAMPLE)),
      fluidRow(column(6, selectInput("abnd_group", "Group", c(NONE, CONDITION, PATIENT))),
               column(6, selectInput("abnd_shape", "Shape", c(NONE, CONDITION, PATIENT)))), hr(),

      selectInput("heat_expr", "Heatmap Expression", NULL, multiple = TRUE),
      selectInput("heat_freq", "Heatmap Frequency", NULL), hr(),

      selectInput("dr_method", "Reduction Method", c("Uniform Manifold Approximation and Projection" = "UMAP", "t-Distributed Stochastic Neighbor Embedding" = "TSNE", "Principal Component Analysis" = "PCA", "Multi-Dimensional Scaling" = "MDS", "Diffusion Map" = "DiffusionMap")),
      fluidRow(column(6, selectInput("dr_color", "Color by", NULL)),
               column(6, selectInput("dr_facet", "Facet by", NULL))),
      numericInput("dr_cells", "Cells Per Sample", 50, 0),
      helpText("All cells will be used if set to 0")),
    mainPanel(tabsetPanel(
      tabPanel("Abundances", plotOutput("plot_abnd")),
      tabPanel("Codes", plotOutput("plot_code")),
      tabPanel("Cluster Expressions", plotOutput("plot_clxp")),
      tabPanel("Heatmap", plotOutput("plot_heat")),
      tabPanel("Reduction", plotOutput("plot_dr")),
      tabPanel("Star", plotOutput("plot_star")) ))))

pane_filter <- tabPanel(
  "Filter",
  sidebarLayout(
    sidebarPanel(
      sliderInput("xp_range", "Quantile Range", 0, 1, c(0, 1)), br(),
      numericInput("files_min", "Minimum cells per sample", 0, 0)),
    mainPanel(tabsetPanel(
        tabPanel("Expression", dataTableOutput("xp_matrix")),
        tabPanel("Scaled Expression", dataTableOutput("xp_scaled")),
        tabPanel("Samples per File", plotOutput("file_plot"), verbatimTextOutput("file_view")) ))))

pane_analyze <- tabPanel(
  "Analyze",
  sidebarLayout(
    sidebarPanel(
      selectInput("mem_channels", "MEM Channel", NULL),
      sliderInput("mem_sample", "MEM Downsample", 0, 1, 1),
      selectInput("mem_colour", "MEM Colour", c("Strain", "Time")),
      selectInput("mem_background", "MEM Contrast Background", NULL, multiple = TRUE), hr(),
      selectInput("sim_channels", "Hilbert Channels", NULL, multiple = TRUE),
      selectInput("sim_colour", "Hilbert Colour", c("Strain", "Time")),
      numericInput("sim_bins", "Hilbert Bins", 3, 20, 5),
      numericInput("sim_dims", "Hilbert Dimensions", 2, 10, 4),
      helpText("Only the first two dimensions will be plotted"), hr(),
      selectInput("ch_channels", "Projection Channels", NULL, multiple = TRUE),
      selectInput("ch_colour", "Projection Colour", NULL),
      selectInput("ch_facet", "Projection Facet", NULL)),
    mainPanel(tabsetPanel(
      tabPanel("MEM Table", dataTableOutput("mem_view")),
      tabPanel("MEM Plot", plotOutput("mem_plot")),
      tabPanel("MEM Contrast", plotOutput("mem_volcano")),
      tabPanel("Hilbert Similarity", dataTableOutput("simh_view")),
      tabPanel("Hilbert Plot", plotOutput("simh_plot")),
      tabPanel("Hilbert Hierarchy", plotOutput("simh_hier")),
      tabPanel("Cosine Similarity", dataTableOutput("simc_view")),
      tabPanel("Cosine Heatmap", plotOutput("simc_heat")),
      tabPanel("Projection", plotOutput("projection")) ))))

ui <- navbarPage("WebCytoMetry", pane_main, pane_cluster, pane_filter, pane_analyze
  , tabPanel("Debug", verbatimTextOutput("debug"))
)
