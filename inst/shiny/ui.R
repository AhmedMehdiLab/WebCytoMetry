library(shiny)

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
      selectInput("cl_channels", "Channels", NULL, multiple = TRUE),
      numericInput("cl_clusters", "Max Clusters", 10, 2, 30),
      selectInput("cl_k", "Clustering", NULL), hr(),
      selectInput("cl_abundance_by", "Abundance Mode", c("Cluster" = "cluster_id", "Sample" = "sample_id")),
      selectInput("heat_channels", "Expression Heatmap Channels", NULL, multiple = TRUE),
      selectInput("heat_by", "Aggregate by", c("Cluster" = "cluster_id", "Sample" = "sample_id", "Both cluster and sample" = "both")), hr(),
      selectInput("dr_method", "Reduction Method", c("Uniform Manifold Approximation and Projection" = "UMAP", "t-Distributed Stochastic Neighbor Embedding" = "TSNE", "Principal Component Analysis" = "PCA", "Multi-Dimensional Scaling" = "MDS", "Diffusion Map" = "DiffusionMap")),
      numericInput("dr_cells", "Cells", 50, 0),
      helpText("All cells will be used if set to 0"),
      selectInput("dr_color", "Color by", NULL),
      selectInput("dr_facet", "Facet by", NULL) ),
    mainPanel(tabsetPanel(
      tabPanel("Abundances", plotOutput("cl_abundance")),
      tabPanel("Codes", plotOutput("cl_code")),
      tabPanel("Cluster Expressions", plotOutput("cl_expr")),
      tabPanel("Expression Heatmap", plotOutput("expr_heat")),
      tabPanel("Frequency Heatmap", plotOutput("freq_heat")),
      tabPanel("Reduction", plotOutput("dr_plot")),
      tabPanel("Star", plotOutput("som_star")),
      tabPanel("Volcano", plotOutput("som_volcano")) ))))

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
