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
        column(4, numericInput("home_trans_c", "c", 0))
      ),
      helpText("Channel values will be transformed using the asinh function")
    ),
    mainPanel(
      verbatimTextOutput("home_info"),
      tabsetPanel(id = "home_view", type = "pills")
    )
  )
)

pane_cluster <- tabPanel(
  "Cluster",
  sidebarLayout(
    sidebarPanel(
      selectInput("cl_channels", "Channels", NULL, multiple = TRUE), hr(),
      selectInput("cl_abundance", "Abundances of", c("Cluster" = "cluster_id", "Sample" = "sample_id")),
      selectInput("heat_channels", "Expression Heatmap Channels", NULL, multiple = TRUE),
      selectInput("heat_by", "Aggregate by", c("Cluster" = "cluster_id", "Sample" = "sample_id", "Both cluster and sample" = "both")), hr(),
      selectInput("dr_method", "Reduction Method", c("Uniform Manifold Approximation and Projection" = "UMAP", "t-Distributed Stochastic Neighbor Embedding" = "TSNE", "Principal Component Analysis" = "PCA", "Multi-Dimensional Scaling" = "MDS", "Diffusion Map" = "DiffusionMap")),
      helpText("Reduction will work on channels specified in the 'Clustering' tab"),
      numericInput("dr_cells", "Cells", 50, 0),
      helpText("All cells will be used if set to 0"),
      selectInput("dr_color", "Color by", NULL),
      selectInput("dr_facet", "Facet by", NULL)
    ),
    mainPanel(tabsetPanel(
        tabPanel("Abundances", plotOutput("cl_abundance")),
        tabPanel("Codes", plotOutput("cl_code")),
        tabPanel("Cluster Expressions", plotOutput("cl_expr")),
        tabPanel("Expression Heatmap", plotOutput("heat_expr")),
        tabPanel("Frequency Heatmap", plotOutput("heat_freq")),
        tabPanel("Reduction", plotOutput("dr_view"))
    ))
  )
)

pane_filter <- tabPanel(
  "Filter",
  sidebarLayout(
    sidebarPanel(
      sliderInput("xp_range", "Quantile Range", 0, 1, c(0, 1)), br(),
      numericInput("filt_cell", "Minimum cells per sample", 0, 0)
    ),
    mainPanel(tabsetPanel(
        tabPanel("Expression", dataTableOutput("xp_matrix")),
        tabPanel("Scaled Expression", dataTableOutput("xp_scaled")),
        tabPanel("Samples per File", plotOutput("filt_view"), verbatimTextOutput("filt_file"))
    ))
  )
)

ui <- navbarPage(
  title = "WebCytoMetry",
  pane_main,
  pane_cluster,
  pane_filter,
  tabPanel("Debug", verbatimTextOutput("debug"))

  #   ),
  #     tabPanel("MEM Plot", plotOutput("mem")),
  #     tabPanel("Similarity Plot", plotOutput("sim")),
  #     tabPanel("Similarity Cluster", plotOutput("simc")),
  #     tabPanel("Similarity Heatmap", plotOutput("simp")),
  #     tabPanel("Channel Order", plotOutput("chano")),
  #     tabPanel("Channel Density", plotOutput("chand")),
  #     tabPanel("Channel Strain", plotOutput("chans")),
  #     tabPanel("Channel Metacluster", plotOutput("chanm")),
  #     tabPanel("Projection Exploration", plotOutput("proje")),
  #     tabPanel("Projection Stars", plotOutput("ps")),
  #     tabPanel("Projection View", plotOutput("px")),
  #     tabPanel("Projection Volcano", plotOutput("pv")),
  #     tabPanel("Projection Contrast", plotOutput("pc"))
  #   ))
  # )
)
