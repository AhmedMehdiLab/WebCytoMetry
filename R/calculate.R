#' Cluster flowSet according to channels
#'
#' @param flow value from \code{\link{import_fcs_path}}
#' @param panel value from \code{\link{import_fcs_path}}
#' @param channels character: channels to consider
#' @param grid_x optional integer: number of columns of self-organising map
#' @param grid_y optional integer: number of rows of self-organising map
#' @param clusters optional integer: number of clusters to generate
#'
#' @return SingleCellExperiment with added cluster information
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' cluster <- do_cluster(flow_item$data, flow_item$panel, flow_item$panel$antigen)
do_cluster <- function(flow, panel, channels, grid_x = 10, grid_y = 10, clusters = 20) {
  CATALYST::prepData(flow, panel, flowCore::pData(flow)) %>% CATALYST::cluster(channels, maxK = clusters)
}

#' Run FlowSOM on flowSet
#'
#' @param flow value from \code{\link{import_fcs_path}}
#' @param channels character: channels to consider
#' @param grid_x optional integer: number of columns of self-organising map
#' @param grid_y optional integer: number of rows of self-organising map
#' @param clusters optional integer: number of clusters to generate
#'
#' @return see \code{\link[FlowSOM]{FlowSOM}}
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' som <- do_som(flow_item$data, flow_item$panel$fcs_colname)
do_som <- function(flow, channels, grid_x = 10, grid_y = 10, clusters = 20) {
  FlowSOM::FlowSOM(flow, colsToUse = channels, nClus = clusters)
}

#' Perform marker enrichment modeling
#'
#' @param expr_anno output of \code{\link{collate_expressions}}
#'
#' @return tibble: MEM values
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' exp_label <- extract_labels(flow_item$data)
#' expr_anno <- collate_expressions(flow_expr, exp_label)
#'
#' mem <- calc_mem(expr_anno)
calc_mem <- function(expr_anno) {
  # calculate summary
  summary <- calc_summary(expr_anno, c("channel", "strain", "time", "file"))

  # calculate medoids
  medoids <- summary %>%
    dplyr::group_by(.data$channel) %>%
    dplyr::summarise(index = cluster::pam(.data$median, 1)$id.med,
                     median = .data$median[.data$index],
                     iqr = .data$iqr[.data$index])

  # calculate marker enrichment modelling
  summary %>%
    dplyr::mutate(
      deviation = .data$median - medoids$median[match(.data$channel, medoids$channel)],
      ratio = .data$iqr / medoids$iqr[match(.data$channel, medoids$channel)],
      mem = asinh((abs(.data$deviation) + 1 / .data$ratio - 1) * sign(.data$deviation))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(mem = .data$mem / max(abs(.data$mem)) * 10)
}

#' Calculate and cosine similarity
#'
#' @param exps_anno output of \code{\link{collate_expressions}}
#'
#' @return matrix: cosine similarities
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' flow_exps <- scale_expressions(flow_expr, c(0, 1))
#' exp_label <- extract_labels(flow_item$data)
#' exps_anno <- collate_expressions(flow_exps, exp_label)
#'
#' sim_cosine <- calc_sim_cosine(exps_anno)
calc_sim_cosine <- function(exps_anno) {
  exps_anno$data %>%
    dplyr::select(dplyr::all_of(exps_anno$channels)) %>%
    as.matrix() %>%
    Radviz::cosine()
}

#' Calculate and plot Hilbert Similarity
#'
#' @param exps_anno output of \code{\link{collate_expressions}}
#' @param channels character: channels of interest
#' @param bins optional integer: number of bins
#' @param dimension optional integer: number of dimensions (> 2)
#'
#' @return
#' \code{hilbert}: tibble: Hilbert similarities
#'
#' \code{dists}: distance matrix
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' flow_exps <- scale_expressions(flow_expr, c(0, 1))
#' exp_label <- extract_labels(flow_item$data)
#' exps_anno <- collate_expressions(flow_exps, exp_label)
#'
#' sim_hilbert <- calc_sim_hilbert(exps_anno, flow_item$panel$antigen)
calc_sim_hilbert <- function(exps_anno, channels, bins = 5, dimension = 4) {
  matrix <- exps_anno$data %>% dplyr::select(dplyr::all_of(channels)) %>% as.matrix()

  # generate and apply cutting points, then extract Hilbert indices
  precut <- matrix %>% hilbertSimilarity::make.cut(bins + 1)
  index <- matrix %>% hilbertSimilarity::do.cut(precut) %>% hilbertSimilarity::do.hilbert(bins)

  # calculate Jensen-Shannon distances and perform Multi-Dimensional Scaling
  dists <- table(index,  exps_anno$data$file) %>% t() %>% hilbertSimilarity::js.dist()
  hilbert <- MASS::isoMDS(dists, k = dimension)$points %>%
    as.data.frame() %>%
    tibble::as_tibble(rownames = "file") %>%
    dplyr::mutate(
      cells = table(exps_anno$data$file)[as.character(.data$file)] %>% as.numeric(),
      strain = stringr::str_extract(.data$file, "[A-Z][0-9]{2}|[A-Z][0-9]{1}") %>% as.factor(),
      time = stringr::str_match(.data$file, "D([0-1]{2})")[, 2] %>% as.integer())

  return(list(hilbert = hilbert, dists = stats::as.dist(dists)))
}

#' Summarize expression data for each channel
#'
#' @param expr_anno output of \code{\link{collate_expressions}}
#' @param grouping columns to group by
#'
#' @return tibble: summarized expression data
#'
#' @examples \dontrun{
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' exp_label <- extract_labels(flow_item$data)
#' expr_anno <- collate_expressions(flow_expr, exp_label)
#'
#' summary <- calc_summary(expr_anno, c("channel", "strain", "time", "file"))
#' }
calc_summary <- function(expr_anno, grouping) {
  # TODO: unsure what this function does, seek clarification
  do_thr_iqr <- function(x) {
    th <- min(x)
    if (th == 0) {
      th <- 1 / min(x[x != 0]) > 10^c(1:10)
      th <- 1 / 10 ^ min(which(th == 0))
    }

    x[x == 0] <- th
    return(x)
  }

  expr_anno$data %>%
    tidyr::pivot_longer(dplyr::all_of(expr_anno$channels), names_to = "channel") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping))) %>%
    dplyr::summarise(mean = mean(.data$value),
                     median = stats::median(.data$value),
                     iqr = stats::IQR(.data$value) %>% do_thr_iqr())
}
