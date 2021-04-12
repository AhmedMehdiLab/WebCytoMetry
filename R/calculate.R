#' Cluster flowSet according to channels
#'
#' @param flow value from \code{\link{import_fcs_path}}
#' @param panel value from \code{\link{import_fcs_path}}
#' @param channels channels to cluster
#' @param clusters optional: number of metaclusters
#'
#' @return clustered flowSet
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#'
#' cluster <- calc_cluster(flow_item$data, flow_item$panel, flow_item$panel$antigen)
calc_cluster <- function(flow, panel, channels, clusters = 20) {
  CATALYST::prepData(flow, panel, flowCore::pData(flow)) %>% CATALYST::cluster(channels, maxK = clusters)
}

#' Perform marker enrichment modeling
#'
#' @param summary output of \code{\link{calc_summary}}
#' @param downsample [0, 1) optional: sampling factor
#' @param seed optional: random seed
#'
#' @return tibble: MEM values
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' exp_label <- extract_labels(flow_item$data)
#'
#' summary <- calc_summary(flow_expr, exp_label, flow_item$panel$antigen[1])
#' mem <- calc_mem(summary)
calc_mem <- function(summary, downsample = 1, seed = 1) {
  set.seed(seed)

  # select reference file as closest to medoid
  reference <- levels(summary$file_name)[cluster::pam(summary$median, 1)$id.med]

  # calculate marker enrichment modelling
  summary %>% dplyr::mutate(
    deviation = median - median[file_name == reference],
    ratio = iqr / iqr[file_name == reference],
    mem = asinh((abs(deviation) + 1 / ratio - 1) * sign(deviation)),
    mem = 10 * mem / max(abs(mem)))
}

#' Calculate and plot cosine similarity
#'
#' @param flow_exps output of \code{\link{scale_expressions}}
#' @param exp_label output of \code{\link{extract_labels}}
#' @param files_inc value from \code{\link{extract_files}}
#'
#' @return tibble: cosine similarities
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' flow_exps <- scale_expressions(flow_expr)
#' exp_label <- extract_labels(flow_item$data)
#' files_inc <- extract_files(exp_label)
#'
#' sim_cosine <- calc_sim_cosine(flow_exps, exp_label, files_inc$files)
calc_sim_cosine <- function(flow_exps, exp_label, files_inc) {
  flow_exps %>%
    dplyr::bind_cols(exp_label, .name_repair = "minimal") %>%
    dplyr::filter(file_name %in% files_inc) %>%
    dplyr::select(-dplyr::all_of(names(exp_label))) %>%
    as.matrix() %>%
    Radviz::cosine()
}

#' Calculate and plot Hilbert Similarity
#'
#' @param flow_exps output of \code{\link{scale_expressions}}
#' @param exp_label output of \code{\link{extract_labels}}
#' @param files_inc value from \code{\link{extract_files}}
#' @param channels channels of interest
#' @param bins optional: number of bins
#' @param dimension optional: number of dimensions (> 2)
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
#' flow_exps <- scale_expressions(flow_expr)
#' exp_label <- extract_labels(flow_item$data)
#' files_inc <- extract_files(exp_label)
#'
#' sim_hilbert <- calc_sim_hilbert(flow_exps, exp_label, files_inc$files, flow_item$panel$antigen)
calc_sim_hilbert <- function(flow_exps, exp_label, files_inc, channels, bins = 10, dimension = 4) {
  # process data
  source <- flow_exps %>%
    dplyr::bind_cols(exp_label, .name_repair = "minimal") %>%
    dplyr::filter(file_name %in% files_inc)

  files <- source$file_name
  matrix <- source %>% dplyr::select(dplyr::all_of(channels)) %>% as.matrix()

  # generate and apply cutting points, then extract Hilbert indices
  precut <- matrix %>% hilbertSimilarity::make.cut(bins + 1)
  index <- matrix %>% hilbertSimilarity::do.cut(precut) %>% hilbertSimilarity::do.hilbert(bins)

  # calculate Jensen-Shannon distances and perform Multi-Dimensional Scaling
  dists <- table(index, files) %>% t() %>% hilbertSimilarity::js.dist()
  hilbert <- MASS::isoMDS(dists, k = dimension)$points %>%
    as.data.frame() %>%
    tibble::as_tibble(rownames = "File") %>%
    dplyr::mutate(
      Cells = table(files)[as.character(File)] %>% as.numeric(),
      Strain = File %>% stringr::str_extract("[A-Z][0-9]{2}|[A-Z][0-9]{1}") %>% as.factor(),
      Time = File %>% stringr::str_extract("[D]{1}[0-1]{2}") %>% as.factor())

  return(list(hilbert = hilbert, dists = stats::as.dist(dists)))
}

#' Run FlowSOM on flowSet
#'
#' @param flow value from \code{\link{import_fcs_path}}
#' @param channels channels to consider
#' @param clusters optional: number of clusters to generate
#'
#' @return FlowSOM output
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#'
#' som <- calc_som(flow_item$data, flow_item$panel$fcs_colname)
calc_som <- function(flow, channels, clusters = 20) {
  FlowSOM::FlowSOM(flow, colsToUse = channels, nClus = clusters)
}

#' Calculate summary
#'
#' @param flow_expr output of \code{\link{extract_expressions}}
#' @param exp_label output of \code{\link{extract_labels}}
#' @param selected selected channel
#'
#' @return tibble: summary
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' exp_label <- extract_labels(flow_item$data)
#'
#' summary <- calc_summary(flow_expr, exp_label, flow_item$panel$antigen[1])
calc_summary <- function(flow_expr, exp_label, selected) {
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

  flow_expr %>%
    dplyr::bind_cols(exp_label, .name_repair = "minimal") %>%
    dplyr::mutate(channel = selected, value = .data[[selected]]) %>%
    dplyr::group_by(channel, experiment, time, file_name) %>%
    dplyr::summarise(median = stats::median(value),
                     iqr = stats::IQR(value) %>% do_thr_iqr()) %>%
    dplyr::ungroup()
}
