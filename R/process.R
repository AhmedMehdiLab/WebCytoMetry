#' Transform parsed FCS data with the arcsinh function
#'
#' Expression data will be transformed with `asinh(a + b * x) + c`
#'
#' @param flow_item output of \code{\link{import_fcs_path}}
#' @param ch_colnames character: channels to transform
#' @param a optional double: parameter a
#' @param b optional double: parameter b
#' @param c optional double: parameter b
#'
#' @return parsed FCS data
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "fcs", "KatJanin", package = "WebCytoMetry"))
#' flow_item <- rescale_item(flow_item, flow_item$panel$fcs_colname)
rescale_item <- function(flow_item, ch_colnames = NULL, a = 1, b = 1, c = 0) {
  if (is.null(ch_colnames)) return(flow_item)

  flow_item$data <- flowCore::transform(
    flow_item$data,
    flowCore::transformList(ch_colnames, flowCore::arcsinhTransform("", a, b, c))
  )
  return(flow_item)
}

#' Perform clustering on parsed FCS data
#'
#' @param flow_item output of \code{\link{import_fcs_path}}
#' @param channels character: channels to consider
#' @param grid_x optional integer: number of columns of self-organising map
#' @param grid_y optional integer: number of rows of self-organising map
#' @param clusters optional integer: number of clusters to generate
#'
#' @return
#' \code{sce} see \code{\link[CATALYST]{cluster}}
#'
#' \code{som} see \code{\link[FlowSOM]{FlowSOM}}
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "fcs", "KatJanin", package = "WebCytoMetry"))
#' flow_item <- rescale_item(flow_item, flow_item$panel$fcs_colname)
#'
#' cluster <- do_cluster(flow_item, flow_item$panel$antigen)
do_cluster <- function(flow_item, channels, grid_x = 10, grid_y = 10, clusters = 20) {
  som_cols <- flow_item$panel$fcs_colname[match(channels, flow_item$panel$antigen)]

  som <- flow_item$data %>% FlowSOM::FlowSOM(colsToUse = som_cols, nClus = clusters, xdim = grid_x, ydim = grid_y)
  sce <- flow_item$data %>%
    CATALYST::prepData(flow_item$panel, flowCore::pData(flow_item$data)) %>%
    CATALYST::cluster(channels, grid_x, grid_y, clusters)

  return(list(sce = sce, som = som))
}

#' Extract expression information from parsed FCS data
#'
#' @param flow_item output of \code{\link{import_fcs_path}} or \code{\link{rescale_item}}
#' @param cluster optional: value from \code{\link{do_cluster}}
#' @param min_cells optional: minimum number of cells per file
#'
#' @return expression information
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "fcs", "KatJanin", package = "WebCytoMetry"))
#' cluster <- do_cluster(flow_item, flow_item$panel$antigen)
#'
#' xprc_info <- collate_expressions(flow_item, cluster$sce)
collate_expressions <- function(flow_item, cluster = NULL, min_cells = 0) {
  # extract expression data
  exp_values <- seq_along(flow_item$data) %>%
    purrr::map_dfr(function(index) flow_item$data[[index]] %>% flowCore::exprs() %>% tibble::as_tibble()) %>%
    dplyr::rename_with(function(name) flow_item$panel$antigen[flow_item$panel$fcs_colname == name])

  # extract expression labels
  exp_labels <- seq_along(flow_item$data) %>%
    purrr::map_dfr(
      function(index)
        flowCore::pData(flow_item$data)[index,] %>%
        tibble::as_tibble() %>%
        dplyr::mutate(dplyr::across(.fns = as.factor)) %>%
        dplyr::slice(rep(1, flowCore::nrow(flow_item$data[[index]])))
    )

  # if clustering is available, attach cluster data to labels
  if (!is.null(cluster)) {
    cluster_codes <- CATALYST::cluster_codes(cluster)
    names(cluster_codes)[1] <- "cluster"

    exp_labels$cluster <- cluster$cluster_id
    exp_labels <- exp_labels %>% dplyr::left_join(cluster_codes)
  }

  # attach labels to expressions
  expr_anno <- exp_values %>%
    dplyr::bind_cols(exp_labels, .name_repair = "minimal") %>%
    dplyr::group_by(.data$strain, .data$file) %>%
    dplyr::filter(dplyr::n() >= min_cells) %>%
    dplyr::ungroup()

  return(list(expr_anno = expr_anno, channels = colnames(exp_values), labels = colnames(exp_labels)))
}

#' Scale expression information
#'
#' @param expr_info output of \code{\link{collate_expressions}}
#' @param range optional [0, 1): quantile range
#'
#' @return scaled expression information
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "fcs", "KatJanin", package = "WebCytoMetry"))
#' expr_info <- collate_expressions(flow_item)
#' exps_info <- rescale_expressions(expr_info)
rescale_expressions <- function(expr_info, range = c(0, 1)) {
  expr_info$expr_anno <- expr_info$expr_anno %>%
    dplyr::mutate(dplyr::across(expr_info$channels, ~ Radviz::do.L(.x, function(x) quantile(x, range))))

  return(expr_info)
}
