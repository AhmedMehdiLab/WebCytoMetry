#' Attach expression labels to expression data, optionally with cluster data
#'
#' @param flow_expx output of \code{\link{extract_expressions}} or \code{\link{scale_expressions}}
#' @param exp_label output of \code{\link{extract_labels}}
#' @param cluster optional: output of \code{\link{do_cluster}}
#' @param min_file_cells optional integer: minimum cells per file
#'
#' @return
#' \code{data} tibble: annotated expression data
#'
#' \code{channels} character: channels
#'
#' \code{labels} character: labels
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' exp_label <- extract_labels(flow_item$data)
#' cluster <- do_cluster(flow_item$data, flow_item$panel, flow_item$panel$antigen)
#'
#' expr_anno <- collate_expressions(flow_expr, exp_label, cluster)
collate_expressions <- function(flow_expx, exp_label, cluster = NULL, min_file_cells = 0) {
  # if clustering is available, attach cluster data to labels
  if (!is.null(cluster)) {
    cluster_codes <- CATALYST::cluster_codes(cluster)
    names(cluster_codes)[1] <- "cluster"

    exp_label$cluster <- cluster$cluster_id
    exp_label <- exp_label %>% dplyr::left_join(cluster_codes)
  }

  # attach labels to expressions
  data <- flow_expx %>%
    dplyr::bind_cols(exp_label, .name_repair = "minimal") %>%
    dplyr::group_by(.data$strain, .data$file) %>%
    dplyr::filter(dplyr::n() >= min_file_cells) %>%
    dplyr::ungroup()

  return(list(data = data, channels = colnames(flow_expx), labels = colnames(exp_label)))
}

#' Extract expression data from flowSet
#'
#' @param flow value from \code{\link{import_fcs_path}}
#' @param panel value from \code{\link{import_fcs_path}}
#'
#' @return tibble: expression matrix
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
extract_expressions <- function(flow, panel) {
  seq_along(flow) %>%
    purrr::map_dfr(function(index) flow[[index]] %>% flowCore::exprs() %>% tibble::as_tibble()) %>%
    dplyr::rename_with(function(name) panel$antigen[panel$fcs_colname == name])
}

#' Extract metadata labels for expression data from flowSet
#'
#' @param flow value from \code{\link{import_fcs_path}}
#'
#' @return tibble: expression labels
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' exp_label <- extract_labels(flow_item$data)
extract_labels <- function(flow) {
  labels <- seq_along(flow) %>%
    purrr::map_dfr(
      function(index)
        flowCore::pData(flow)[index,] %>%
        tibble::as_tibble() %>%
        dplyr::mutate(dplyr::across(.fns = as.factor)) %>%
        dplyr::slice(rep(1, flowCore::nrow(flow[[index]])))
      )
}

#' Scale expression data
#'
#' @param flow_expr output of \code{\link{extract_expressions}}
#' @param range integer[2]: quantile range for scaling
#'
#' @return tibble: scaled expression matrix
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' flow_exps <- scale_expressions(flow_expr, c(0, 1))
scale_expressions <- function(flow_expr, range) {
  flow_expr %>% dplyr::transmute(dplyr::across(.fns = ~ Radviz::do.L(.x, function(x) quantile(x, range))))
}
