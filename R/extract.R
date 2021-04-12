collate_expressions <- function(flow_expx, exp_label, cluster = NULL, minimum = 0, downsample = 1) {
  flow_exps %>%
    dplyr::bind_cols(exp_label, .name_repair = "minimal")
  labels$cluster <- cluster()$cluster_id
  labels %>% dplyr::left_join(CATALYST::cluster_codes(cluster()), by = c("cluster" = "som100"))
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

#' Filter and plot experiment-file pairs based on cell count
#'
#' @param exp_label output of \code{\link{extract_expressions}}
#' @param minimum optional: minimum number of cells per experiment-file pair
#' @param colour optional: bar colour
#'
#' @return
#' \code{file} vector: file names
#'
#' \code{plot} ggplot2: plot
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' exp_label <- extract_labels(flow_item$data)
#' files_inc <- extract_files(exp_label)
extract_files <- function(exp_label, minimum = 0, colour = "Strain") {
  counts <- exp_label %>%
    dplyr::select(Strain = experiment, File = file_name) %>%
    dplyr::group_by(Strain, File) %>%
    dplyr::summarise(Cells = dplyr::n())

  file <- counts %>% dplyr::filter(Cells >= minimum) %>% dplyr::pull(File) %>% as.character()
  plot <- counts %>%
    ggplot2::ggplot(ggplot2::aes(x = File, y = Cells, fill = !!rlang::sym(colour))) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = minimum, colour = "black") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))

  return(list(files = file, plot = plot))
}

extract_files2 <- function(exp_label, minimum = 0, colour = "experiment") {
  counts <- exp_label %>%
    dplyr::select(.data$experiment, .data$file_name) %>%
    dplyr::group_by(.data$experiment, .data$file_name) %>%
    dplyr::summarise(cells = dplyr::n())

  file <- counts %>% dplyr::filter(.data$cells >= minimum) %>% dplyr::pull(.data$file_name) %>% as.character()
  plot <- counts %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$file_name, y = .data$cells, fill = !!rlang::sym(colour))) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = minimum, colour = "black") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))

  return(list(files = file, plot = plot))
}

#' Extract metadata labels for expression data from flowSet
#'
#' @param flow value from \code{\link{import_fcs_path}}
#'
#' @return tibble: metadata labels
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' exp_label <- extract_labels(flow_item$data)
extract_labels <- function(flow) {
  seq_along(flow) %>% purrr::map_dfr(
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
#' @param range optional: quantile range for scaling
#'
#' @return tibble: scaled expression matrix
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' flow_exps <- scale_expressions(flow_expr)
scale_expressions <- function(flow_expr, range = c(0, 1)) {
  flow_expr %>% dplyr::transmute(dplyr::across(.fns = ~ Radviz::do.L(.x, function(x) quantile(x, range))))
}
