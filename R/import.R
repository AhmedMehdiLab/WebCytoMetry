#' Import FCS data from FCS folder
#'
#' A FCS folder is a folder with a `fcs` sub-directory containing `.fcs` files,
#' optionally along with `labels.csv` (with headings `File+Name` and
#' `Experiment`) and `markers.csv` (with `marker_label` and `marker_name`).
#'
#' @param fcs_path path: path to FCS folder
#' @param marker_index optional integer: index of flowFrame to use as panel source
#'
#' @importFrom magrittr %>%
#' @return parsed FCS data
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
import_fcs_path <- function(fcs_path, marker_index = 1) {
  # process FCS files
  fcs_data <- suppressWarnings(flowCore::read.flowSet(path = file.path(fcs_path, "fcs")))

  # assign metadata
  fcs_names <- flowCore::pData(fcs_data)$name
  flowCore::pData(fcs_data)$file_name <- fcs_names
  flowCore::pData(fcs_data)$condition <- fcs_names
  flowCore::pData(fcs_data)$sample_id <- fcs_names
  flowCore::pData(fcs_data)$patient_id <- fcs_names
  flowCore::pData(fcs_data)$file <- stringr::str_match(fcs_names, "^(.*)\\..*$")[, 2]
  flowCore::pData(fcs_data)$strain <- stringr::str_extract(fcs_names, "[A-Z][0-9]{2}|[A-Z][0-9][a-z][A-Z]{3}|[A-Z][0-9][A-Z]{2}|[A-Z][0-9]{1}")
  flowCore::pData(fcs_data)$time <- 10

  # process labels
  label_path <- file.path(fcs_path, "labels.csv")
  if (file.exists(label_path)) {
    labels <- readr::read_csv(label_path, col_types = "cc")

    for (i in seq_len(nrow(labels))) {
      match <- which(fcs_names == labels$`File+Name`[i])
      flowCore::pData(fcs_data)$strain[match] <- labels$Experiment[i]
    }
  }

  # process markers
  marker_path <- file.path(fcs_path, "markers.csv")
  markers <- if (file.exists(marker_path)) readr::read_csv(marker_path, col_types = "cc")
    else tibble::tibble(marker_label = character(), marker_name = character())
  markers$internal <- make.names (markers$marker_name)

  # extract parameter names
  panel <- flowCore::pData(flowCore::parameters(fcs_data[[marker_index]])) %>%
    dplyr::select(fcs_colname = .data$name, antigen = .data$desc, external = .data$desc) %>%
    tibble::tibble()
  for (i in seq_len(nrow(panel))) {
    match <- which(panel$antigen == markers$marker_label[i])
    panel$antigen[match] <- markers$internal[i]
    panel$external[match] <- markers$marker_name[i]
  }

  return(list(data = fcs_data, panel = panel))
}

#' Import all FCS data from FCS folder sub-directories
#'
#' See \code{\link{import_fcs_path}} for a description of FCS folders
#'
#' @param fcs_root path: root directory of FCS folders
#'
#' @return list: parsed FCS data
#' @export
#'
#' @examples
#' flow_data <- import_fcs_root(system.file("extdata", package = "WebCytoMetry"))
import_fcs_root <- function(fcs_root) {
  flow_data <- list()

  # iterate through sub-directories
  for (fcs_dir in list.dirs(fcs_root, full.names = FALSE, recursive = FALSE))
    flow_data[[fcs_dir]] <- import_fcs_path(file.path(fcs_root, fcs_dir))

  return(flow_data)
}
