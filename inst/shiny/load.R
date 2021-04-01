# constants
TISS_PTN <- "[A-Z]{3}|[A-Z][a-z][A-Z]{2}|[a-z]{6}"
EXPR_PTN <- "[A-Z][0-9]{2}|[A-Z][0-9][a-z][A-Z]{3}|[A-Z][0-9][A-Z]{2}|[A-Z][0-9]{1}"
STIM_PTN <- "[0-9]{6}"
SUBJ_PTN <- "[0-9]{3}"

flow_load <- function(flow_path) {
  flow_raws <- list()
  
  for (flow_dir in list.dirs(flow_path, full.names = F, recursive = F)) {
    fcs_files <- list.files(file.path(flow_path, flow_dir, "fcs"), full.names = T)
    path_labs <- file.path(flow_path, flow_dir, "labels.csv")
    path_mark <- file.path(flow_path, flow_dir, "markers.csv")
    
    data <- ncdfFlow::read.ncdfFlowSet(fcs_files)
    data <- ncdfFlow::as.flowSet(data)
    labs <- if (file.exists(path_labs)) readr::read_csv(path_labs)
    mark <- if (file.exists(path_mark)) readr::read_csv(path_mark)
    flow_raws[[flow_dir]] <- flow_init(data, labs, mark)
  }
  
  return(flow_raws)
}

flow_init <- function(raw_set, labels, markers, transform = TRUE) {
  # extract metadata labels
  flow_names <- flowCore::pData(raw_set)$name
  flowCore::pData(raw_set)$file_name <- flow_names
  flowCore::pData(raw_set)$condition <- flow_names
  flowCore::pData(raw_set)$sample_id <- flow_names
  flowCore::pData(raw_set)$patient_id <- NA
  flowCore::pData(raw_set)$Experiment <- stringr::str_extract(flow_names, EXPR_PTN)
  flowCore::pData(raw_set)$Stimulation <- stringr::str_extract(flow_names, STIM_PTN)
  flowCore::pData(raw_set)$Subject <- stringr::str_extract(flow_names, SUBJ_PTN)
  flowCore::pData(raw_set)$Tissue <- stringr::str_extract(flow_names, TISS_PTN)
  flowCore::pData(raw_set)$time <- 10
  
  # attach experiment labels
  if (!is.null(labels))
    for (i in seq_len(nrow(labels))) {
      match <- which(flow_names == labels$`File+Name`[i])
      flowCore::pData(raw_set)$Experiment[match] <- labels$Experiment[i]
    }
  
  # rename markers
  if (!is.null(markers))
    for (i in seq_along(raw_set)) {
      desc <- flowCore::pData(flowCore::parameters(raw_set[[i]]))$desc
      
      for (j in seq_len(nrow(markers))) {
        match <- which(desc == markers$marker_label[j])
        desc[match] <- markers$marker_name[j]
      }
      
      flowCore::pData(flowCore::parameters(raw_set[[i]]))$desc <- desc
    }
  
  return(list(data = raw_set,
              pane = tibble::tibble(fcs_colname = flowCore::colnames(raw_set),
                                    antigen = flowCore::pData(flowCore::parameters(raw_set[[1]]))$desc)))
}