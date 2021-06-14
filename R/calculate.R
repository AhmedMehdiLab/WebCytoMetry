#' Summarize expression data for each channel
#'
#' @param expr_info output of \code{\link{collate_expressions}}
#' @param grouping columns to group by
#'
#' @return tibble: summarized expression data
#'
#' @examples \dontrun{
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' expr_info <- collate_expressions(flow_item)
#'
#' summary <- calc_summary(expr_info, c("channel", "strain", "time", "file"))
#' }
calc_summary <- function(expr_info, grouping) {
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

  expr_info$expr_anno %>%
    tidyr::pivot_longer(dplyr::all_of(expr_info$channels), names_to = "channel") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping))) %>%
    dplyr::summarise(mean = mean(.data$value),
                     median = stats::median(.data$value),
                     iqr = stats::IQR(.data$value) %>% do_thr_iqr())
}

#' Calculate marker enrichment modeling values
#'
#' @param expr_info output of \code{\link{collate_expressions}}
#'
#' @return tibble: MEM values
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' expr_info <- collate_expressions(flow_item)
#'
#' mem <- calc_mem(expr_info)
calc_mem <- function(expr_info) {
  # calculate summary
  summary <- calc_summary(expr_info, c("channel", "strain", "time", "file"))

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

#' Calculate cosine similarity
#'
#' @param exps_info output of \code{\link{rescale_expressions}}
#'
#' @return matrix: cosine similarities
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' expr_info <- collate_expressions(flow_item)
#' exps_info <- rescale_expressions(expr_info)
#'
#' sim_cosine <- calc_sim_cosine(exps_info)
calc_sim_cosine <- function(exps_info) {
  exps_info$expr_anno %>%
    dplyr::select(dplyr::all_of(exps_info$channels)) %>%
    as.matrix() %>%
    Radviz::cosine()
}

#' Calculate Hilbert Similarity
#'
#' @param exps_info output of \code{\link{collate_expressions}}
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
#' expr_info <- collate_expressions(flow_item)
#' exps_info <- rescale_expressions(expr_info)
#'
#' sim_hilbert <- calc_sim_hilbert(exps_info, flow_item$panel$antigen)
calc_sim_hilbert <- function(exps_info, channels, bins = 5, dimension = 4) {
  matrix <- exps_info$expr_anno %>% dplyr::select(dplyr::all_of(channels)) %>% as.matrix()

  # generate and apply cutting points, then extract Hilbert indices
  precut <- matrix %>% hilbertSimilarity::make.cut(bins + 1)
  index <- matrix %>% hilbertSimilarity::do.cut(precut) %>% hilbertSimilarity::do.hilbert(bins)

  # calculate Jensen-Shannon distances and perform Multi-Dimensional Scaling
  dists <- table(index,  exps_info$expr_anno$file) %>% t() %>% hilbertSimilarity::js.dist()
  hilbert <- MASS::isoMDS(dists, k = dimension)$points %>%
    as.data.frame() %>%
    tibble::as_tibble(rownames = "file") %>%
    dplyr::mutate(
      cells = table(exps_info$expr_anno$file)[as.character(.data$file)] %>% as.numeric(),
      strain = stringr::str_extract(.data$file, "[A-Z][0-9]{2}|[A-Z][0-9]{1}") %>% as.factor(),
      time = stringr::str_match(.data$file, "D([0-1]{2})")[, 2] %>% as.integer()
    )

  return(list(hilbert = hilbert, dists = stats::as.dist(dists)))
}

#' Compute GLMM ADMB
#'
#' @param xprc_info output of \code{\link{collate_expressions}}
#' @param meta "meta{n}": selected metaclustering method
#' @param formula character: formula expression
#' @param min_cells optional integer: minimum cells per cluster
#'
#' @return list: GLMM ADMB results per cluster
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' cluster <- do_cluster(flow_item, flow_item$panel$antigen)
#' xprc_info <- collate_expressions(flow_item, cluster$sce)
#'
#' glmm_list <- compute_glmm(xprc_info, "meta6", "n ~ strain + offset(logTotal)")
compute_glmm <- function(xprc_info, meta, formula = NULL, min_cells = 0) {
  if (is.null(formula)) return(list())

  # compute cluster sizes per file
  clusters <- xprc_info$expr_anno %>%
    dplyr::group_by(.data$strain, .data$time, .data$file, .data$cluster) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::filter(.data$n >= min_cells) %>%
    dplyr::group_by(.data$strain, .data$time, .data$file) %>%
    dplyr::mutate(total = sum(.data$n), logTotal = log(.data$total))

  # iterate over clusters
  sort(unique(xprc_info$expr_anno$cluster)) %>%
    purrr::set_names() %>%
    purrr::map(function(current_cluster) {
      # model changes in size using GLMM assuming distribution of counts follows a negative binomial distribution
      subset <- clusters %>% dplyr::filter(.data$cluster == current_cluster)
      result <- try(glmmADMB::glmmadmb(eval(parse(text = formula)), subset, family = "nbinom", admb.opts = glmmADMB::admbControl(shess = FALSE, noinit = FALSE)))

      if (class(result) == "try-error") return(NULL)
      if (nrow(summary(result)$coefficients) == 1) return(NULL)
      return(list(tags = c("cluster" = current_cluster, "channel" = NA, "metacluster" = unique(xprc_info$expr_anno[[meta]][xprc_info$expr_anno$cluster == current_cluster])),
                  item = result))
    }) %>%
    purrr::compact()
}

#' Calculate lmer
#'
#' @param xprc_info output of \code{\link{collate_expressions}}
#' @param meta "meta{n}": selected metaclustering method
#' @param formula character: formula expression
#' @param min_files optional: minimum files per cluster
#'
#' @return list: lmer results per cluster per channel
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' cluster <- do_cluster(flow_item, flow_item$panel$antigen)
#' xprc_info <- collate_expressions(flow_item, cluster$sce)
#'
#' lmer_list <- compute_lmer(xprc_info, "meta6", "median ~ strain + (1|strain)")
compute_lmer <- function(xprc_info, meta, formula, min_files = 7) {
  if (is.null(formula)) return(list())

  # remove combinations of channels and clusters where 2/3 of the values are null
  data <- xprc_info %>%
    calc_summary(c("strain", "file", "cluster", "channel")) %>%
    dplyr::group_by(.data$strain, .data$channel, .data$cluster) %>%
    dplyr::mutate(fzero = sum(.data$median == 0) / length(.data$median)) %>%
    dplyr::group_by(.data$channel, .data$cluster) %>%
    dplyr::filter(min(.data$fzero) <= 2/3) %>%
    dplyr::ungroup()

  # iterate over clusters
  sort(unique(data$cluster)) %>%
    purrr::set_names() %>%
    purrr::map(function(current_cluster) {
      # iterate over channels
      xprc_info$channels %>%
        purrr::set_names() %>%
        purrr::map(function(current_channel) {
          # model change in signal intensity using a linear model
          subset <- data %>% dplyr::filter(.data$cluster == current_cluster, .data$channel == current_channel)
          result <- try(lme4::lmer(eval(parse(text = formula)), subset))

          if (class(result) == "try-error") return(NULL)
          if (nrow(subset) <= min_files - 1) return(NULL)
          return(list(tags = c("cluster" = current_cluster, "channel" = current_channel, "metacluster" = unique(xprc_info$expr_anno[[meta]][xprc_info$expr_anno$cluster == current_cluster])),
                      item = result))
        })
    }) %>%
    purrr::flatten() %>%
    purrr::compact()
}

#' Collate computation results
#'
#' @param results output of \code{\link{compute_glmm}} or \code{\link{compute_lmer}}
#'
#' @return tibble: results
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' cluster <- do_cluster(flow_item, flow_item$panel$antigen)
#' xprc_info <- collate_expressions(flow_item, cluster$sce)
#'
#' glmm_list <- compute_glmm(xprc_info, "meta6", "n ~ strain + offset(logTotal)")
#' lmer_list <- compute_lmer(xprc_info, "meta6", "median ~ strain + (1|strain)")
#'
#' glmm_data <- compute_results(glmm_list)
#' lmer_data <- compute_results(lmer_list)
compute_results <- function(results) {
  results %>%
    purrr::map_dfr(function(result) {
      summ <- summary(multcomp::glht(result$item, c("strainE7=0")), test = multcomp::adjusted("none"))

      # extract results
      tibble::tibble(
        contrast = names(summ$test$coefficients),
        estimate = summ$test$coefficients,
        std.err = summ$test$sigma,
        p.value = summ$test$pvalues,
        dplyr::bind_rows(result$tags)
      )
    }) %>%
    dplyr::mutate(
      # adjust results
      p.value = ifelse(.data$p.value == 0, min(.data$p.value[.data$p.value != 0]) * 0.5, .data$p.value),
      adj.p.v = length(unique(.data$cluster)) * .data$p.value,
      log.p.v = -log10(.data$p.value),
      loga.pv = -log10(.data$adj.p.v)
    )
}
