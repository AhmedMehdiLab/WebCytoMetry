#' Project expressions into a two-dimensional circular plot
#'
#' @param exps_anno output of \code{\link{collate_expressions}}
#' @param channels character: selected channels
#'
#' @return see \code{\link[Radviz]{do.radviz}}
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' flow_exps <- scale_expressions(flow_expr, c(0, 1))
#' exp_label <- extract_labels(flow_item$data)
#' exps_anno <- collate_expressions(flow_exps, exp_label)
#'
#' plot_radviz <- plot_radviz(exps_anno, flow_item$panel$antigen)
#' plot(plot_radviz)
plot_radviz <- function(exps_anno, channels) {
  springs <- channels %>%
    Radviz::make.S() %>%
    Radviz::do.optimRadviz(calc_sim_cosine(exps_anno)) %>%
    Radviz::get.optim() %>%
    Radviz::make.S()

  exps_anno$data %>% Radviz::do.radviz(springs)
}

#' Compute GLMM ADMB
#'
#' @param expn_anno output of \code{\link{collate_expressions}}
#' @param meta "meta{n}": selected metaclustering method
#' @param min_cluster_cells optional integer: minimum cells per cluster
#'
#' @return list: GLMM ADMB results per cluster
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' exp_label <- extract_labels(flow_item$data)
#' cluster <- do_cluster(flow_item$data, flow_item$panel, flow_item$panel$antigen)
#' expn_anno <- collate_expressions(NULL, exp_label, cluster)
#'
#' glmm_list <- compute_glmm(expn_anno, "meta5")
compute_glmm <- function(expn_anno, meta, min_cluster_cells = 0) {
  # compute cluster sizes per file
  clusters <- expn_anno$data %>%
    dplyr::group_by(.data$strain, .data$time, .data$file, .data$cluster) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::filter(.data$n >= min_cluster_cells) %>%
    dplyr::group_by(.data$strain, .data$time, .data$file) %>%
    dplyr::mutate(total = sum(.data$n), logTotal = log(.data$total))

  # iterate over clusters
  sort(unique(expn_anno$data$cluster)) %>%
    purrr::set_names() %>%
    purrr::map(function(current_cluster) {
      # model changes in size using GLMM assuming distribution of counts follows a negative binomial distribution
      subset <- clusters %>% dplyr::filter(.data$cluster == current_cluster)
      result <- try(glmmADMB::glmmadmb(n ~ strain + offset(logTotal), subset, family = "nbinom", admb.opts = glmmADMB::admbControl(shess = FALSE, noinit = FALSE)))

      if (class(result) == "try-error") return(NULL)
      if (nrow(summary(result)$coefficients) == 1) return(NULL)
      return(list(tags = c("cluster" = current_cluster, "metacluster" = unique(expn_anno$data[[meta]][expn_anno$data$cluster == current_cluster])),
                  item = result))
    }) %>%
    purrr::compact()
}

#' Calculate lmer
#'
#' @param expr_anno output of \code{\link{collate_expressions}}
#' @param meta "meta{n}": selected metaclustering method
#' @param min_files optional: minimum files per cluster
#'
#' @return list: lmer results per cluster per channel
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' exp_label <- extract_labels(flow_item$data)
#' cluster <- do_cluster(flow_item$data, flow_item$panel, flow_item$panel$antigen)
#' expr_anno <- collate_expressions(flow_expr, exp_label, cluster)
#'
#' lmer_list <- compute_lmer(expr_anno, "meta6")
compute_lmer <- function(expr_anno, meta, min_files = 7) {
  # remove combinations of channels and clusters where 2/3 of the values are null
  data <- expr_anno %>%
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
      expr_anno$channels %>%
        purrr::set_names() %>%
        purrr::map(function(current_channel) {
          # model change in signal intensity using a linear model
          subset <- data %>% dplyr::filter(.data$cluster == current_cluster, .data$channel == current_channel)
          result <- try(lme4::lmer(median ~ strain + (1|strain), subset))

          if (class(result) == "try-error") return(NULL)
          if (nrow(subset) <= min_files - 1) return(NULL)
          return(list(tags = c("cluster" = current_cluster, "channel" = current_channel, "metacluster" = unique(expr_anno$data[[meta]][expr_anno$data$cluster == current_cluster])),
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
#' flow_expr <- extract_expressions(flow_item$data, flow_item$panel)
#' exp_label <- extract_labels(flow_item$data)
#' cluster <- do_cluster(flow_item$data, flow_item$panel, flow_item$panel$antigen)
#' expr_anno <- collate_expressions(flow_expr, exp_label, cluster)
#'
#' lmer_list <- compute_lmer(expr_anno, "meta6")
#' results <- collate_results(lmer_list)
collate_results <- function(results) {
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
      adj.p.v = length(levels(.data$cluster)) * .data$p.value,
      log.p.v = -log10(.data$p.value),
      loga.pv = -log10(.data$adj.p.v)
    )
}
