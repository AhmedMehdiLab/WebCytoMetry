#' Project expressions into a two-dimensional circular plot
#'
#' @param flow_exps output of \code{\link{scale_expressions}}
#' @param exp_label output of \code{\link{extract_labels}}
#' @param sim_cosine output of \code{\link{calc_sim_cosine}}
#' @param channels selected channels
#' @param colour optional: dot color
#' @param facet optional: plot facet
#'
#' @return ggplot2: projection
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
#' plot_radviz <- plot_radviz(flow_exps, exp_label, sim_cosine, flow_item$panel$antigen)
plot_radviz <- function(flow_exps, exp_label, sim_cosine, channels, colour = "experiment", facet = "experiment") {
  springs <- channels %>% Radviz::make.S() %>% Radviz::do.optimRadviz(sim_cosine) %>% Radviz::get.optim() %>% Radviz::make.S()

  flow_exps %>%
    dplyr::bind_cols(exp_label, .name_repair = "minimal") %>%
    Radviz::do.radviz(springs) %>%
    plot() +
    ggplot2::geom_point(ggplot2::aes(color = !!rlang::sym(colour))) +
    ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(facet)))
}

#' Compute GLMM ADMB
#'
#' @param metaclust annotated metaclusters
#' @param files_inc output of \code{\link{extract_files}}
#' @param selected selected metaclustering
#' @param minimum optional: minimum cells per cluster
#'
#' @return tibble: GLMM ADMB results
#' @export
#'
#' @examples \dontrun{
#' TODO
#' }
calc_glmmadmb <- function(metaclust, files_inc, selected, minimum = 50) {
  # compute cluster sizes per file
  clusters <- metaclust %>%
    dplyr::filter(file_name %in% files_inc) %>%
    dplyr::group_by(experiment, time, file_name, cluster) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::filter(n >= minimum) %>% # filter out clusters with low cell numbers
    dplyr::group_by(experiment, time, file_name) %>%
    dplyr::summarise(total = sum(n), logTotal = log(total))

  # for each cluster, model changes in size using GLMM assuming distribution of
  # counts follows a negative binomial distribution
  levels(metaclust$cluster) %>%
    purrr::map_dfr(function(current) {
      result <- try(glmmADMB::glmmadmb(
        n ~ experiment + stats::offset(logTotal),
        clusters %>% dplyr::filter(cluster == current),
        family = "nbinom",
        admb.opts = glmmADMB::admbControl(shess = FALSE, noinit = FALSE)))

      if (class(result) == "try-error") NULL
      else if (nrow(summary(result)$coefficients) == 1) NULL
      else {
        summ <- summary(multcomp::glht(result, c("StrainE7=0")), test = multcomp::adjusted("none"))
        tibble::tibble(
          Contrast = names(summ$test$coefficients),
          Estimate = summ$test$coefficients,
          Std.Error = summ$test$sigma,
          pValue = summ$test$pvalues,
          Cluster = current,
          Meta = unique(metaclust[[selected]][metaclust$cluster == current]))
      }
    }) %>%
    dplyr::mutate(
      pValue = ifelse(pValue == 0, min(pValue[pValue != 0]) * 0.5, pValue),
      adjPValue = length(levels(Cluster)) * pValue,
      logPValue = -log10(pValue),
      logAdjPValue = -log10(adjPValue)
    )
}

#' Calculate contrast
#'
#' @param summary output of \code{\link{calc_summary}}
#' @param others other channels to consider
#'
#' @return tibble: contrast
#' @export
#'
#' @examples
#' \dontrun{
#' TODO
#' }
calc_contrast <- function(summary, others) {
  data <- summary %>%
    dplyr::group_by(experiment, channel, cluster) %>%
    dplyr::summarise(fzero = sum(median == 0) / length(median)) %>%
    dplyr::group_by(channel, cluster) %>%
    dplyr::summarise(fzero = min(fzero)) %>%
    dplyr::filter(fzero <= 2/3) # remove combinations of channels and clusters where 2/3 of the values are null

  # model change in signal intensity using a linear model
  levels(data$cluster) %>% purrr::map_dfr(function(current) {
    others %>% purrr::map_dfr(function(other) {
      part <- data %>% dplyr::filter(cluster == current, channel = other)
      if (nrow(part) > 6) {
        mod <- lme4::lmer(Median ~ Strain + (1|Strain), part)
        summ <- summary(multcomp::glht(result, c("StrainE7=0")), test = multcomp::adjusted("none"))
        tibble::tibble(
          Contrast = names(summ$test$coefficients),
          Estimate = summ$test$coefficients,
          Std.Error = summ$test$sigma,
          pValue = summ$test$pvalues,
          Cluster = current,
          Meta = unique(metaclust[[selected]][metaclust$cluster == current]))
      }
    })
  }) %>%
    dplyr::mutate(
      pValue = ifelse(pValue == 0, min(pValue[pValue != 0]) * 0.5, pValue),
      adjPValue = length(levels(Cluster)) * pValue,
      logPValue = -log10(pValue),
      logAdjPValue = -log10(adjPValue)
    )
}
