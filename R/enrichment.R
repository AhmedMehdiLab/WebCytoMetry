#' Perform marker enrichment modeling
#'
#' TODO: unsure what this function does, seek clarification
#'
#' @param expressions value from \code{\link{get_expressions}} or
#' \code{\link{scale_expressions}}
#' @param files_inc value from \code{\link{filter_files}}
#' @param meta_freq value from \code{\link{create_meta_freq}}
#' @param channel selected channel
#' @param downsample [0, 1) optional: down-sample data
#' @param seed optional: random seed
#'
#' @return MEM results
#' @export
#'
#' @examples \dontrun{0}

calc_mem <- function(flow_expr, files_inc, channel, downsample = 1, seed = 1) {
  set.seed(seed)

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

  # prepare data
  data <- flow_expr$expressions %>%
    dplyr::bind_cols(flow_expr$labels, .name_repair = "minimal") %>%
    dplyr::filter(File %in% files_inc) %>%
    dplyr::slice_sample(prop = downsample) %>%
    dplyr::transmute(Channel = channel, Value = .data[[channel]], File = file_name, Strain = experiment, Time = time) %>%
    dplyr::group_by(Channel, Strain, Time, File)
  browser()
  # calculate medians and IQR
  median <- data %>% dplyr::summarise(Median = stats::median(Value))
  iqr <- data %>% dplyr::summarise(IQR = stats::IQR(Value) %>% do_thr_iqr())

  # reference file as closest file to medoid
  pam <- median %>%
    dplyr::ungroup() %>%
    dplyr::select(-Strain, -Time) %>%
    tidyr::pivot_wider(names_from = Channel, values_from = Median) %>%
    dplyr::select(-File) %>%
    cluster::pam(1)

  reference <- levels(median$File)[pam$id.med]

  # calculate median deviation from reference
  diff <- median %>%
    dplyr::group_by(Channel) %>%
    dplyr::mutate(Deviation = Median - Median[File == reference]) %>%
    dplyr::filter(File != reference)

  # calculate IQR ratio to reference
  ratio <- iqr %>%
    dplyr::group_by(Channel) %>%
    dplyr::mutate(Ratio = IQR / IQR[File == reference]) %>%
    dplyr::filter(File != reference)

  # calculate marker enrichment modeling
  mem <- diff %>%
    dplyr::left_join(ratio, by = c("Channel", "Strain", "Time", "File")) %>%
    dplyr::mutate(MEM = asinh((abs(Deviation) + 1 / Ratio - 1) * sign(Deviation))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(MEM = 10 * MEM / max(abs(MEM)),
                  Channel = stats::reorder(Channel, abs(MEM), mean),
                  Channel = factor(Channel, levels = rev(levels(Channel))))

  return(mem)
}

#' Title
#'
#' @param live_mem
#' @param colour
#'
#' @return
#' @export
#'
#' @examples
flow_plot_mem <- function(live_mem, colour = "Strain") {
  ggplot2::ggplot(live_mem, ggplot2::aes(x = Channel, y = MEM, color = !!sym(colour), group = File)) +
    ggplot2::geom_point() + ggplot2::geom_line() +
    ggplot2::theme_light(base_size=16) + ggplot2::ylim(-25,25) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
