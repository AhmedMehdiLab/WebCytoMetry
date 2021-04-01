do_thr_iqr <- function(x) {
  th <- min(x)
  if (th == 0) {
    th <- 1 / min(x[x != 0]) > 10^c(1:10)
    th <- 1 / 10 ^ min(which(th == 0))
  }
  
  x[x == 0] <- th
  return(x)
}

# calculate marker enrichment modeling
flow_calc_mem <- function(live_mat, files_filt, flow_annos, channels, sample = 0.05) {
  library(tidyverse)
  files <- flow_annos$file
  experiments <- flow_annos$exps
  time <- flow_annos$time
  files.include <- files_filt
  
  # TODO gather is superseded
  live.df <-
    data.frame(
      live_mat,
      File = files,
      Strain = experiments,
      Time = time,
      id = seq(nrow(live_mat))
    ) %>%
    dplyr::sample_frac(sample) %>%
    tidyr::gather(key = 'Channel',
                  value = 'value',
                  tidyselect::one_of(make.names(channels))) %>%
    tibble()
  
  ## compute median per file
  live.median <- live.df %>%
    dplyr::filter(File %in% files.include) %>%
    group_by(Channel, Strain, Time, File) %>%
    summarize(value_median = median(value))
  
  ## compute IQR per file
  live.iqr <- live.df %>%
    dplyr::filter(File %in% files.include) %>%
    group_by(Channel, Strain, Time, File) %>%
    summarize(value_iqr = IQR(value)) %>%
    mutate(value_iqr = do_thr_iqr(value_iqr))
  
  ## define reference file as the file closest to medoid over all files
  live.pam <- live.median %>%
    ungroup() %>%
    select(-Strain, -Time) %>%
    spread(Channel, value_median) %>%
    select(-File) %>%
    cluster::pam(., 1)
  
  ref.live.file <- levels(live.median$File)[live.pam$id.med]
  
  ## compute median difference to reference
  live.diff <- live.median %>%
    group_by(Channel) %>%
    mutate(value_diff = value_median - value_median[File == ref.live.file]) %>%
    dplyr::filter(File != ref.live.file)
  
  ## compute IQR ratio to reference
  live.ratio <- live.iqr %>%
    group_by(Channel) %>%
    mutate(value_ratio = value_iqr / value_iqr[File == ref.live.file]) %>%
    dplyr::filter(File != ref.live.file)
  
  ## compute MEM
  live.mem <- live.diff %>%
    left_join(live.ratio,
              by = c("Channel", "Strain", "Time", "File")) %>%
    mutate(
      MEM = abs(value_diff) + 1 / value_ratio - 1,
      MEM = MEM * sign(value_diff),
      MEM = asinh(MEM)
    ) %>%
    ungroup() %>%
    mutate(
      MEM = 10 * MEM / max(abs(MEM)),
      Channel = reorder(Channel, abs(MEM), mean),
      Channel = factor(Channel, levels = rev(levels(Channel)))
    )
  
  return(list(df = live.df, mem = live.mem))
}

flow_plot_mem <- function(live_mem, colour = "Strain") {
  ggplot2::ggplot(live_mem, ggplot2::aes(x = Channel, y = MEM, color = !!sym(colour), group = File)) +
  ggplot2::geom_point() + ggplot2::geom_line() +
  ggplot2::theme_light(base_size=16) + ggplot2::ylim(-25,25) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
