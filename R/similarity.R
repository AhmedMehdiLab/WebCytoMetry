#' Title
#'
#' @param flow_anno
#' @param live_scaled
#' @param files_filt
#' @param channels
#' @param colour
#' @param bins
#'
#' @return
#' @export
#'
#' @examples
flow_similarity <- function(flow_anno, live_scaled, files_filt, channels, colour = "Strain", bins = 11) {
  files <- flow_anno$file
  files.include <- files_filt
  live.scaled <- live_scaled
  filter <- channels

  ceil.all.cuts <-
    make.cut(live.scaled[files %in% files.include, filter], n = bins + 1)
  ceil.all.hc <-
    do.hilbert(do.cut(live.scaled[files %in% files.include, filter], ceil.all.cuts), bins)
  ceil.all.table <-
    table(ceil.all.hc, droplevels(files[files %in% files.include]))

  ceil.all.dist <- js.dist(t(ceil.all.table))

  ceil.all.proj <- MASS::isoMDS(d = ceil.all.dist, k = 4)$points
  rownames(ceil.all.proj) <- labels(ceil.all.dist)
  colnames(ceil.all.proj) <- c('X1', 'X2', 'X3', 'X4')

  ceil.all.proj.df <- data.frame(ceil.all.proj) %>%
    tbl_df() %>%
    rename(X = X1, Y = X2) %>%
    mutate(
      File = factor(labels(ceil.all.dist)),
      N_cells = as.numeric(table(files)[as.character(File)]),
      Strain = str_extract(File, '[A-Z][0-9]{2}|[A-Z][0-9]{1}'),
      Strain = factor(Strain),
      Time = str_extract(File, '[D]{1}[0-1]{2}'),
      Time = factor(Time)
    )

  plot <- ceil.all.proj.df %>%
    dplyr::filter(File %in% files.include) %>%
    ggplot(aes(x = X, y = Y)) +
    geom_point(aes(color = !!sym(colour), size = N_cells), col = 'grey60') +
    xlab('Dim. 1') +
    ylab('Dim. 2') +
    theme_light(base_size = 16) +
    theme(aspect.ratio = 1)

  live.sim <- Radviz::cosine(live.scaled[files %in% files.include,])
  return(list(dist = ceil.all.dist, plot = plot, sim = live.sim))
}
