
# project and visualize
# TODO provide metas
flow_proj_plot <- function(flow_annos, live_scaled, pheno_S, colour = "Strain", facet = "Strain", clusts = NULL, metas = NULL) {
  library(Radviz)
  library(tidyverse)
  files <- flow_annos[[1]]
  experiments <- flow_annos[[2]]
  time <- flow_annos[[3]]
  
  clusts <- if (is.null(clusts)) NA else clusts
  metas <- if (is.null(metas)) NA else metas
  
  colnames(live_scaled) <- make.names(colnames(live_scaled))
  pheno.rv <- do.radviz(live_scaled, pheno_S)
  #browser()
  pheno.rv.proj <-
    ggplot(data = bind_cols(
      data.frame(pheno.rv$proj$data),
      data.frame(
        Time = time,
        Strain = experiments,
        Cluster = clusts,
        Individual = files,
        Metacluster = metas
      )
    ) %>% sample_frac(0.05)) +
    geom_text(
      data = data.frame(springs(pheno.rv),
                        Channel = factor(
                          rownames(springs(pheno.rv)),
                          levels = rownames(springs(pheno.rv))
                        )),
      aes(x = X1, y = X2, label = Channel),
      col = 'orangered4'
    ) +
    scale_x_continuous(expand = c(0.3, 0.05)) +
    coord_equal() +
    theme_light(base_size = 18) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_rect(fill = "grey30")
    )
  #browser()
  densityx <- pheno.rv.proj +
    stat_density2d(
      aes(x = rx, y = ry, fill = stat(density) ^ 0.1),
      geom = "tile",
      contour = FALSE,
      n = 350
    ) +
    geom_text(
      data = data.frame(springs(pheno.rv),
                        Channel = factor(
                          rownames(springs(pheno.rv)),
                          levels = rownames(springs(pheno.rv)))),
      mapping = aes(
        x = X1,
        y = X2,
        label = Channel,
        color = !!sym(colour)
      ),
      col = 'orangered4'
    ) +
    scale_fill_continuous(low = "white", high = "dodgerblue4") +
    guides(fill = FALSE) +
    facet_grid( ~ .data[[facet]])
  
  strain <- pheno.rv.proj +
    geom_density2d(aes(x = rx, y = ry, color = !!sym(colour)), size = 2) +
    geom_point(aes(x = rx, y = ry, color = !!sym(colour))) +
    facet_grid( ~ .data[[facet]])
  
  #browser()
  meta <- pheno.rv.proj +
    geom_point(aes(x = rx, y = ry, color=Metacluster), size=0.5)+
    facet_grid(~ .data[[facet]])
  
  indiv <- pheno.rv.proj +
    geom_point(aes(x = rx, y = ry, color=Metacluster), size=0.5)+
    facet_wrap(~ Individual)
  
  return(list(density = densityx, strain = strain, meta = meta, indiv = indiv))
}

# compare experiments through 1D projections
flow_proj_explore <- function(live.df, pheno_S, files_filt) {
  files.include <- files_filt
  live.df %>%
    dplyr::filter(Channel %in% rownames(pheno_S),
                  File %in% files.include) %>%
    ggplot2::ggplot(ggplot2::aes(x = Channel, y = value)) +
    cytofan::geom_fan() +
    ggplot2::facet_grid(Strain ~ Time) +
    ggplot2::theme_light(base_size = 16) +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
