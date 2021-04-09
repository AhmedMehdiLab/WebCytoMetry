#' Collate expression data for each channel
#'
#' Extracts expression data from flowSet annotated with metadata labels
#'
#' @param flow value from \code{\link{import_fcs_path}}
#' @param panel value from \code{\link{import_fcs_path}}
#'
#' @return
#' \code{expressions} tibble: expression matrix
#'
#' \code{labels} tibble: metadata labels
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- get_expressions(flow_item$data, flow_item$panel)
get_expressions <- function(flow, panel) {
  labels <- flow %>% seq_along() %>% purrr::map_dfr(
    function(index)
      flowCore::pData(flow)[index,] %>%
      tibble::as_tibble() %>%
      dplyr::mutate(dplyr::across(.fns = as.factor)) %>%
      dplyr::slice(rep(1, flowCore::nrow(flow[[index]])))
  )

  expressions <- flow %>% seq_along() %>%
    purrr::map_dfr(function(index) flow[[index]] %>% flowCore::exprs() %>% tibble::as.tibble()) %>%
    dplyr::rename_with(function(name) panel$antigen[panel$fcs_colname == name])

  return(list(labels = labels, expressions = expressions))
}

#' Scale expression data matrix
#'
#' @param expressions output of \code{\link{create_expression_matrix}}
#' @param range optional: quantile range for scaling
#'
#' @return tibble: scaled expression matrix
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- get_expressions(flow_item$data, flow_item$panel)
#' flow_exps <- scale_expressions(flow_expr$expressions)
scale_expressions <- function(expressions, range = c(0, 1)) {
  expressions %>% dplyr::transmute(dplyr::across(.fns = ~ Radviz::do.L(.x, function(x) quantile(x, range))))
}

#' Filter and plot experiment-file pairs based on cell count
#'
#' @param labels output of \code{\link{get_expressions}}
#' @param minimum minimum number of cells per experiment-file pair
#'
#' @return
#' \code{file} vector: file names
#' \code{plot} ggplot2: plot
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "KatJanin", package = "WebCytoMetry"))
#' flow_expr <- get_expressions(flow_item$data, flow_item$panel)
#' files_inc <- filter_files(flow_expr$labels, 5000)
filter_files <- function(labels, minimum) {
  counts <- labels %>%
    dplyr::select(Strain = experiment, File = file_name) %>%
    dplyr::group_by(Strain, File) %>%
    dplyr::summarise(Cells = dplyr::n())

  file <- counts %>% dplyr::filter(Cells >= minimum) %>% dplyr::pull(File) %>% as.character()
  plot <- counts %>%
    ggplot2::ggplot(ggplot2::aes(x = File, y = Cells, fill = Strain)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = minimum, colour = "black") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))

  return(list(files = file, plot = plot))
}

#' Title
#'
#' @param live_scaled
#' @param live_sim
#' @param channels
#' @param max
#' @param top
#'
#' @return
#' @export
#'
#' @examples
flow_order_channels <- function(live_scaled, live_sim, channels, max = 1, top = 20) {
  library(Radviz)

  pheno.S <- make.S(channels)
  rownames(live_sim) <- make.names(rownames(live_sim))
  colnames(live_sim) <- make.names(colnames(live_sim))

  pheno.optim <- lapply(seq_len(max), function(i) {
    optim.cur.cells <- do.optim(pheno.S, live_sim, top = top)
    return(list(
      best = in.da(make.S(get.optim(optim.cur.cells)), live_sim),
      springs = get.optim(optim.cur.cells)
    ))
  })
  pheno.optim.best <- lapply(pheno.optim, function(x)
    x$best)

  pheno.plot <- plot(sort(unlist(pheno.optim.best)))
  pheno.S <-
    make.S(pheno.optim[[which.max(pheno.optim.best)]]$springs)
  return(list(plot = pheno.plot, S = pheno.S))
}

#' Title
#'
#' @param liveSet
#' @param functional
#' @param phenotypic
#' @param flow_type
#' @param gridSize
#' @param nClus
#'
#' @return
#' @export
#'
#' @examples
flow_start_flowSOM <- function(liveSet, functional, phenotypic, flow_type = "phenotypic", gridSize = 10, nClus = 6) {
  library(FlowSOM)
  markers <- flowCore::pData(flowCore::parameters(flowWorkspace::gs_cyto_data(liveSet)[[1]]))
  markers$type <- "NA"
  markers$type_fun <- "NA"

  markers$type[markers$name %in% functional] <- "functional"
  markers$type[markers$name %in% phenotypic] <- "phenotypic"
  #browser()
  liveSOM <- FlowSOM::FlowSOM(
    flowWorkspace::gs_cyto_data(liveSet),
    colsToUse = with(markers, name[type %in% flow_type]),
    nClus = nClus,
    xdim = gridSize,
    ydim = gridSize
  )
  clusters <- liveSOM$FlowSOM$map$mapping[, 1]
  metas <- liveSOM$metaclustering[clusters]

  return(list(som = liveSOM, clust = clusters, meta = metas))
}

#' Title
#'
#' @param live.mat
#' @param files.include
#' @param som
#' @param flow_annos
#' @param channels
#' @param metas
#' @param clusters
#'
#' @return
#' @export
#'
#' @examples
flow_FlowSOM <- function(live.mat, files.include, som, flow_annos, channels, metas, clusters) {
  files <- flow_annos[[1]]
  experiments <- factor(flow_annos[[2]])
  time <- flow_annos[[3]]
  liveSOM <- som
  set1 <- RColorBrewer::brewer.pal(9, "Set1")
  # visualize FlowSOM tree with star plot
  plot_stars <- PlotStars(liveSOM$FlowSOM,
                          backgroundValues = liveSOM$metaclustering,backgroundColor= setNames(colorRampPalette(set1[c(1:nlevels(metas))])(nlevels(metas)),
                                                                                              levels(metas)))
  ## same trick using tSNE
  # downsample to a manageable number of cells, eg 100000
  # also sample equally from each condition
  n <- 200

  tsub <- lapply(levels(experiments),function(cur.stim) {
    sample(which(experiments==cur.stim & files %in% files.include),floor(n/nlevels(experiments)))
  })
  tsub <- unlist(tsub)
  #browser()
  proj <- tsne::tsne(live.mat[tsub,channels])
  proj <- data.frame(proj,
                     Strain=experiments[tsub],
                     Cluster=clusters[tsub],
                     Meta=as.numeric(metas[tsub]))
  meta.cols <- c(colorRampPalette(set1[c(1:nlevels(metas))])(nlevels(metas)), levels(metas))
  metaxx <- proj %>%
    ggplot(aes(color=factor(Meta)))+
    geom_point(aes(x=X1,y=X2),size=2)+
    facet_grid(~Strain)+
    scale_color_manual(values=meta.cols,guide_colorbar(title="Metacluster"))+
    xlab(NULL)+
    ylab(NULL)+
    theme_light(base_size = 16)

  cluster.df <- data.frame(File=files,
                           Strain=experiments,
                           Time=time,
                           Cluster=clusters) %>%
    dplyr::filter(File %in% files.include) %>%
    group_by(Strain,Time,File,Cluster) %>%
    count() %>%
    group_by(Strain,Time,File) %>%
    mutate(total=sum(n),
           logTotal=log(total))

  nlim <- 50
  cluster.df <- cluster.df %>% dplyr::filter(n>=nlim)
  cluster.models=NULL
  cluster.models <- lapply(sort(unique(clusters)),function(cur.cluster) {
    mod <- try(glmmADMB::glmmadmb(n ~ Strain + offset(logTotal) ,
                                  data=cluster.df %>%
                                    dplyr::filter(Cluster==cur.cluster) %>%
                                    droplevels(.),
                                  zeroInflation=FALSE,
                                  family="nbinom",
                                  admb.opts=glmmADMB::admbControl(shess=FALSE,noinit=FALSE)),
               silent=TRUE)
    if(class(mod)=="try-error") {
      return(NA)
    } else {
      return(mod)
    }
  })
  names(cluster.models) <- sort(unique(clusters))
  #browser()

  cluster.cont <- lapply(sort(unique(clusters)),function(cur.cluster) {
    a=which(names(cluster.models)==cur.cluster)
    cur.model <- cluster.models[[a]]
    if(is.na(cur.model)) {
      return(NULL)
    } else if (nrow(summary(cluster.models[[a]])$coefficients)==1){
      return(NULL)
    }else{

      summ           <- summary(multcomp::glht(cur.model,
                                               linfct = c("StrainE7=0")),
                                test = multcomp::adjusted("none"))
      summ           <- data.frame(Contrast=names(summ$test$coefficients),
                                   Estimate=summ$test$coefficients,
                                   Std.Error=summ$test$sigma,
                                   pValue=summ$test$pvalues,
                                   Cluster=cur.cluster,
                                   Meta=unique(metas[clusters==cur.cluster]))
    }

    return(summ)
  })
  cluster.cont <- bind_rows(cluster.cont)

  meta.cols <- setNames(colorRampPalette(set1[c(1:nlevels(metas))])(nlevels(metas)),
                        levels(metas))

  volcano <- cluster.cont %>%
    mutate(pValue=ifelse(pValue==0,min(pValue[pValue!=0])*0.5,pValue),
           adjPValue=length(cluster.models)*pValue,
           logPValue=-log10(pValue),
           logAdjPValue=-log10(adjPValue)) %>%
    ggplot(aes(x=Estimate,
               y=logPValue,
               color=Meta))+
    geom_point(size=3)+
    scale_color_manual(values=meta.cols)+
    theme_light(base_size = 16)

  func.df <- data.frame(live.mat[,channels],
                        File=files,
                        Strain=experiments,
                        Time=time,
                        Cluster=clusters,
                        Meta=metas) %>%
    inner_join(cluster.df %>%
                 dplyr::select(Strain,File,Cluster, Time) %>% #
                 distinct(),
               c("File", "Strain", "Cluster")) %>%
    tidyr::gather(key='Channel',
                  value='value',
                  one_of(channels)) %>%
    group_by(Strain,File,Cluster,Channel) %>%
    summarize(value=median(value))

  func.df <- func.df %>%
    inner_join(func.df %>%
                 group_by(Strain,Channel,Cluster) %>%
                 summarize(fzero=sum(value==0)/length(value)) %>%
                 group_by(Channel,Cluster) %>%
                 summarize(fzero=min(fzero)) %>%
                 dplyr::filter(fzero<=6/9) %>%
                 dplyr::select(Channel,Cluster),
               by = c("Cluster", "Channel"))
  #browser()
  func.cont <- lapply(sort(unique(func.df$Cluster)),function(cur.cluster) {
    #cat('Processing',cur.cluster,'\n')
    summs <- lapply(channels, function(cur.ch) {
      #     cat('\t',cur.ch,'\n')
      df <- func.df %>%
        dplyr::filter(Cluster==cur.cluster,
                      Channel==cur.ch) %>%
        droplevels(.)
      # df
      if(nrow(df)>6) {
        #print((df))
        mod <- lme4::lmer(value ~ Strain + (1|Strain),
                          data=df)
        # print(mod)
        summ           <- summary(multcomp::glht(mod,
                                                 linfct = c("StrainE7=0")),
                                  test = multcomp::adjusted("none"))
        summ           <- data.frame(Contrast=names(summ$test$coefficients),
                                     Estimate=summ$test$coefficients,
                                     Std.Error=summ$test$sigma,
                                     pValue=summ$test$pvalues,
                                     Channel=cur.ch,
                                     Cluster=cur.cluster,
                                     Meta=unique(metas[clusters==cur.cluster]))
        return(summ)
      }
    })
    summs <- summs[lapply(summs,class)=='data.frame']
    summs <- bind_rows(summs)
    return(summs)
  })
  func.cont <- bind_rows(func.cont)

  plot_cont <- func.cont %>%
    mutate(pValue=ifelse(pValue==0,min(pValue[pValue!=0])*0.5,pValue),
           adjPValue=length(cluster.models)*pValue,
           logPValue=-log10(pValue),
           logAdjPValue=-log10(adjPValue)) %>%
    ggplot(aes(x=Estimate,
               y=logPValue,
               color=Channel))+
    geom_point(size=3)+
    theme_light(base_size = 14)

  return(list(stars = plot_stars, volcano = volcano, contrast = plot_cont, meta = metaxx))
}
