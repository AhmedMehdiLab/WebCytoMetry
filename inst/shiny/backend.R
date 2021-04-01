library(tidyverse)
library(hilbertSimilarity)
library(Radviz)

# assign channels to 'type' or 'state'
flow_panel <- function(flow_prep, type, state) {
  panel <- flow_prep
  panel$marker_class <- "none"
  panel$marker_class[panel$fcs_colname %in% type] <- "type"
  panel$marker_class[panel$fcs_colname %in% state] <- "state"
  return(panel)
}

# rows for single cells, cols for channels, scaled based on quantile
flow_matrix <- function(gate_set, index = 1,
                        scale_fun = function(x) quantile(x, c(0, 0.975)),
                        channels = NULL) {
  cyto_data <- flowWorkspace::gs_cyto_data(gate_set)[[index]]
  cyto_meta <- flowCore::pData(flowCore::parameters(cyto_data))
  
  live_mat <- flowWorkspace::lapply(gate_set, function(set) flowCore::exprs(flowWorkspace::gh_pop_get_data(set)))
  live_mat <- do.call(rbind, live_mat)
  
  colnames(live_mat) <- cyto_meta$desc
  live_mat <- live_mat[, !is.na(cyto_meta$desc)]
  if (!is.null(channels)) live_mat <- live_mat[, channels]
  
  live_scaled <- apply(live_mat, 2, Radviz::do.L, fun=scale_fun)
  return(list(matrix = live_mat, scaled = live_scaled))
}

# extract cell-level annotations
flow_annos <- function(gate_set) {
  get_anno <- function (set, cond) {
    unlist(lapply(
      flowCore::sampleNames(set), function(file) {
        rep(flowCore::pData(set)[file, cond],
            nrow(flowWorkspace::gh_pop_get_data(set[[file]])))
      }))
  }
  
  # get annotations
  files <- factor(get_anno(gate_set, "name"))
  experiments <- get_anno(gate_set, "Experiment")
  time <- factor(get_anno(gate_set, "time"))
  
  return(list(file = files, exps = experiments, time = time))
}

# plot experiments and sample sizes
flow_plot_strain <- function(flow_anno) {
  library(magrittr)
  tibble::tibble(strain = flow_anno$exps,
                 file = flow_anno$file) %>%
    dplyr::group_by(strain, file) %>%
    dplyr::summarise(N = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(order = dplyr::dense_rank(N)) %>%
    ggplot2::ggplot(ggplot2::aes(x = order, y = N, color = strain)) +
    ggplot2::geom_point() +
    ggplot2::xlab("Samples") +
    ggplot2::ylab("Cells")
}

# filter experiments based on sample size
flow_filter <- function(flow_anno, limit) {
  library(magrittr)
  tibble::tibble(strain = flow_anno$exps,
                 file = flow_anno$file) %>%
  dplyr::group_by(strain, file) %>%
  dplyr::summarise(N = n()) %>%
  dplyr::filter(N > limit) %>%
  dplyr::pull(file) %>%
  as.character()
}

# optimize channel order
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

# run FlowSOM
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

# identify populations
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