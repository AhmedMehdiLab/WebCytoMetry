#' Project expressions into a two-dimensional circular plot
#'
#' @param xpsc_info output of \code{\link{rescale_expressions}}
#' @param channels character: selected channels
#'
#' @return see \code{\link[Radviz]{do.radviz}}
#' @export
#'
#' @examples
#' flow_item <- import_fcs_path(system.file("extdata", "fcs", "KatJanin", package = "WebCytoMetry"))
#' cluster <- do_cluster(flow_item, flow_item$panel$antigen)
#' xprc_info <- collate_expressions(flow_item, cluster$sce)
#' xpsc_info <- rescale_expressions(xprc_info)
#'
#' radviz <- prep_radviz(xpsc_info, flow_item$panel$antigen)
#' plot(radviz)
prep_radviz <- function(xpsc_info, channels) {
  springs <- channels %>%
    Radviz::make.S() %>%
    Radviz::do.optimRadviz(calc_sim_cosine(xpsc_info)) %>%
    Radviz::get.optim() %>%
    Radviz::make.S()

  xpsc_info$expr_anno %>% Radviz::do.radviz(springs)
}
