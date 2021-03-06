% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate.R
\name{calc_mem}
\alias{calc_mem}
\title{Calculate marker enrichment modeling values}
\usage{
calc_mem(expr_info)
}
\arguments{
\item{expr_info}{output of \code{\link{collate_expressions}}}
}
\value{
tibble: MEM values
}
\description{
Calculate marker enrichment modeling values
}
\examples{
flow_item <- import_fcs_path(system.file("extdata", "fcs", "KatJanin", package = "WebCytoMetry"))
expr_info <- collate_expressions(flow_item)

mem <- calc_mem(expr_info)
}
