% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process.R
\name{do_cluster}
\alias{do_cluster}
\title{Perform clustering on parsed FCS data}
\usage{
do_cluster(flow_item, channels, grid_x = 10, grid_y = 10, clusters = 20)
}
\arguments{
\item{flow_item}{output of \code{\link{import_fcs_path}}}

\item{channels}{character: channels to consider}

\item{grid_x}{optional integer: number of columns of self-organising map}

\item{grid_y}{optional integer: number of rows of self-organising map}

\item{clusters}{optional integer: number of clusters to generate}
}
\value{
\code{sce} see \code{\link[CATALYST]{cluster}}

\code{som} see \code{\link[FlowSOM]{FlowSOM}}
}
\description{
Perform clustering on parsed FCS data
}
\examples{
flow_item <- import_fcs_path(system.file("extdata", "fcs", "KatJanin", package = "WebCytoMetry"))
flow_item <- rescale_item(flow_item, flow_item$panel$fcs_colname)

cluster <- do_cluster(flow_item, flow_item$panel$antigen)
}
