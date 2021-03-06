% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate.R
\name{calc_summary}
\alias{calc_summary}
\title{Summarize expression data for each channel}
\usage{
calc_summary(expr_info, grouping)
}
\arguments{
\item{expr_info}{output of \code{\link{collate_expressions}}}

\item{grouping}{columns to group by}
}
\value{
tibble: summarized expression data
}
\description{
Summarize expression data for each channel
}
\examples{
\dontrun{
flow_item <- import_fcs_path(system.file("extdata", "fcs", "KatJanin", package = "WebCytoMetry"))
expr_info <- collate_expressions(flow_item)

summary <- calc_summary(expr_info, c("channel", "strain", "time", "file"))
}
}
