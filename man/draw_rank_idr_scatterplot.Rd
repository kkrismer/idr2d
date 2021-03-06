% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{draw_rank_idr_scatterplot}
\alias{draw_rank_idr_scatterplot}
\title{Create scatterplot of IDR values}
\usage{
draw_rank_idr_scatterplot(
  df,
  remove_na = TRUE,
  xlab = "rank in replicate 1",
  ylab = "rank in replicate 2",
  log_idr = FALSE,
  title = "rank - IDR dependence",
  color_gradient = c("rainbow", "default"),
  alpha = 1,
  max_points_shown = 2500
)
}
\arguments{
\item{df}{part of output of \code{\link{estimate_idr}}, a data frame with at
least the following named columns:
\tabular{rl}{
  \code{rank} \tab integer; rank of the peak, established by
  value column, ascending order\cr
  \code{rep_rank} \tab integer; rank of corresponding
  replicate peak.\cr
  \code{idr} \tab IDR of the peak and the
  corresponding peak in the other replicate.
}}

\item{remove_na}{logical; should NA values be removed?}

\item{xlab}{character; x axis label}

\item{ylab}{character; y axis label}

\item{log_idr}{logical; use logarithmized IDRs for colors to better
distinguish highly significant IDRs}

\item{title}{character; plot title}

\item{color_gradient}{character; either "rainbow" or "default"}

\item{alpha}{numeric; transparency of dots, from 0.0 - 1.0, where 1.0 is
completely opaque; default is 1.0}

\item{max_points_shown}{integer; default is 2500}
}
\value{
ggplot2 object; IDR rank scatterplot
}
\description{
Creates diagnostic plots to visualize the results of
\code{\link{estimate_idr}}.
}
\examples{
idr_results <- estimate_idr1d(idr2d:::chipseq$rep1_df,
                              idr2d:::chipseq$rep2_df,
                              value_transformation = "log")
draw_rank_idr_scatterplot(idr_results$rep1_df)

}
