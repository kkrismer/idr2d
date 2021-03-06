% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{estimate_idr1d}
\alias{estimate_idr1d}
\title{Estimates IDR for Genomic Peak Data}
\usage{
estimate_idr1d(
  rep1_df,
  rep2_df,
  value_transformation = c("identity", "additive_inverse", "multiplicative_inverse",
    "log", "log_additive_inverse"),
  ambiguity_resolution_method = c("overlap", "midpoint", "value"),
  remove_nonstandard_chromosomes = TRUE,
  max_factor = 1.5,
  jitter_factor = 1e-04,
  max_gap = -1L,
  mu = 0.1,
  sigma = 1,
  rho = 0.2,
  p = 0.5,
  eps = 0.001,
  max_iteration = 30,
  local_idr = TRUE
)
}
\arguments{
\item{rep1_df}{data frame of observations (i.e., genomic peaks) of
replicate 1, with at least the following columns (position of columns
matter, column names are irrelevant):
\tabular{rll}{
  column 1:  \tab \code{chr} \tab character; genomic location of peak -
  chromosome (e.g., \code{"chr3"})\cr
  column 2:  \tab \code{start} \tab integer; genomic location of peak -
  start coordinate\cr
  column 3:  \tab \code{end} \tab integer; genomic location of peak -
  end coordinate\cr
  column 4:  \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
  to rank the interactions
}}

\item{rep2_df}{data frame of observations (i.e., genomic peaks) of
replicate 2, with the following columns (position of columns
matter, column names are irrelevant):
\tabular{rll}{
  column 1:  \tab \code{chr} \tab character; genomic location of peak -
  chromosome (e.g., \code{"chr3"})\cr
  column 2:  \tab \code{start} \tab integer; genomic location of peak -
  start coordinate\cr
  column 3:  \tab \code{end} \tab integer; genomic location of peak -
  end coordinate\cr
  column 4:  \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
  to rank the interactions
}}

\item{value_transformation}{the values in \code{x} have to be transformed in
a way such that when ordered in descending order, more significant
interactions end up on top of the list. If the values in \code{x} are
p-values, \code{"log_additive_inverse"} is recommended. The following
transformations are supported:
\tabular{rl}{
  \code{"identity"} \tab no transformation is performed on \code{x}\cr
  \code{"additive_inverse"} \tab \code{x. = -x}\cr
  \code{"multiplicative_inverse"} \tab \code{x. = 1 / x}\cr
  \code{"log"} \tab \code{x. = log(x)}. Note: zeros are replaced by
  \code{.Machine$double.xmin}\cr
  \code{"log_additive_inverse"} \tab \code{x. = -log(x)}, recommended if
  \code{x} are p-values. Note: zeros are replaced by
  \code{.Machine$double.xmin}
}

either \code{"ascending"} (more significant
interactions have lower value in \code{value} column) or \code{"descending"}
(more significant interactions have higher value in \code{value} column)}

\item{ambiguity_resolution_method}{defines how ambiguous assignments
(when one interaction in replicate 1 overlaps with multiple interactions in
replicate 2 or vice versa)
are resolved. Available methods:
\tabular{rl}{
  \code{"value"} \tab interactions are prioritized by ascending or descending
  \code{value} column (see \code{sorting_direction}), e.g., if two
  interactions in replicate 1 overlap with one interaction in replicate 2,
  the interaction from replicate 1 is chosen which has a lower (if
  \code{sorting_direction} is \code{"ascending"}) or higher (if
  \code{"descending"}) value \cr
  \code{"overlap"} \tab the interaction pair is chosen which has the highest
  relative overlap, i.e., overlap in nucleotides of replicate 1 interaction
  anchor A and replicate 2 interaction anchor A,
  plus replicate 1 interaction anchor B and replicate 2 interaction anchor B,
  normalized by their lengths\cr
  \code{"midpoint"} \tab the interaction pair is chosen which has the
  smallest
  distance between their anchor midpoints, i.e., distance from midpoint of
  replicate 1 interaction anchor A to midpoint of
  replicate 2 interaction anchor A, plus distance from midpoint of
  replicate 1 interaction anchor B to midpoint of
  replicate 2 interaction anchor B
}}

\item{remove_nonstandard_chromosomes}{removes peaks containing
genomic locations on non-standard chromosomes using
\code{\link[GenomeInfoDb:seqlevels-wrappers]{keepStandardChromosomes}}
(default is TRUE)}

\item{max_factor}{numeric; controls the replacement values for \code{Inf}
and \code{-Inf}. \code{Inf} are replaced by \code{max(x) * max_factor} and
\code{-Inf} are replaced by \code{min(x) / max_factor}.}

\item{jitter_factor}{numeric; controls the magnitude of the noise that
is added to \code{x}. This is done to break ties in \code{x}. Set
\code{jitter_factor = NULL} for no jitter.}

\item{max_gap}{integer; maximum gap in nucleotides allowed between two
anchors for them to be considered as overlapping
(defaults to -1, i.e., overlapping anchors)}

\item{mu}{ a starting value for the mean of the reproducible component.}

\item{sigma}{ a starting value for the standard deviation of the
    reproducible component.}

\item{rho}{ a starting value for the correlation coefficient of the
    reproducible component.}

\item{p}{ a starting value for the proportion of reproducible component.}

\item{eps}{ Stopping criterion. Iterations stop when the increment of
    log-likelihood is < eps*log-likelihood, Default=0.001.}

\item{max_iteration}{integer; maximum number of iterations for
IDR estimation (defaults to 30)}

\item{local_idr}{see \code{\link[idr:est.IDR]{est.IDR}}}
}
\value{
List with three components, (\code{rep1_df}, \code{rep2_df},
and \code{analysis_type}) containing the interactions from input
data frames \code{rep1_df} and \code{rep2_df} with
the following additional columns:
\tabular{rll}{
  column 1: \tab \code{chr} \tab character; genomic location of peak -
  chromosome (e.g., \code{"chr3"})\cr
  column 2: \tab \code{start} \tab integer; genomic location of peak  -
  start coordinate\cr
  column 3: \tab \code{end} \tab integer; genomic location of peak -
  end coordinate\cr
  column 4: \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
  to rank the peaks\cr
  column 5: \tab \code{rep_value} \tab numeric; value of corresponding
  replicate peak. If no corresponding peak was found, \code{rep_value} is set
  to \code{NA}.\cr
  column 6: \tab \code{rank} \tab integer; rank of the peak, established by
  value column, ascending order\cr
  column 7: \tab \code{rep_rank} \tab integer; rank of corresponding
  replicate peak. If no corresponding peak was found, \code{rep_rank} is
  set to \code{NA}.\cr
  column 8: \tab \code{idx} \tab integer; peak index, primary key\cr
  column 9: \tab \code{rep_idx} \tab integer; specifies the index of the
  corresponding peak in the other replicate (foreign key). If no
  corresponding peak was found, \code{rep_idx} is set to \code{NA}.\cr
  column 10: \tab \code{idr} \tab IDR of the peak and the
  corresponding peak in the other replicate. If no corresponding
  peak was found, \code{idr} is set to \code{NA}.
}
}
\description{
This method estimates Irreproducible Discovery Rates (IDR) for peaks in
replicated ChIP-seq experiments.
}
\examples{
idr_results <- estimate_idr1d(idr2d:::chipseq$rep1_df,
                              idr2d:::chipseq$rep2_df,
                              value_transformation = "log")
summary(idr_results)

}
\references{
Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011) Measuring
reproducibility of high-throughput experiments. Annals of Applied
Statistics, Vol. 5, No. 3, 1752-1779.
}
