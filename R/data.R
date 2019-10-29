#' Example Genomic Interaction Data Set
#'
#' This object contains genomic interactions on chromosomes 1 to 5, which could
#' be the results of HiC or ChIA-PET experiments, done in duplicates.
#'
#' @format A list with two components, the data frames \code{rep1_df} and
#' \code{rep2_df}, which have the following seven columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chr_a} \tab character; genomic location of anchor A -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2: \tab \code{start_a} \tab integer; genomic location of anchor A -
#'   start coordinate\cr
#'   column 3: \tab \code{end_a} \tab integer; genomic location of anchor A -
#'   end coordinate\cr
#'   column 4: \tab \code{chr_b} \tab character; genomic location of anchor B -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 5: \tab \code{start_b} \tab integer; genomic location of anchor B -
#'   start coordinate\cr
#'   column 6: \tab \code{end_b} \tab integer; genomic location of anchor B -
#'   end coordinate\cr
#'   column 7: \tab \code{fdr} \tab numeric; False Discovery Rate -
#'   significance of interaction
#' }
"chiapet"

#' Example Genomic Peak Data Set
#'
#' This object contains genomic peaks from two replicate ChIP-seq experiments.
#'
#' @format A list with two components, the data frames \code{rep1_df} and
#' \code{rep2_df}, which have the following four columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chr} \tab character; genomic location of peak -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2: \tab \code{start} \tab integer; genomic location of peak  -
#'   start coordinate\cr
#'   column 3: \tab \code{end} \tab integer; genomic location of peak -
#'   end coordinate\cr
#'   column 4: \tab \code{value} \tab numeric; heuristic used
#'   to rank the peaks\cr
#' }
"chipseq"

#' Example HiC data set
#'
#' This object contains data from a HiC contact map of human chromosome 1 and
#' a resolution of 2.5 * 10^6, extracted from GEO series GSE71831.
#'
#' @format A list with two components, the data frames \code{rep1_df} and
#' \code{rep2_df}, which have the following four columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chr} \tab character; genomic location of block -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2: \tab \code{region1} \tab integer; genomic location of block  -
#'   coordinate A\cr
#'   column 3: \tab \code{region2} \tab integer; genomic location of block -
#'   coordinate B\cr
#'   column 4: \tab \code{value} \tab numeric; heuristic used
#'   to rank blocks, in this case: number of reads\cr
#' }
"hic"
