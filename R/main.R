
#' @title Finds One-to-One Correspondence between Peaks from Replicate 1
#' and 2
#'
#' @description
#' This method establishes a bijective assignment between peaks from
#' replicate 1 and 2. A peak in replicate 1 is assigned to a
#' peak in replicate 2 if and only if (1) they overlap (or the gap between the
#' peaks is less than or equal to \code{max_gap}), and (2) there is no other
#' peak in
#' replicate 2 that overlaps with the peak in replicate 1 and has a
#' lower \emph{ambiguity resolution value}.
#'
#' @inheritParams establish_overlap1d
#'
#' @return Data frames \code{rep1_df} and \code{rep2_df} with
#' the following columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chr} \tab character; genomic location of peak -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2: \tab \code{start} \tab integer; genomic location of peak  -
#'   start coordinate\cr
#'   column 3: \tab \code{end} \tab integer; genomic location of peak -
#'   end coordinate\cr
#'   column 4: \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the peaks\cr
#'   column 5: \tab \code{rep_value} \tab numeric; value of corresponding
#'   replicate peak. If no corresponding peak was found, \code{rep_value} is set
#'   to \code{NA}.\cr
#'   column 6: \tab \code{rank} \tab integer; rank of the peak, established by
#'   value column, ascending order\cr
#'   column 7: \tab \code{rep_rank} \tab integer; rank of corresponding
#'   replicate peak. If no corresponding peak was found, \code{rep_rank} is
#'   set to \code{NA}.\cr
#'   column 8: \tab \code{idx} \tab integer; peak index, primary key\cr
#'   column 9: \tab \code{rep_idx} \tab integer; specifies the index of the
#'   corresponding peak in the other replicate (foreign key). If no
#'   corresponding peak was found, \code{rep_idx} is set to \code{NA}.
#' }
#'
#' @examples
#' rep1_df <- idr2d:::chipseq$rep1_df
#' rep1_df$value <- preprocess(rep1_df$value, "log")
#'
#' rep2_df <- idr2d:::chipseq$rep2_df
#' rep2_df$value <- preprocess(rep2_df$value, "log")
#'
#' mapping <- establish_bijection1d(rep1_df, rep2_df)
#'
#' @export
establish_bijection1d <- function(rep1_df, rep2_df,
                                  ambiguity_resolution_method = c("overlap",
                                                                  "midpoint",
                                                                  "value"),
                                  max_gap = 0L) {
    return(establish_bijection(rep1_df, rep2_df, "IDR1D",
                               ambiguity_resolution_method))
}

#' @title Finds One-to-One Correspondence between Interactions from Replicate 1
#' and 2
#'
#' @description
#' This method establishes a bijective assignment between interactions from
#' replicate 1 and 2. An interaction in replicate 1 is assigned to an
#' interaction in replicate 2 if and only if (1) both anchors of the
#' interactions
#' overlap (or the gap between anchor A/B in replicate 1 and 2 is less than
#' or equal to \code{max_gap}), and (2) there is no other interaction in
#' replicate 2 that overlaps with the interaction in replicate 1 and has a
#' lower \emph{ambiguity resolution value}.
#'
#' @inheritParams establish_overlap2d
#'
#' @return Data frames \code{rep1_df} and \code{rep2_df} with
#' the following columns:
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
#'   column 7: \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the interactions\cr
#'   column 8: \tab \code{"rep_value"} \tab numeric; value of corresponding
#'   replicate interaction. If no corresponding interaction was found,
#'   \code{rep_value} is set to \code{NA}.\cr
#'   column 9: \tab \code{"rank"} \tab integer; rank of the interaction,
#'   established by value column, ascending order\cr
#'   column 10: \tab \code{"rep_rank"} \tab integer; rank of corresponding
#'   replicate interaction. If no corresponding interaction was found,
#'   \code{rep_rank} is set to \code{NA}.\cr
#'   column 11: \tab \code{"idx"} \tab integer; interaction index,
#'   primary key\cr
#'   column 12: \tab \code{"rep_idx"} \tab integer; specifies the index of the
#'   corresponding interaction in the other replicate (foreign key). If no
#'   corresponding interaction was found, \code{rep_idx} is set to \code{NA}.
#' }
#'
#' @examples
#' rep1_df <- idr2d:::chiapet$rep1_df
#' rep1_df$fdr <- preprocess(rep1_df$fdr, "log_additive_inverse")
#'
#' rep2_df <- idr2d:::chiapet$rep2_df
#' rep2_df$fdr <- preprocess(rep2_df$fdr, "log_additive_inverse")
#'
#' mapping <- establish_bijection2d(rep1_df, rep2_df)
#'
#' @export
establish_bijection2d <- function(rep1_df, rep2_df,
                                  ambiguity_resolution_method = c("overlap",
                                                                  "midpoint",
                                                                  "value"),
                                  max_gap = 0L) {
    return(establish_bijection(rep1_df, rep2_df, "IDR2D",
                               ambiguity_resolution_method))
}

#' @title Finds One-to-One Correspondence between Peaks or interactions
#' from Replicate 1 and 2
#'
#' @description
#' This method establishes a bijective assignment between observations
#' (genomic peaks in case of ChIP-seq, genomic interactions in case of
#' ChIA-PET, HiChIP, and HiC) from
#' replicate 1 and 2. An observation in replicate 1 is assigned to an
#' observation in replicate 2 if and only if (1) the observation loci in both
#' replicates overlap (or the gap between them is less than
#' or equal to \code{max_gap}), and (2) there is no other observation in
#' replicate 2 that overlaps with the observation in replicate 1 and has a
#' lower \emph{ambiguity resolution value}.
#'
#'@param rep1_df data frame of observations (i.e., genomic peaks or genomic
#' interactions) of
#' replicate 1. If \code{analysis_type} is IDR1D, the columns of \code{rep1_df}
#' are described in \code{\link{establish_bijection1d}}, otherwise in
#' \code{\link{establish_bijection2d}}
#' @param rep2_df data frame of observations (i.e., genomic peaks or genomic
#' interactions) of replicate 2. Same columns as \code{rep1_df}.
#' @param analysis_type "IDR2D" for genomic interaction data sets,
#' "IDR1D" for genomic peak data sets
#' @param ambiguity_resolution_method defines how ambiguous assignments
#' (when one interaction or peak in replicate 1 overlaps with
#' multiple interactions or peaks in replicate 2 or vice versa)
#' are resolved. For available methods, see \code{\link{establish_overlap1d}} or
#' \code{\link{establish_overlap2d}}, respectively.
#'
#' @inheritParams determine_anchor_overlap
#'
#' @return See \code{\link{establish_bijection1d}} or
#' \code{\link{establish_bijection2d}}, respectively.
#'
#' @examples
#' rep1_df <- idr2d:::chipseq$rep1_df
#' rep1_df$value <- preprocess(rep1_df$value, "log")
#'
#' rep2_df <- idr2d:::chipseq$rep2_df
#' rep2_df$value <- preprocess(rep2_df$value, "log")
#'
#' mapping <- establish_bijection(rep1_df, rep2_df, analysis_type = "IDR1D")
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' @importFrom dplyr select
#' @importFrom futile.logger flog.warn
#' @importFrom magrittr "%>%"
#' @export
establish_bijection <- function(rep1_df, rep2_df,
                                analysis_type = c("IDR1D", "IDR2D"),
                                ambiguity_resolution_method = c("overlap",
                                                                "midpoint",
                                                                "value"),
                                max_gap = 0L) {
    # avoid CRAN warnings
    rep1_idx <- rep2_idx <- arv <- value <- NULL

    # argument handling
    analysis_type <- match.arg(analysis_type, choices = c("IDR1D", "IDR2D"))

    if (analysis_type == "IDR1D") {
        columns <- c("chr", "start", "end", "value")
    } else if (analysis_type == "IDR2D") {
        columns <- c("chr_a", "start_a", "end_a",
                     "chr_b", "start_b", "end_b", "value")
    } else {
        stop("unknown analysis type")
    }

    # remove irrelevant columns, rename relevant columns
    rep1_df <- rep1_df[, seq_len(length(columns))]
    rep2_df <- rep2_df[, seq_len(length(columns))]
    colnames(rep1_df) <- columns
    colnames(rep2_df) <- columns

    if (nrow(rep1_df) > 0 && nrow(rep2_df) > 0) {
        # shuffle to break preexisting order
        rep1_df <- rep1_df[sample.int(nrow(rep1_df)), ]
        rep2_df <- rep2_df[sample.int(nrow(rep2_df)), ]

        # sort by value column, descending order (higher values are better)
        rep1_df <- dplyr::arrange(rep1_df, dplyr::desc(value))
        rep2_df <- dplyr::arrange(rep2_df, dplyr::desc(value))

        # add idx column
        rep1_df$idx <- seq_len(nrow(rep1_df))
        rep2_df$idx <- seq_len(nrow(rep2_df))

        # add rank column
        rep1_df$rank <- seq_len(nrow(rep1_df))
        rep2_df$rank <- seq_len(nrow(rep2_df))

        if (analysis_type == "IDR1D") {
            pairs_df <- establish_overlap1d(rep1_df, rep2_df,
                                            ambiguity_resolution_method,
                                            max_gap)
        } else if (analysis_type == "IDR2D") {
            pairs_df <- establish_overlap2d(rep1_df, rep2_df,
                                            ambiguity_resolution_method,
                                            max_gap)
        }

        top_pairs_df <- pairs_df %>% dplyr::group_by(rep1_idx) %>%
            dplyr::slice(which.min(arv))
        top_pairs_df <- top_pairs_df %>% dplyr::group_by(rep2_idx) %>%
            dplyr::slice(which.min(arv))

        if (nrow(top_pairs_df) != length(unique(top_pairs_df$rep1_idx)) &&
            nrow(top_pairs_df) != length(unique(top_pairs_df$rep2_idx))) {
            futile.logger::flog.warn("ambiguous interaction assignments")
        }

        rep1_df$rep_idx <- as.integer(NA)
        rep2_df$rep_idx <- as.integer(NA)

        rep1_df$rep_idx[top_pairs_df$rep1_idx] <- top_pairs_df$rep2_idx
        rep2_df$rep_idx[top_pairs_df$rep2_idx] <- top_pairs_df$rep1_idx

        rep1_df$rep_rank <- as.integer(NA)
        rep2_df$rep_rank <- as.integer(NA)

        rep1_df$rep_rank[top_pairs_df$rep1_idx] <-
            rep2_df$rank[top_pairs_df$rep2_idx]
        rep2_df$rep_rank[top_pairs_df$rep2_idx] <-
            rep1_df$rank[top_pairs_df$rep1_idx]

        rep1_df$rep_value <- as.numeric(NA)
        rep2_df$rep_value <- as.numeric(NA)

        rep1_df$rep_value[top_pairs_df$rep1_idx] <-
            rep2_df$value[top_pairs_df$rep2_idx]
        rep2_df$rep_value[top_pairs_df$rep2_idx] <-
            rep1_df$value[top_pairs_df$rep1_idx]
    } else {
        rep1_df$idx <- integer(0)
        rep2_df$idx <- integer(0)
        rep1_df$rep_idx <- integer(0)
        rep2_df$rep_idx <- integer(0)
        rep1_df$rank <- integer(0)
        rep2_df$rank <- integer(0)
        rep1_df$rep_rank <- integer(0)
        rep2_df$rep_rank <- integer(0)
        rep1_df$rep_value <- numeric(0)
        rep2_df$rep_value <- numeric(0)
    }

    columns <- c(columns, "rep_value", "rank", "rep_rank", "idx", "rep_idx")
    if (analysis_type == "IDR1D") {
        rep1_df <- dplyr::select(rep1_df, columns)
        rep2_df <- dplyr::select(rep2_df, columns)
    } else if (analysis_type == "IDR2D") {
        rep1_df <- dplyr::select(rep1_df, columns)
        rep2_df <- dplyr::select(rep2_df, columns)
    }

    return(list(rep1_df = rep1_df, rep2_df = rep2_df))
}

#' @title Prepares Data for IDR Analysis
#'
#' @description
#' This method removes invalid values, establishes the correct
#' ranking, and breaks ties prior to IDR analysis.
#'
#' \code{Inf} and \code{-Inf} are replaced by \code{max(x) * max_factor}
#' and \code{min(x) / max_factor}, respectively.
#'
#' \code{NA} values in \code{x} are replaced by \code{mean(x)}.
#'
#' All values in \code{x} are transformed using the transformation specified
#' in \code{value_transformation}.
#'
#' Lastly, a small amount of noise is added to \code{x} to break ties. The
#' magnitude of the noise is controlled by \code{jitter_factor}.
#'
#' @param x numeric vector of values
#' @param value_transformation the values in \code{x} have to be transformed in
#' a way such that when ordered in descending order, more significant
#' interactions end up on top of the list. If the values in \code{x} are
#' p-values, \code{"log_additive_inverse"} is recommended. The following
#' transformations are supported:
#' \tabular{rl}{
#'   \code{"identity"} \tab no transformation is performed on \code{x}\cr
#'   \code{"additive_inverse"} \tab \code{x. = -x}\cr
#'   \code{"multiplicative_inverse"} \tab \code{x. = 1 / x}\cr
#'   \code{"log"} \tab \code{x. = log(x)}. Note: zeros are replaced by
#'   \code{.Machine$double.xmin}\cr
#'   \code{"log_additive_inverse"} \tab \code{x. = -log(x)}, recommended if
#'   \code{x} are p-values. Note: zeros are replaced by
#'   \code{.Machine$double.xmin}
#' }
#'
#' either \code{"ascending"} (more significant
#' interactions have lower value in \code{value} column) or \code{"descending"}
#' (more significant interactions have higher value in \code{value} column)
#' @param max_factor numeric; controls the replacement values for \code{Inf}
#' and \code{-Inf}. \code{Inf} are replaced by \code{max(x) * max_factor} and
#' \code{-Inf} are replaced by \code{min(x) / max_factor}.
#' @param jitter_factor numeric; controls the magnitude of the noise that
#' is added to \code{x}. This is done to break ties in \code{x}.
#'
#' @return numeric vector; transformed and stripped values of \code{x}, ready
#' for IDR analysis
#'
#' @examples
#' rep1_df <- idr2d:::chiapet$rep1_df
#' rep1_df$fdr <- preprocess(rep1_df$fdr, "log_additive_inverse")
#'
#' @importFrom futile.logger flog.warn
#' @importFrom stringr str_trim
#' @export
preprocess <- function(x, value_transformation = c("identity",
                                                   "additive_inverse",
                                                   "multiplicative_inverse",
                                                   "log",
                                                   "log_additive_inverse"),
                       max_factor = 1.5,
                       jitter_factor = 0.0001) {
    # argument handling
    value_transformation <- match.arg(value_transformation,
                                      choices = c("identity",
                                                  "additive_inverse",
                                                  "multiplicative_inverse",
                                                  "log",
                                                  "log_additive_inverse"))

    # replace Inf with real values
    max_value <- max(x[x != Inf], na.rm = TRUE)
    x[x == Inf] <- max_value * max_factor

    # replace -Inf with real values
    min_value <- min(x[x != -Inf], na.rm = TRUE)
    x[x == -Inf] <- min_value / max_factor

    # replace NA with mean
    x[is.na(x)] <- mean(x, na.rm = TRUE)

    if (value_transformation == "additive_inverse") {
        x <- (-1) * x
    } else if (value_transformation == "multiplicative_inverse") {
        x <- 1 / x
    } else if (value_transformation == "log") {
        x[x <= .Machine$double.xmin] <- .Machine$double.xmin
        x <- log(x)
    } else if (value_transformation == "log_additive_inverse") {
        x[x <= .Machine$double.xmin] <- .Machine$double.xmin
        x <- (-1) * log(x)
    }

    if (length(unique(x)) < 10) {
        futile.logger::flog.warn("low complexity value distribution")
    }

    # add jitter to break ties
    invisible(tryCatch({
        x <- jitter(x, factor = jitter_factor)
    },
    error = function(e) {
        futile.logger::flog.warn(stringr::str_trim(e))
    }))

    return(x)
}

#' @title Estimates IDR for Genomic Peak Data
#'
#' @description
#' This method estimates Irreproducible Discovery Rates (IDR) for peaks in
#' replicated ChIP-seq experiments.
#'
#' @references
#' Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011) Measuring
#' reproducibility of high-throughput experiments. Annals of Applied
#' Statistics, Vol. 5, No. 3, 1752-1779.
#'
#' @param remove_nonstandard_chromosomes removes peaks containing
#' genomic locations on non-standard chromosomes using
#' \code{\link[GenomeInfoDb:seqlevels-wrappers]{keepStandardChromosomes}}
#' (default is TRUE)
#' @param max_iteration integer; maximum number of iterations for
#' IDR estimation (defaults to 30)
#' @param local_idr see \code{\link[idr:est.IDR]{est.IDR}}
#'
#' @inheritParams establish_bijection1d
#' @inheritParams idr::est.IDR
#' @inheritParams preprocess
#'
#' @return List with three components, (\code{rep1_df}, \code{rep2_df},
#' and \code{analysis_type}) containing the interactions from input
#' data frames \code{rep1_df} and \code{rep2_df} with
#' the following additional columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chr} \tab character; genomic location of peak -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2: \tab \code{start} \tab integer; genomic location of peak  -
#'   start coordinate\cr
#'   column 3: \tab \code{end} \tab integer; genomic location of peak -
#'   end coordinate\cr
#'   column 4: \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the peaks\cr
#'   column 5: \tab \code{rep_value} \tab numeric; value of corresponding
#'   replicate peak. If no corresponding peak was found, \code{rep_value} is set
#'   to \code{NA}.\cr
#'   column 6: \tab \code{rank} \tab integer; rank of the peak, established by
#'   value column, ascending order\cr
#'   column 7: \tab \code{rep_rank} \tab integer; rank of corresponding
#'   replicate peak. If no corresponding peak was found, \code{rep_rank} is
#'   set to \code{NA}.\cr
#'   column 8: \tab \code{idx} \tab integer; peak index, primary key\cr
#'   column 9: \tab \code{rep_idx} \tab integer; specifies the index of the
#'   corresponding peak in the other replicate (foreign key). If no
#'   corresponding peak was found, \code{rep_idx} is set to \code{NA}.\cr
#'   column 10: \tab \code{idr} \tab IDR of the peak and the
#'   corresponding peak in the other replicate. If no corresponding
#'   peak was found, \code{idr} is set to \code{NA}.
#' }
#'
#' @examples
#' idr_results <- estimate_idr1d(idr2d:::chipseq$rep1_df,
#'                               idr2d:::chipseq$rep2_df,
#'                               value_transformation = "log")
#' summary(idr_results)
#'
#' @export
estimate_idr1d <- function(rep1_df, rep2_df,
                           value_transformation = c("identity",
                                                    "additive_inverse",
                                                    "multiplicative_inverse",
                                                    "log",
                                                    "log_additive_inverse"),
                           ambiguity_resolution_method = c("overlap",
                                                           "midpoint",
                                                           "value"),
                           remove_nonstandard_chromosomes = TRUE,
                           max_factor = 1.5, jitter_factor = 0.0001,
                           max_gap = 0L,
                           mu = 0.1, sigma = 1.0, rho = 0.2, p = 0.5,
                           eps = 0.001, max_iteration = 30, local_idr = TRUE) {
    return(estimate_idr(rep1_df, rep2_df, "IDR1D",
                        value_transformation,
                        ambiguity_resolution_method,
                        remove_nonstandard_chromosomes,
                        max_factor, jitter_factor, max_gap,
                        mu, sigma, rho, p,
                        eps, max_iteration, local_idr))
}

#' @title Estimates IDR for Genomic Interaction Data
#'
#' @description
#' This method estimates Irreproducible Discovery Rates (IDR) between
#' two replicates of experiments identifying genomic interactions, such as
#' HiC, ChIA-PET, and HiChIP.
#'
#' @references
#' Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011) Measuring
#' reproducibility of high-throughput experiments. Annals of Applied
#' Statistics, Vol. 5, No. 3, 1752-1779.
#'
#' @param remove_nonstandard_chromosomes removes interactions
#' containing
#' genomic locations on non-standard chromosomes using
#' \code{\link[GenomeInfoDb:seqlevels-wrappers]{keepStandardChromosomes}}
#' (default is TRUE)
#' @param max_iteration integer; maximum number of iterations for
#' IDR estimation (defaults to 30)
#' @param local_idr see \code{\link[idr:est.IDR]{est.IDR}}
#'
#' @inheritParams establish_bijection2d
#' @inheritParams idr::est.IDR
#' @inheritParams preprocess
#'
#' @return List with three components, (\code{rep1_df}, \code{rep2_df},
#' and \code{analysis_type}) containing the interactions from input
#' data frames \code{rep1_df} and \code{rep2_df} with
#' the following additional columns:
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
#'   column 7: \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the interactions\cr
#'   column 8: \tab \code{"rep_value"} \tab numeric; value of corresponding
#'   replicate interaction. If no corresponding interaction was found,
#'   \code{rep_value} is set to \code{NA}.\cr
#'   column 9: \tab \code{"rank"} \tab integer; rank of the interaction,
#'   established by value column, ascending order\cr
#'   column 10: \tab \code{"rep_rank"} \tab integer; rank of corresponding
#'   replicate interaction. If no corresponding interaction was found,
#'   \code{rep_rank} is set to \code{NA}.\cr
#'   column 11: \tab \code{"idx"} \tab integer; interaction index,
#'   primary key\cr
#'   column 12: \tab \code{"rep_idx"} \tab integer; specifies the index of the
#'   corresponding interaction in the other replicate (foreign key). If no
#'   corresponding interaction was found, \code{rep_idx} is set to \code{NA}.\cr
#'   \code{idr} \tab IDR of the interaction and the
#'   corresponding interaction in the other replicate. If no corresponding
#'   interaction was found, \code{idr} is set to \code{NA}.
#' }
#'
#' @examples
#' idr_results <- estimate_idr2d(idr2d:::chiapet$rep1_df,
#'                               idr2d:::chiapet$rep2_df,
#'                               value_transformation = "log_additive_inverse")
#' summary(idr_results)
#'
#' @export
estimate_idr2d <- function(rep1_df, rep2_df,
                           value_transformation = c("identity",
                                                    "additive_inverse",
                                                    "multiplicative_inverse",
                                                    "log",
                                                    "log_additive_inverse"),
                           ambiguity_resolution_method = c("overlap",
                                                           "midpoint",
                                                           "value"),
                           remove_nonstandard_chromosomes = TRUE,
                           max_factor = 1.5, jitter_factor = 0.0001,
                           max_gap = 0L,
                           mu = 0.1, sigma = 1.0, rho = 0.2, p = 0.5,
                           eps = 0.001, max_iteration = 30, local_idr = TRUE) {
    return(estimate_idr(rep1_df, rep2_df, "IDR2D",
                        value_transformation,
                        ambiguity_resolution_method,
                        remove_nonstandard_chromosomes,
                        max_factor, jitter_factor,
                        max_gap,
                        mu, sigma, rho, p,
                        eps, max_iteration, local_idr))
}

#' @title Estimates IDR for Genomic Peaks or Genomic Interactions
#' @references
#' Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011) Measuring
#' reproducibility of high-throughput experiments. Annals of Applied
#' Statistics, Vol. 5, No. 3, 1752-1779.
#'
#' @param remove_nonstandard_chromosomes removes peaks and interactions
#' containing
#' genomic locations on non-standard chromosomes using
#' \code{\link[GenomeInfoDb:seqlevels-wrappers]{keepStandardChromosomes}}
#' (default is TRUE)
#' @param max_iteration integer; maximum number of iterations for
#' IDR estimation (defaults to 30)
#' @param local_idr see \code{\link[idr:est.IDR]{est.IDR}}
#'
#' @inheritParams idr::est.IDR
#' @inheritParams preprocess
#' @inheritParams establish_bijection
#'
#'
#' @return See \code{\link{estimate_idr1d}} or
#' \code{\link{estimate_idr2d}}, respectively.
#'
#' @examples
#' idr_results <- estimate_idr(idr2d:::chiapet$rep1_df,
#'                             idr2d:::chiapet$rep2_df,
#'                             analysis_type = "IDR2D",
#'                             value_transformation = "log_additive_inverse")
#' summary(idr_results)
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom futile.logger flog.warn
#' @importFrom stringr str_trim
#' @importFrom idr est.IDR
#' @export
estimate_idr <- function(rep1_df, rep2_df, analysis_type = "IDR2D",
                         value_transformation = c("identity",
                                                  "additive_inverse",
                                                  "multiplicative_inverse",
                                                  "log",
                                                  "log_additive_inverse"),
                         ambiguity_resolution_method = c("overlap",
                                                         "midpoint",
                                                         "value"),
                         remove_nonstandard_chromosomes = TRUE,
                         max_factor = 1.5, jitter_factor = 0.0001,
                         max_gap = 0L,
                         mu = 0.1, sigma = 1.0, rho = 0.2, p = 0.5,
                         eps = 0.001, max_iteration = 30, local_idr = TRUE) {
    # avoid CRAN warnings
    rep2_idx <- rep1_value <- rep2_value <- idr <- NULL

    # argument handling
    analysis_type <- match.arg(analysis_type, choices = c("IDR1D", "IDR2D"))

    if (analysis_type == "IDR1D") {
        columns <- c("chr", "start", "end", "value")
    } else if (analysis_type == "IDR2D") {
        columns <- c("chr_a", "start_a", "end_a",
                     "chr_b", "start_b", "end_b", "value")
    } else {
        stop("unknown analysis type")
    }

    # remove irrelevant columns, rename relevant columns
    rep1_df <- rep1_df[, seq_len(length(columns))]
    rep2_df <- rep2_df[, seq_len(length(columns))]
    colnames(rep1_df) <- columns
    colnames(rep2_df) <- columns

    if (remove_nonstandard_chromosomes) {
        if (analysis_type == "IDR1D") {
            rep1_df <- remove_nonstandard_chromosomes1d(rep1_df)
            rep2_df <- remove_nonstandard_chromosomes1d(rep2_df)
        } else if (analysis_type == "IDR2D") {
            rep1_df <- remove_nonstandard_chromosomes2d(rep1_df)
            rep2_df <- remove_nonstandard_chromosomes2d(rep2_df)
        }
    }

    rep1_df$value <- preprocess(rep1_df$value,
                                value_transformation,
                                max_factor = max_factor,
                                jitter_factor = jitter_factor)

    rep2_df$value <- preprocess(rep2_df$value,
                                value_transformation,
                                max_factor = max_factor,
                                jitter_factor = jitter_factor)

    mapping <- establish_bijection(rep1_df, rep2_df, analysis_type,
                                   ambiguity_resolution_method,
                                   max_gap = max_gap)

    if (nrow(mapping$rep1_df) > 0 && nrow(mapping$rep2_df) > 0) {
        rep1_df <- mapping$rep1_df
        rep2_df <- mapping$rep2_df
        idx_df <- data.frame(
            rep1_idx = rep1_df$idx,
            rep2_idx = rep1_df$rep_idx,
            rep1_value = rep1_df$value,
            rep2_value = rep1_df$rep_value
        )

        idx_df <- dplyr::filter(idx_df,
                                !is.na(rep2_idx) & !is.infinite(rep2_idx))

        if (nrow(idx_df) > 0) {
            idr_matrix <- as.matrix(dplyr::select(idx_df,
                                                  rep1_value,
                                                  rep2_value))

            if (length(unique(idr_matrix[, 1])) < 10 ||
                length(unique(idr_matrix[, 2])) < 10) {
                futile.logger::flog.warn("low complexity data set")
            }

            invisible(tryCatch({
                idr_results <- idr::est.IDR(idr_matrix, mu, sigma, rho, p,
                                            eps = eps,
                                            max.ite = max_iteration)
                if (local_idr) {
                    idx_df$idr <- idr_results$idr
                } else {
                    idx_df$idr <- idr_results$IDR
                }
            },
            error = function(e) {
                idx_df$idr <- as.numeric(NA)
                futile.logger::flog.warn(stringr::str_trim(e))
            }))
        } else {
            idx_df$idr <- numeric(0)
        }

        if (nrow(rep1_df) > 0) {
            rep1_df$idr <- as.numeric(NA)
            if (nrow(idx_df) > 0) {
                rep1_df$idr[idx_df$rep1_idx] <- idx_df$idr
                rep1_df <- dplyr::arrange(rep1_df, idr)
            }
        } else {
            rep1_df$idr <- numeric(0)
        }

        if (nrow(rep2_df) > 0) {
            rep2_df$idr <- as.numeric(NA)
            if (length(idx_df$idr) > 0) {
                rep2_df$idr[idx_df$rep2_idx] <- idx_df$idr
                rep2_df <- dplyr::arrange(rep2_df, idr)
            }
        } else {
            rep2_df$idr <- numeric(0)
        }
    } else {
        if (nrow(rep1_df) > 0) {
            rep1_df$idr <- as.numeric(NA)
        } else {
            rep1_df$idr <- numeric(0)
        }

        if (nrow(rep2_df) > 0) {
            rep2_df$idr <- as.numeric(NA)
        } else {
            rep2_df$idr <- numeric(0)
        }
    }

    columns <- c(columns, "rep_value", "rank", "rep_rank",
                 "idx", "rep_idx", "idr")
    if (analysis_type == "IDR1D") {
        rep1_df <- dplyr::select(rep1_df, columns)
        rep2_df <- dplyr::select(rep2_df, columns)
    } else if (analysis_type == "IDR2D") {
        rep1_df <- dplyr::select(rep1_df, columns)
        rep2_df <- dplyr::select(rep2_df, columns)
    }

    res <- list(rep1_df = rep1_df, rep2_df = rep2_df,
                analysis_type = analysis_type)
    class(res) <- c("idr2d_result", class(res))
    return(res)
}

#' @importFrom methods is
#' @importFrom dplyr filter
#' @export
summary.idr2d_result <- function(object, ...) {
    idr <- NULL

    stopifnot(methods::is(object, "idr2d_result"))

    num_rep_df <- dplyr::filter(object$rep1_df, !is.na(idr))

    num_sig_df <- dplyr::filter(num_rep_df, idr < 0.05)
    num_high_sig_df <- dplyr::filter(num_rep_df, idr < 0.01)

    res <- list(
        analysis_type = object$analysis_type,
        rep1_num_interactions = nrow(object$rep1_df),
        rep2_num_interactions = nrow(object$rep2_df),
        num_reproducible_interactions = nrow(num_rep_df),
        num_significant_interactions = nrow(num_sig_df),
        num_highly_significant_interactions = nrow(num_high_sig_df)
    )

    class(res) <- c("idr2d_result_summary", class(res))
    return(res)
}

#' @importFrom methods is
#' @export
print.idr2d_result_summary <- function(x, ...) {
    stopifnot(methods::is(x, "idr2d_result_summary"))

    cat("analysis type: ", x$analysis_type, "\n", sep = "")

    if (x$analysis_type == "IDR2D HiC") {
        cat("number of blocks: ",
            x$num_blocks, "\n", sep = "")
        cat("number of blocks with significant IDR (IDR < 0.05): ",
            x$num_significant_blocks, "\n", sep = "")
        cat("number of blocks with highly significant IDR (IDR < 0.01): ",
            x$num_highly_significant_blocks, "\n", sep = "")
        cat("percentage of blocks with significant IDR (IDR < 0.05): ",
            round(x$num_significant_blocks / x$num_blocks * 100, digits = 2),
            " %\n", sep = "")
    } else {
        cat("number of interactions in replicate 1: ",
            x$rep1_num_interactions, "\n", sep = "")
        cat("number of interactions in replicate 2: ",
            x$rep2_num_interactions, "\n", sep = "")
        cat("number of reproducible interactions: ",
            x$num_reproducible_interactions, "\n", sep = "")
        cat("number of interactions with significant IDR (IDR < 0.05): ",
            x$num_significant_interactions, "\n", sep = "")
        cat("number of interactions with highly significant IDR (IDR < 0.01): ",
            x$num_highly_significant_interactions, "\n", sep = "")
        cat("percentage of interactions with significant IDR (IDR < 0.05): ",
            round(x$num_significant_interactions /
                max(x$rep1_num_interactions, x$rep2_num_interactions) * 100,
                digits = 2),
            " %\n", sep = "")
    }
}
