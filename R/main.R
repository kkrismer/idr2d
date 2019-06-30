
#' @title Finds One-to-One Correspondence between Peaks from Replicate 1
#' and 2
#'
#' @description
#' This method establishes a bijective assignment between peaks from
#' replicate 1 and 2. A peak in replicate 1 is assigned to a
#' peak in replicate 2 if and only if (1) they overlap (or the gap between the
#' peaks is less than or equal to \code{max.gap}), and (2) there is no other
#' peak in
#' replicate 2 that overlaps with the peak in replicate 1 and has a
#' lower \emph{ambiguity resolution value}.
#'
#' @inheritParams overlap1D
#'
#' @return Data frames \code{rep1.df} and \code{rep2.df} with
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
#'   column 5: \tab \code{rep.value} \tab numeric; value of corresponding
#'   replicate peak If no corresponding peak was found, \code{rep.value} is set
#'   to \code{NA}.\cr
#'   column 6: \tab \code{rank} \tab integer; rank of the peak, established by
#'   value column, ascending order\cr
#'   column 7: \tab \code{rep.rank} \tab integer; rank of corresponding
#'   replicate peak. If no corresponding peak was found, \code{rep.rank} is
#'   set to \code{NA}.\cr
#'   column 8: \tab \code{idx} \tab integer; peak index, primary key\cr
#'   column 9: \tab \code{rep.idx} \tab integer; specifies the index of the
#'   corresponding peak in the other replicate (foreign key). If no
#'   corresponding peak was found, \code{rep.idx} is set to \code{NA}.
#' }
#'
#' @examples
#' rep1.df <- idr2d:::chipseq$rep1.df
#' rep1.df$fdr <- preprocess(rep1.df$fdr, "log.additive.inverse")
#'
#' rep2.df <- idr2d:::chipseq$rep2.df
#' rep2.df$fdr <- preprocess(rep2.df$fdr, "log.additive.inverse")
#'
#' mapping <- establishBijection1D(rep1.df, rep2.df)
#'
#' @export
establishBijection1D <- function(rep1.df, rep2.df,
                                 ambiguity.resolution.method = c("overlap",
                                                                 "midpoint",
                                                                 "value"),
                                 max.gap = 0L) {
    return(establishBijection(rep1.df, rep2.df, "IDR1D",
                              ambiguity.resolution.method))
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
#' or equal to \code{max.gap}), and (2) there is no other interaction in
#' replicate 2 that overlaps with the interaction in replicate 1 and has a
#' lower \emph{ambiguity resolution value}.
#'
#' @inheritParams overlap2D
#'
#' @return Data frames \code{rep1.df} and \code{rep2.df} with
#' the following columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chr.a} \tab character; genomic location of anchor A -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2: \tab \code{start.a} \tab integer; genomic location of anchor A -
#'   start coordinate\cr
#'   column 3: \tab \code{end.a} \tab integer; genomic location of anchor A -
#'   end coordinate\cr
#'   column 4: \tab \code{chr.b} \tab character; genomic location of anchor B -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 5: \tab \code{start.b} \tab integer; genomic location of anchor B -
#'   start coordinate\cr
#'   column 6: \tab \code{end.b} \tab integer; genomic location of anchor B -
#'   end coordinate\cr
#'   column 7: \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the interactions\cr
#'   column 8: \tab \code{"rep.value"} \tab numeric; value of corresponding
#'   replicate interaction. If no corresponding interaction was found,
#'   \code{rep.value} is set to \code{NA}.\cr
#'   column 9: \tab \code{"rank"} \tab integer; rank of the interaction,
#'   established by value column, ascending order\cr
#'   column 10: \tab \code{"rep.rank"} \tab integer; rank of corresponding
#'   replicate interaction. If no corresponding interaction was found,
#'   \code{rep.rank} is set to \code{NA}.\cr
#'   column 11: \tab \code{"idx"} \tab integer; interaction index,
#'   primary key\cr
#'   column 12: \tab \code{"rep.idx"} \tab integer; specifies the index of the
#'   corresponding interaction in the other replicate (foreign key). If no
#'   corresponding interaction was found, \code{rep.idx} is set to \code{NA}.
#' }
#'
#' @examples
#' rep1.df <- idr2d:::chiapet$rep1.df
#' rep1.df$fdr <- preprocess(rep1.df$fdr, "log.additive.inverse")
#'
#' rep2.df <- idr2d:::chiapet$rep2.df
#' rep2.df$fdr <- preprocess(rep2.df$fdr, "log.additive.inverse")
#'
#' mapping <- establishBijection2D(rep1.df, rep2.df)
#'
#' @export
establishBijection2D <- function(rep1.df, rep2.df,
                                 ambiguity.resolution.method = c("overlap",
                                                                 "midpoint",
                                                                 "value"),
                                 max.gap = 0L) {
    return(establishBijection(rep1.df, rep2.df, "IDR2D",
                              ambiguity.resolution.method))
}

#' @title Finds One-to-One Correspondence between Peaks or interactions
#' from Replicate 1 and 2
#'
#' @param analysis.type "IDR2D" for genomic interaction data sets,
#' "IDR1D" for genomic peak data sets
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' @importFrom dplyr select
#' @importFrom futile.logger flog.warn
#' @importFrom magrittr "%>%"
#' @export
establishBijection <- function(rep1.df, rep2.df,
                               analysis.type = c("IDR1D", "IDR2D"),
                               ambiguity.resolution.method = c("overlap",
                                                               "midpoint",
                                                               "value"),
                               max.gap = 0L) {
    # avoid CRAN warnings
    rep1.idx <- rep2.idx <- arv <- NULL
    chr <- start <- end <- NULL
    chr.a <- start.a <- end.a <- chr.b <- start.b <- end.b <- NULL
    value <- rep.value <- rank <- rep.rank <- idx <- rep.idx <- NULL

    # argument handling
    analysis.type <- match.arg(analysis.type, choices = c("IDR1D", "IDR2D"))

    if (analysis.type == "IDR1D") {
        columns <- c("chr", "start", "end", "value")
    } else if (analysis.type == "IDR2D") {
        columns <- c("chr.a", "start.a", "end.a",
                     "chr.b", "start.b", "end.b", "value")
    } else {
        stop("unknown analysis type")
    }

    # remove irrelevant columns, rename relevant columns
    rep1.df <- rep1.df[, seq_len(length(columns))]
    rep2.df <- rep2.df[, seq_len(length(columns))]
    colnames(rep1.df) <- columns
    colnames(rep2.df) <- columns

    if (nrow(rep1.df) > 0 && nrow(rep2.df) > 0) {
        # shuffle to break preexisting order
        rep1.df <- rep1.df[sample.int(nrow(rep1.df)), ]
        rep2.df <- rep2.df[sample.int(nrow(rep2.df)), ]

        # sort by value column, descending order (higher values are better)
        rep1.df <- dplyr::arrange(rep1.df, dplyr::desc(value))
        rep2.df <- dplyr::arrange(rep2.df, dplyr::desc(value))

        # add idx column
        rep1.df$idx <- seq_len(nrow(rep1.df))
        rep2.df$idx <- seq_len(nrow(rep2.df))

        # add rank column
        rep1.df$rank <- seq_len(nrow(rep1.df))
        rep2.df$rank <- seq_len(nrow(rep2.df))

        if (analysis.type == "IDR1D") {
            pairs.df <- overlap1D(rep1.df, rep2.df,
                                  ambiguity.resolution.method, max.gap)
        } else if (analysis.type == "IDR2D") {
            pairs.df <- overlap2D(rep1.df, rep2.df,
                                  ambiguity.resolution.method, max.gap)
        }

        top.pairs.df <- pairs.df %>% dplyr::group_by(rep1.idx) %>%
            dplyr::slice(which.min(arv))
        top.pairs.df <- top.pairs.df %>% dplyr::group_by(rep2.idx) %>%
            dplyr::slice(which.min(arv))

        if (nrow(top.pairs.df) != length(unique(top.pairs.df$rep1.idx)) &&
            nrow(top.pairs.df) != length(unique(top.pairs.df$rep2.idx))) {
            futile.logger::flog.warn("ambiguous interaction assignments")
        }

        rep1.df$rep.idx <- as.integer(NA)
        rep2.df$rep.idx <- as.integer(NA)

        rep1.df$rep.idx[top.pairs.df$rep1.idx] <- top.pairs.df$rep2.idx
        rep2.df$rep.idx[top.pairs.df$rep2.idx] <- top.pairs.df$rep1.idx

        rep1.df$rep.rank <- as.integer(NA)
        rep2.df$rep.rank <- as.integer(NA)

        rep1.df$rep.rank[top.pairs.df$rep1.idx] <-
            rep2.df$rank[top.pairs.df$rep2.idx]
        rep2.df$rep.rank[top.pairs.df$rep2.idx] <-
            rep1.df$rank[top.pairs.df$rep1.idx]

        rep1.df$rep.value <- as.numeric(NA)
        rep2.df$rep.value <- as.numeric(NA)

        rep1.df$rep.value[top.pairs.df$rep1.idx] <-
            rep2.df$value[top.pairs.df$rep2.idx]
        rep2.df$rep.value[top.pairs.df$rep2.idx] <-
            rep1.df$value[top.pairs.df$rep1.idx]
    } else {
        rep1.df$idx <- integer(0)
        rep2.df$idx <- integer(0)
        rep1.df$rep.idx <- integer(0)
        rep2.df$rep.idx <- integer(0)
        rep1.df$rank <- integer(0)
        rep2.df$rank <- integer(0)
        rep1.df$rep.rank <- integer(0)
        rep2.df$rep.rank <- integer(0)
        rep1.df$rep.value <- numeric(0)
        rep2.df$rep.value <- numeric(0)
    }

    if (analysis.type == "IDR1D") {
        rep1.df <- dplyr::select(rep1.df, chr, start, end,
                                 value, rep.value,
                                 rank, rep.rank,
                                 idx, rep.idx)
        rep2.df <- dplyr::select(rep2.df, chr, start, end,
                                 value, rep.value,
                                 rank, rep.rank,
                                 idx, rep.idx)
    } else if (analysis.type == "IDR2D") {
        rep1.df <- dplyr::select(rep1.df, chr.a, start.a, end.a,
                                 chr.b, start.b, end.b,
                                 value, rep.value,
                                 rank, rep.rank,
                                 idx, rep.idx)
        rep2.df <- dplyr::select(rep2.df, chr.a, start.a, end.a,
                                 chr.b, start.b, end.b,
                                 value, rep.value,
                                 rank, rep.rank,
                                 idx, rep.idx)
    }

    return(list(rep1.df = rep1.df, rep2.df = rep2.df))
}

#' @title Prepares Data for IDR Analysis
#'
#' @description
#' This method removes invalid values, establishes the correct
#' ranking, and breaks ties prior to IDR analysis.
#'
#' \code{Inf} and \code{-Inf} are replaced by \code{max(x) * max.factor}
#' and \code{min(x) / max.factor}, respectively.
#'
#' \code{NA} values in \code{x} are replaced by \code{mean(x)}.
#'
#' All values in \code{x} are transformed using the transformation specified
#' in \code{value.transformation}.
#'
#' Lastly, a small amount of noise is added to \code{x} to break ties. The
#' magnitude of the noise is controlled by \code{jitter.factor}.
#'
#' @param x numeric vector of values
#' @param value.transformation the values in \code{x} have to be transformed in
#' a way such that when ordered in descending order, more significant
#' interactions end up on top of the list. If the values in \code{x} are
#' p-values, \code{"log.additive.inverse"} is recommended. The following
#' transformations are supported:
#' \tabular{rl}{
#'   \code{"identity"} \tab no transformation is performed on \code{x}\cr
#'   \code{"additive.inverse"} \tab \code{x = -x}\cr
#'   \code{"multiplicative.inverse"} \tab \code{x = 1 / x}\cr
#'   \code{"log"} \tab \code{x = log(x)}. Note: zeros are replaced by
#'   \code{.Machine$double.xmin}\cr
#'   \code{"log.additive.inverse"} \tab \code{x = -log(x)}, recommended if
#'   \code{x} are p-values. Note: zeros are replaced by
#'   \code{.Machine$double.xmin}
#' }
#'
#' either \code{"ascending"} (more significant
#' interactions have lower value in \code{value} column) or \code{"descending"}
#' (more significant interactions have higher value in \code{value} column)
#' @param max.factor numeric; controls the replacement values for \code{Inf}
#' and \code{-Inf}. \code{Inf} are replaced by \code{max(x) * max.factor} and
#' \code{-Inf} are replaced by \code{min(x) / max.factor}.
#' @param jitter.factor numeric; controls the magnitude of the noise that
#' is added to \code{x}. This is done to break ties in \code{x}.
#'
#' @return numeric vector; transformed and stripped values of \code{x}, ready
#' for IDR analysis
#'
#' @examples
#' rep1.df <- idr2d:::chiapet$rep1.df
#' rep1.df$fdr <- preprocess(rep1.df$fdr, "log.additive.inverse")
#'
#' @importFrom futile.logger flog.warn
#' @importFrom stringr str_trim
#' @export
preprocess <- function(x, value.transformation = c("identity",
                                                   "additive.inverse",
                                                   "multiplicative.inverse",
                                                   "log",
                                                   "log.additive.inverse"),
                       max.factor = 1.5,
                       jitter.factor = 0.0001) {
    # argument handling
    value.transformation <- match.arg(value.transformation,
                                      choices = c("identity",
                                                  "additive.inverse",
                                                  "multiplicative.inverse",
                                                  "log",
                                                  "log.additive.inverse"))

    # replace Inf with real values
    max.value <- max(x[x != Inf], na.rm = TRUE)
    x[x == Inf] <- max.value * max.factor

    # replace -Inf with real values
    min.value <- min(x[x != -Inf], na.rm = TRUE)
    x[x == -Inf] <- min.value / max.factor

    # replace NA with mean
    x[is.na(x)] <- mean(x, na.rm = TRUE)

    if (value.transformation == "additive.inverse") {
        x <- (-1) * x
    } else if (value.transformation == "multiplicative.inverse") {
        x <- 1 / x
    } else if (value.transformation == "log") {
        x[x <= .Machine$double.xmin] <- .Machine$double.xmin
        x <- log(x)
    } else if (value.transformation == "log.additive.inverse") {
        x[x <= .Machine$double.xmin] <- .Machine$double.xmin
        x <- (-1) * log(x)
    }

    if (length(unique(x)) < 10) {
        futile.logger::flog.warn("low complexity value distribution")
    }

    # add jitter to break ties
    invisible(tryCatch({
        x <- jitter(x, factor = jitter.factor)
    }, error = function(e) {
        futile.logger::flog.warn(stringr::str_trim(e))
    }))

    return(x)
}

#' @title Estimates IDR for Genomic Peak Data
#'
#' @references
#' Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011) Measuring
#' reproducibility of high-throughput experiments. Annals of Applied
#' Statistics, Vol. 5, No. 3, 1752-1779.
#'
#' @inheritParams estimateIDR
#'
#' @return List with two components (\code{rep1.df} and \code{rep1.df})
#' containing the peaks from input data frames \code{rep1.df} and
#' \code{rep2.df} with
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
#'   column 5: \tab \code{rep.value} \tab numeric; value of corresponding
#'   replicate peak If no corresponding peak was found, \code{rep.value} is set
#'   to \code{NA}.\cr
#'   column 6: \tab \code{rank} \tab integer; rank of the peak, established by
#'   value column, ascending order\cr
#'   column 7: \tab \code{rep.rank} \tab integer; rank of corresponding
#'   replicate peak. If no corresponding peak was found, \code{rep.rank} is
#'   set to \code{NA}.\cr
#'   column 8: \tab \code{idx} \tab integer; peak index, primary key\cr
#'   column 9: \tab \code{rep.idx} \tab integer; specifies the index of the
#'   corresponding peak in the other replicate (foreign key). If no
#'   corresponding peak was found, \code{rep.idx} is set to \code{NA}.\cr
#'   column 10: \tab \code{idr} \tab IDR of the peak and the
#'   corresponding peak in the other replicate. If no corresponding
#'   peak was found, \code{idr} is set to \code{NA}.
#' }
#'
#' @examples
#' # TODO
#'
#' @export
estimateIDR1D <- function(rep1.df, rep2.df,
                          value.transformation = c("identity",
                                                   "additive.inverse",
                                                   "multiplicative.inverse",
                                                   "log",
                                                   "log.additive.inverse"),
                          ambiguity.resolution.method = c("overlap",
                                                          "midpoint",
                                                          "value"),
                          remove.nonstandard.chromosomes = TRUE,
                          max.factor = 1.5, jitter.factor = 0.0001,
                          max.gap = 0L,
                          mu = 0.1, sigma = 1.0, rho = 0.2, p = 0.5,
                          eps = 0.001, max.iteration = 30, local.idr = TRUE) {
    return(estimateIDR(rep1.df, rep2.df, "IDR1D",
                       value.transformation,
                       ambiguity.resolution.method,
                       remove.nonstandard.chromosomes,
                       max.factor, jitter.factor, max.gap,
                       mu, sigma, rho, p,
                       eps, max.iteration, local.idr))
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
#' @inheritParams estimateIDR
#'
#' @return List with two components (\code{rep1.df} and \code{rep1.df})
#' containing the interactions from input data frames \code{rep1.df} and
#' \code{rep2.df} with
#' the following additional columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chr.a} \tab character; genomic location of anchor A -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2: \tab \code{start.a} \tab integer; genomic location of anchor A -
#'   start coordinate\cr
#'   column 3: \tab \code{end.a} \tab integer; genomic location of anchor A -
#'   end coordinate\cr
#'   column 4: \tab \code{chr.b} \tab character; genomic location of anchor B -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 5: \tab \code{start.b} \tab integer; genomic location of anchor B -
#'   start coordinate\cr
#'   column 6: \tab \code{end.b} \tab integer; genomic location of anchor B -
#'   end coordinate\cr
#'   column 7: \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the interactions\cr
#'   column 8: \tab \code{"rep.value"} \tab numeric; value of corresponding
#'   replicate interaction. If no corresponding interaction was found,
#'   \code{rep.value} is set to \code{NA}.\cr
#'   column 9: \tab \code{"rank"} \tab integer; rank of the interaction,
#'   established by value column, ascending order\cr
#'   column 10: \tab \code{"rep.rank"} \tab integer; rank of corresponding
#'   replicate interaction. If no corresponding interaction was found,
#'   \code{rep.rank} is set to \code{NA}.\cr
#'   column 11: \tab \code{"idx"} \tab integer; interaction index,
#'   primary key\cr
#'   column 12: \tab \code{"rep.idx"} \tab integer; specifies the index of the
#'   corresponding interaction in the other replicate (foreign key). If no
#'   corresponding interaction was found, \code{rep.idx} is set to \code{NA}.\cr
#'   \code{idr} \tab IDR of the interaction and the
#'   corresponding interaction in the other replicate. If no corresponding
#'   interaction was found, \code{idr} is set to \code{NA}.
#' }
#'
#' @examples
#' idr.results <- estimateIDR2D(idr2d:::chiapet$rep1.df,
#'                            idr2d:::chiapet$rep2.df,
#'                            value.transformation = "log.additive.inverse")
#'
#' @export
estimateIDR2D <- function(rep1.df, rep2.df,
                          value.transformation = c("identity",
                                                   "additive.inverse",
                                                   "multiplicative.inverse",
                                                   "log",
                                                   "log.additive.inverse"),
                          ambiguity.resolution.method = c("overlap",
                                                          "midpoint",
                                                          "value"),
                          remove.nonstandard.chromosomes = TRUE,
                          max.factor = 1.5, jitter.factor = 0.0001,
                          max.gap = 0L,
                          mu = 0.1, sigma = 1.0, rho = 0.2, p = 0.5,
                          eps = 0.001, max.iteration = 30, local.idr = TRUE) {
    return(estimateIDR(rep1.df, rep2.df, "IDR2D",
                       value.transformation,
                       ambiguity.resolution.method,
                       remove.nonstandard.chromosomes,
                       max.factor, jitter.factor,
                       max.gap,
                       mu, sigma, rho, p,
                       eps, max.iteration, local.idr))
}

#' @title Estimates IDR for Genomic Peaks or Genomic Interactions
#' @references
#' Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011) Measuring
#' reproducibility of high-throughput experiments. Annals of Applied
#' Statistics, Vol. 5, No. 3, 1752-1779.
#'
#' @param remove.nonstandard.chromosomes removes peaks and interactions
#' containing
#' genomic locations on non-standard chromosomes using
#' \code{\link[GenomeInfoDb:seqlevels-wrappers]{keepStandardChromosomes}}
#' (default is TRUE)
#' @param max.iteration integer; maximum number of iterations for
#' IDR estimation (defaults to 30)
#' @param local.idr TODO
#'
#' @inheritParams idr::est.IDR
#' @inheritParams preprocess
#' @inheritParams establishBijection
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom futile.logger flog.warn
#' @importFrom stringr str_trim
#' @importFrom idr est.IDR
#' @export
estimateIDR <- function(rep1.df, rep2.df, analysis.type = "IDR2D",
                        value.transformation = c("identity",
                                                 "additive.inverse",
                                                 "multiplicative.inverse",
                                                 "log",
                                                 "log.additive.inverse"),
                        ambiguity.resolution.method = c("overlap",
                                                        "midpoint",
                                                        "value"),
                        remove.nonstandard.chromosomes = TRUE,
                        max.factor = 1.5, jitter.factor = 0.0001,
                        max.gap = 0L,
                        mu = 0.1, sigma = 1.0, rho = 0.2, p = 0.5,
                        eps = 0.001, max.iteration = 30, local.idr = TRUE) {
    # avoid CRAN warnings
    rep2.idx <- rep1.value <- rep2.value <- NULL
    chr <- start <- end <- NULL
    chr.a <- start.a <- end.a <- chr.b <- start.b <- end.b <- NULL
    value <- rep.value <- rank <- rep.rank <- idx <- rep.idx <- idr <- NULL

    # argument handling
    analysis.type <- match.arg(analysis.type, choices = c("IDR1D", "IDR2D"))

    if (analysis.type == "IDR1D") {
        columns <- c("chr", "start", "end", "value")
    } else if (analysis.type == "IDR2D") {
        columns <- c("chr.a", "start.a", "end.a",
                     "chr.b", "start.b", "end.b", "value")
    } else {
        stop("unknown analysis type")
    }

    # remove irrelevant columns, rename relevant columns
    rep1.df <- rep1.df[, seq_len(length(columns))]
    rep2.df <- rep2.df[, seq_len(length(columns))]
    colnames(rep1.df) <- columns
    colnames(rep2.df) <- columns

    if (remove.nonstandard.chromosomes) {
        if (analysis.type == "IDR1D") {
            rep1.df <- removeNonstandardChromosomes1D(rep1.df)
            rep2.df <- removeNonstandardChromosomes1D(rep2.df)
        } else if (analysis.type == "IDR2D") {
            rep1.df <- removeNonstandardChromosomes2D(rep1.df)
            rep2.df <- removeNonstandardChromosomes2D(rep2.df)
        }
    }

    rep1.df$value <- preprocess(rep1.df$value,
                                value.transformation,
                                max.factor = max.factor,
                                jitter.factor = jitter.factor)

    rep2.df$value <- preprocess(rep2.df$value,
                                value.transformation,
                                max.factor = max.factor,
                                jitter.factor = jitter.factor)

    mapping <- establishBijection(rep1.df, rep2.df, analysis.type,
                                  ambiguity.resolution.method,
                                  max.gap = max.gap)

    if (nrow(mapping$rep1.df) > 0 && nrow(mapping$rep2.df) > 0) {
        rep1.df <- mapping$rep1.df
        rep2.df <- mapping$rep2.df
        idx.df <- data.frame(
            rep1.idx = rep1.df$idx,
            rep2.idx = rep1.df$rep.idx,
            rep1.value = rep1.df$value,
            rep2.value = rep1.df$rep.value
        )

        idx.df <- dplyr::filter(idx.df,
                                !is.na(rep2.idx) & !is.infinite(rep2.idx))

        if (nrow(idx.df) > 0) {
            idr.matrix <- as.matrix(dplyr::select(idx.df,
                                                  rep1.value,
                                                  rep2.value))

            if (length(unique(idr.matrix[, 1])) < 10 ||
                length(unique(idr.matrix[, 2])) < 10) {
                futile.logger::flog.warn("low complexity data set")
            }

            invisible(tryCatch({
                idr.results <- idr::est.IDR(idr.matrix, mu, sigma, rho, p,
                                            eps = eps,
                                            max.ite = max.iteration)
                if (local.idr) {
                    idx.df$idr <- idr.results$idr
                } else {
                    idx.df$idr <- idr.results$IDR
                }
            }, error = function(e) {
                idx.df$idr <- as.numeric(NA)
                futile.logger::flog.warn(stringr::str_trim(e))
            }))
        } else {
            idx.df$idr <- numeric(0)
        }

        if (nrow(rep1.df) > 0) {
            rep1.df$idr <- as.numeric(NA)
            if (nrow(idx.df) > 0) {
                rep1.df$idr[idx.df$rep1.idx] <- idx.df$idr
                rep1.df <- dplyr::arrange(rep1.df, idr)
            }
        } else {
            rep1.df$idr <- numeric(0)
        }

        if (nrow(rep2.df) > 0) {
            rep2.df$idr <- as.numeric(NA)
            if (length(idx.df$idr) > 0) {
                rep2.df$idr[idx.df$rep2.idx] <- idx.df$idr
                rep2.df <- dplyr::arrange(rep2.df, idr)
            }
        } else {
            rep2.df$idr <- numeric(0)
        }
    } else {
        if (nrow(rep1.df) > 0) {
            rep1.df$idr <- as.numeric(NA)
        } else {
            rep1.df$idr <- numeric(0)
        }

        if (nrow(rep2.df) > 0) {
            rep2.df$idr <- as.numeric(NA)
        } else {
            rep2.df$idr <- numeric(0)
        }
    }

    if (analysis.type == "IDR1D") {
        rep1.df <- dplyr::select(rep1.df, chr, start, end,
                                 value, rep.value,
                                 rank, rep.rank,
                                 idx, rep.idx, idr)
        rep2.df <- dplyr::select(rep2.df, chr, start, end,
                                 value, rep.value,
                                 rank, rep.rank,
                                 idx, rep.idx, idr)
    } else if (analysis.type == "IDR2D") {
        rep1.df <- dplyr::select(rep1.df, chr.a, start.a, end.a,
                                 chr.b, start.b, end.b,
                                 value, rep.value,
                                 rank, rep.rank,
                                 idx, rep.idx, idr)
        rep2.df <- dplyr::select(rep2.df, chr.a, start.a, end.a,
                                 chr.b, start.b, end.b,
                                 value, rep.value,
                                 rank, rep.rank,
                                 idx, rep.idx, idr)
    }

    return(list(rep1.df = rep1.df, rep2.df = rep2.df))
}
