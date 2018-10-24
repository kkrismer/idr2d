
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
#' @inheritParams overlap
#'
#' @return Original data frames \code{rep1.df} and \code{rep2.df} with
#' additional column \code{replicate.idx}, which specifies the index of the
#' corresponding interaction in the other replicate. If no corresponding
#' interaction was found, \code{replicate.idx} is set to \code{NA}.
#'
#' @examples
#' rep1.df <- idr2d:::chiapet$rep1.df
#' rep1.df$fdr <- preprocess(rep1.df$fdr, "multiplicative.inverse")
#'
#' rep2.df <- idr2d:::chiapet$rep2.df
#' rep2.df$fdr <- preprocess(rep2.df$fdr, "multiplicative.inverse")
#'
#' mapping <- establishBijection(rep1.df, rep2.df)
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' @importFrom futile.logger flog.info
#' @importFrom futile.logger flog.warn
#' @importFrom magrittr "%>%"
#' @export
establishBijection <- function(rep1.df, rep2.df,
                               ambiguity.resolution.method = c("value",
                                                               "overlap",
                                                               "midpoint",
                                                               "expansion"),
                               max.gap = 1000L) {
    # avoid CRAN warnings
    rep1.idx <- rep2.idx <- arv <- NULL

    ambiguity.resolution.method <- match.arg(ambiguity.resolution.method,
                                             choices = c("value",
                                                         "overlap",
                                                         "midpoint",
                                                         "expansion"))

    if (nrow(rep1.df) > 0 && nrow(rep2.df) > 0) {
        # shuffle to break preexisting order
        rep1.df <- rep1.df[sample.int(nrow(rep1.df)), ]
        rep2.df <- rep2.df[sample.int(nrow(rep2.df)), ]

        # sort by value column
        rep1.df <- dplyr::arrange(rep1.df, rep1.df[, 7])
        rep2.df <- dplyr::arrange(rep2.df, rep2.df[, 7])

        pairs.df <- overlap(rep1.df, rep2.df,
                            ambiguity.resolution.method, max.gap)

        if (ambiguity.resolution.method == "expansion") {
            rep1.df <- data.frame(
                chrA = rep1.df$chrA[pairs.df$rep1.idx],
                startA = rep1.df$startA[pairs.df$rep1.idx],
                endA = rep1.df$endA[pairs.df$rep1.idx],
                chrB = rep1.df$chrB[pairs.df$rep1.idx],
                startB = rep1.df$startB[pairs.df$rep1.idx],
                endB = rep1.df$endB[pairs.df$rep1.idx],
                countAB = rep1.df$countAB[pairs.df$rep1.idx],
                countA = rep1.df$countA[pairs.df$rep1.idx],
                countB = rep1.df$countB[pairs.df$rep1.idx],
                lengthA = rep1.df$lengthA[pairs.df$rep1.idx],
                lengthB = rep1.df$lengthB[pairs.df$rep1.idx],
                distance = rep1.df$distance[pairs.df$rep1.idx],
                pp = rep1.df$pp[pairs.df$rep1.idx],
                pvalue = rep1.df$pvalue[pairs.df$rep1.idx],
                fdr = rep1.df$fdr[pairs.df$rep1.idx],
                replicate.idx = seq_len(nrow(pairs.df))
            )
            rep2.df <- data.frame(
                chrA = rep2.df$chrA[pairs.df$rep2.idx],
                startA = rep2.df$startA[pairs.df$rep2.idx],
                endA = rep2.df$endA[pairs.df$rep2.idx],
                chrB = rep2.df$chrB[pairs.df$rep2.idx],
                startB = rep2.df$startB[pairs.df$rep2.idx],
                endB = rep2.df$endB[pairs.df$rep2.idx],
                countAB = rep2.df$countAB[pairs.df$rep2.idx],
                countA = rep2.df$countA[pairs.df$rep2.idx],
                countB = rep2.df$countB[pairs.df$rep2.idx],
                lengthA = rep2.df$lengthA[pairs.df$rep2.idx],
                lengthB = rep2.df$lengthB[pairs.df$rep2.idx],
                distance = rep2.df$distance[pairs.df$rep2.idx],
                pp = rep2.df$pp[pairs.df$rep2.idx],
                pvalue = rep2.df$pvalue[pairs.df$rep2.idx],
                fdr = rep2.df$fdr[pairs.df$rep2.idx],
                replicate.idx = seq_len(nrow(pairs.df))
            )
        } else {
            top.pairs.df <- pairs.df %>% dplyr::group_by(rep1.idx) %>%
                dplyr::slice(which.min(arv))
            top.pairs.df <- top.pairs.df %>% dplyr::group_by(rep2.idx) %>%
                dplyr::slice(which.min(arv))

            futile.logger::flog.info("top pairs selected")

            # replicated interaction ID

            rep1.df$replicate.idx <- NA
            rep2.df$replicate.idx <- NA

            if (nrow(top.pairs.df) != length(unique(top.pairs.df$rep1.idx)) &&
                nrow(top.pairs.df) != length(unique(top.pairs.df$rep2.idx))) {
                futile.logger::flog.warn("ambiguous interaction assignments")
            }

            rep1.df$replicate.idx[top.pairs.df$rep1.idx] <-
                top.pairs.df$rep2.idx
            rep2.df$replicate.idx[top.pairs.df$rep2.idx] <-
                top.pairs.df$rep1.idx
        }

        futile.logger::flog.info("pairs assigned")
    } else {
        rep1.df$replicate.idx <- integer(0)
        rep2.df$replicate.idx <- integer(0)
    }

    futile.logger::flog.info("overlap established")
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
#' rep1.df$fdr <- preprocess(rep1.df$fdr, "multiplicative.inverse")
#'
#' @importFrom idr est.IDR
#' @export
preprocess <- function(x, value.transformation = c("identity",
                                                   "additive.inverse",
                                                   "multiplicative.inverse",
                                                   "log",
                                                   "log.additive.inverse"),
                       max.factor = 1.5,
                       jitter.factor = 0.0001) {
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

    # add jitter to break ties
    x <- jitter(x, factor = jitter.factor)

    return(x)
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
#' @param max.iteration integer; maximum number of iterations for
#' IDR estimation (defaults to 30)
#' @inheritParams idr::est.IDR
#' @inheritParams preprocess
#' @inheritParams establishBijection
#'
#' @return Original data frames \code{rep1.df} and \code{rep2.df} with
#' additional column \code{idr}, which holds the IDR of the interaction and the
#' corresponding interaction in the other replicate. If no corresponding
#' interaction was found, \code{idr} is set to \code{NA}.
#'
#' @examples
#' idr.df <- estimateIDR(idr2d:::chiapet$rep1.df,
#'                       idr2d:::chiapet$rep2.df,
#'                       value.transformation = "multiplicative.inverse")
#'
#' @importFrom futile.logger flog.info
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @export
estimateIDR <- function(rep1.df, rep2.df,
                        value.transformation = c("identity",
                                                 "additive.inverse",
                                                 "multiplicative.inverse",
                                                 "log",
                                                 "log.additive.inverse"),
                        ambiguity.resolution.method = c("value",
                                                        "overlap",
                                                        "midpoint",
                                                        "expansion"),
                        max.factor = 1.5,
                        jitter.factor = 0.0001,
                        mu = 0.1, sigma = 1.0, rho = 0.2, p = 0.5,
                        eps = 0.001, max.iteration = 30) {
    # avoid CRAN warnings
    rep2.idx <- rep1.value <- rep2.value <- idr <- NULL

    rep1.df[, 7] <- preprocess(rep1.df[, 7], value.transformation,
                                max.factor = max.factor,
                                jitter.factor = jitter.factor)

    rep2.df[, 7] <- preprocess(rep2.df[, 7], value.transformation,
                                max.factor = max.factor,
                                jitter.factor = jitter.factor)

    mapping <- establishBijection(rep1.df, rep2.df, ambiguity.resolution.method)

    idx.df <- data.frame(
        rep1.idx = seq_len(nrow(mapping$rep1.df)),
        rep2.idx = mapping$rep1.df$replicate.idx,
        rep1.value = mapping$rep1.df[, 7],
        rep2.value = mapping$rep2.df[mapping$rep1.df$replicate.idx, 7]
    )

    idx.df <- dplyr::filter(idx.df, !is.na(rep2.idx) & !is.infinite(rep2.idx))

    idr.matrix <- as.matrix(dplyr::select(idx.df, rep1.value, rep2.value))

    idr.results <- idr::est.IDR(idr.matrix, mu, sigma, rho, p, eps = eps,
                                max.ite = max.iteration)

    # TODO check IDR vs idr?
    idx.df$idr <- idr.results$IDR

    if (nrow(rep1.df) > 0) {
        rep1.df$idr <- as.numeric(NA)
        futile.logger::flog.info(nrow(rep1.df))
        futile.logger::flog.info(nrow(idx.df))

        rep1.df$idr[idx.df$rep1.idx] <- idx.df$idr
        rep1.df <- dplyr::arrange(rep1.df, idr)
    } else {
        rep1.df$idr <- numeric(0)
    }

    if (nrow(rep2.df) > 0) {
        rep2.df$idr <- as.numeric(NA)
        futile.logger::flog.info(nrow(rep2.df))
        futile.logger::flog.info(nrow(idx.df))

        rep2.df$idr[idx.df$rep2.idx] <- idx.df$idr
        rep2.df <- dplyr::arrange(rep2.df, idr)
    } else {
        rep2.df$idr <- numeric(0)
    }

    return(list(rep1.df = rep1.df, rep2.df = rep2.df))
}
