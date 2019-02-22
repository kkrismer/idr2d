
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
#' rep1.df$fdr <- preprocess(rep1.df$fdr, "log.additive.inverse")
#'
#' rep2.df <- idr2d:::chiapet$rep2.df
#' rep2.df$fdr <- preprocess(rep2.df$fdr, "log.additive.inverse")
#'
#' mapping <- establishBijection(rep1.df, rep2.df)
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr slice
#' @importFrom futile.logger flog.warn
#' @importFrom magrittr "%>%"
#' @export
establishBijection <- function(rep1.df, rep2.df,
                               ambiguity.resolution.method = c("value",
                                                               "overlap",
                                                               "midpoint"),
                               max.gap = 1000L) {
    # avoid CRAN warnings
    rep1.idx <- rep2.idx <- arv <- NULL

    ambiguity.resolution.method <- match.arg(ambiguity.resolution.method,
                                             choices = c("value",
                                                         "overlap",
                                                         "midpoint"))

    if (nrow(rep1.df) > 0 && nrow(rep2.df) > 0) {
        # shuffle to break preexisting order
        rep1.df <- rep1.df[sample.int(nrow(rep1.df)), ]
        rep2.df <- rep2.df[sample.int(nrow(rep2.df)), ]

        # sort by value column
        rep1.df <- dplyr::arrange(rep1.df, rep1.df[, 7])
        rep2.df <- dplyr::arrange(rep2.df, rep2.df[, 7])

        pairs.df <- overlap(rep1.df, rep2.df,
                            ambiguity.resolution.method, max.gap)

        top.pairs.df <- pairs.df %>% dplyr::group_by(rep1.idx) %>%
            dplyr::slice(which.min(arv))
        top.pairs.df <- top.pairs.df %>% dplyr::group_by(rep2.idx) %>%
            dplyr::slice(which.min(arv))

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
    } else {
        rep1.df$replicate.idx <- integer(0)
        rep2.df$replicate.idx <- integer(0)
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
#' @param local.idr TODO
#' @inheritParams idr::est.IDR
#' @inheritParams preprocess
#' @inheritParams establishBijection
#'
#' @return List with two components (\code{rep1.df} and \code{rep1.df})
#' containing the original data frames \code{rep1.df} and \code{rep2.df} with
#' additional column \code{idr}, which holds the IDR of the interaction and the
#' corresponding interaction in the other replicate. If no corresponding
#' interaction was found, \code{idr} is set to \code{NA}.
#'
#' @examples
#' idr.results <- estimateIDR(idr2d:::chiapet$rep1.df,
#'                            idr2d:::chiapet$rep2.df,
#'                            value.transformation = "log.additive.inverse")
#'
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom futile.logger flog.warn
#' @importFrom stringr str_trim
#' @export
estimateIDR <- function(rep1.df, rep2.df,
                        value.transformation = c("identity",
                                                 "additive.inverse",
                                                 "multiplicative.inverse",
                                                 "log",
                                                 "log.additive.inverse"),
                        ambiguity.resolution.method = c("value",
                                                        "overlap",
                                                        "midpoint"),
                        max.factor = 1.5,
                        jitter.factor = 0.0001,
                        max.gap = 1000L,
                        mu = 0.1, sigma = 1.0, rho = 0.2, p = 0.5,
                        eps = 0.001, max.iteration = 30, local.idr = TRUE) {
    # avoid CRAN warnings
    rep2.idx <- rep1.value <- rep2.value <- idr <- NULL

    rep1.df[, 7] <- preprocess(rep1.df[, 7], value.transformation,
                                max.factor = max.factor,
                                jitter.factor = jitter.factor)

    rep2.df[, 7] <- preprocess(rep2.df[, 7], value.transformation,
                                max.factor = max.factor,
                                jitter.factor = jitter.factor)

    mapping <- establishBijection(rep1.df, rep2.df, ambiguity.resolution.method,
                                  max.gap = max.gap)

    if (nrow(mapping$rep1.df) > 0 && nrow(mapping$rep2.df) > 0) {
        rep1.df <- mapping$rep1.df
        rep2.df <- mapping$rep2.df
        idx.df <- data.frame(
            rep1.idx = seq_len(nrow(rep1.df)),
            rep2.idx = rep1.df$replicate.idx,
            rep1.value = rep1.df[, 7],
            rep2.value = rep2.df[rep1.df$replicate.idx, 7]
        )
        rep1.df$replicate.idx <- NULL
        rep2.df$replicate.idx <- NULL

        idx.df <- dplyr::filter(idx.df, !is.na(rep2.idx) & !is.infinite(rep2.idx))

        if (nrow(idx.df) > 0) {
            idr.matrix <- as.matrix(dplyr::select(idx.df, rep1.value, rep2.value))

            if (length(unique(idr.matrix[, 1])) < 10 ||
                length(unique(idr.matrix[, 2])) < 10) {
                futile.logger::flog.warn("low complexity data set")
            }

            invisible(tryCatch({
                idr.results <- idr::est.IDR(idr.matrix, mu, sigma, rho, p,
                                            eps = eps,
                                            max.ite = max.iteration)
                # TODO check IDR vs idr?
                if(local.idr) {
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
            rep1.df$idx <- as.numeric(NA)
            rep1.df$rep.idx <- as.numeric(NA)
            rep1.df$idr <- as.numeric(NA)
            if (length(idx.df$idr) > 0) {
                rep1.df$idx[idx.df$rep1.idx] <- idx.df$rep1.idx
                rep1.df$rep.idx[idx.df$rep1.idx] <- idx.df$rep2.idx
                rep1.df$idr[idx.df$rep1.idx] <- idx.df$idr
                rep1.df <- dplyr::arrange(rep1.df, idr)
            }
        } else {
            rep1.df$idx <- numeric(0)
            rep1.df$rep.idx <- numeric(0)
            rep1.df$idr <- numeric(0)
        }

        if (nrow(rep2.df) > 0) {
            rep2.df$idx <- as.numeric(NA)
            rep2.df$rep.idx <- as.numeric(NA)
            rep2.df$idr <- as.numeric(NA)
            if (length(idx.df$idr) > 0) {
                rep2.df$idx[idx.df$rep2.idx] <- idx.df$rep2.idx
                rep2.df$rep.idx[idx.df$rep2.idx] <- idx.df$rep1.idx
                rep2.df$idr[idx.df$rep2.idx] <- idx.df$idr
                rep2.df <- dplyr::arrange(rep2.df, idr)
            }

        } else {
            rep2.df$idx <- numeric(0)
            rep2.df$rep.idx <- numeric(0)
            rep2.df$idr <- numeric(0)
        }
    } else {
        if (nrow(rep1.df) > 0) {
            rep1.df$idx <- as.numeric(NA)
            rep1.df$rep.idx <- as.numeric(NA)
            rep1.df$idr <- as.numeric(NA)
        } else {
            rep1.df$idx <- numeric(0)
            rep1.df$rep.idx <- numeric(0)
            rep1.df$idr <- numeric(0)
        }

        if (nrow(rep2.df) > 0) {
            rep2.df$idx <- as.numeric(NA)
            rep2.df$rep.idx <- as.numeric(NA)
            rep2.df$idr <- as.numeric(NA)
        } else {
            rep2.df$idx <- numeric(0)
            rep2.df$rep.idx <- numeric(0)
            rep2.df$idr <- numeric(0)
        }
    }

    return(list(rep1.df = rep1.df, rep2.df = rep2.df))
}
