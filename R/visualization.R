
#' @title Create histogram of IDR values
#'
#' @description
#' Creates diagnostic plots to visualize the results of \code{estimateIDR}.
#'
#'
#' @param df part of output of \code{estimateIDR}, a data frame with at
#' least the following named columns:
#' \tabular{rl}{
#'   \code{idr} \tab IDR of the peak and the
#'   corresponding peak in the other replicate.
#' }
#' @param remove.na logical; should NA values be removed?
#' @param xlab character; x axis label
#' @param ylab character; y axis label
#' @param title character; plot title
#'
#' @return ggplot2 object; IDR distribution histogram
#'
#' @examples
#' idr.results <- estimateIDR1D(idr2d:::chipseq$rep1.df,
#'                              idr2d:::chipseq$rep2.df,
#'                              value.transformation = "log")
#' IDRDistributionHistogram(idr.results$rep1.df)
#'
#' @importFrom stats complete.cases
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @export
IDRDistributionHistogram <- function(df, remove.na = TRUE,
                                     xlab = "IDR",
                                     ylab = "density",
                                     title = "IDR value distribution") {
    # avoid CRAN warnings
    idr <- NULL

    if (remove.na) {
        df <- dplyr::filter(df, stats::complete.cases(idr))
    }

    g <- ggplot2::ggplot(df, ggplot2::aes(x = idr)) +
        ggplot2::geom_density() +
        ggplot2::scale_x_continuous(limits = c(0, 1.0)) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::labs(x = xlab, y = ylab, title = title)

    return(g)
}

#' @title Create scatterplot of IDR values
#'
#' @description
#' Creates diagnostic plots to visualize the results of \code{estimateIDR}.
#'
#' @param df part of output of \code{estimateIDR}, a data frame with at
#' least the following named columns:
#' \tabular{rl}{
#'   \code{rank} \tab integer; rank of the peak, established by
#'   value column, ascending order\cr
#'   \code{rep.rank} \tab integer; rank of corresponding
#'   replicate peak.\cr
#'   \code{idr} \tab IDR of the peak and the
#'   corresponding peak in the other replicate.
#' }
#' @param remove.na logical; should NA values be removed?
#' @param xlab character; x axis label
#' @param ylab character; y axis label
#' @param title character; plot title
#' @param color.gradient character; either "rainbow" or "default"
#' @param max.points.shown integer; default is 2500
#'
#' @return ggplot2 object; IDR rank scatterplot
#'
#' @examples
#' idr.results <- estimateIDR1D(idr2d:::chipseq$rep1.df,
#'                              idr2d:::chipseq$rep2.df,
#'                              value.transformation = "log")
#' rankIDRScatterplot(idr.results$rep1.df)
#'
#' @importFrom dplyr filter
#' @importFrom stats complete.cases
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom grDevices rainbow
#' @export
rankIDRScatterplot <- function(df, remove.na = TRUE,
                               xlab = "rank in replicate 1",
                               ylab = "rank in replicate 2",
                               title = "rank - IDR dependence",
                               color.gradient = c("rainbow", "default"),
                               max.points.shown = 2500) {
    # avoid CRAN warnings
    rank <- rep.rank <- idr <- NULL

    # argument handling
    color.gradient <- match.arg(color.gradient,
                                choices = c("rainbow", "default"))

    if (remove.na) {
        df <- dplyr::filter(df, stats::complete.cases(idr))
    }

    if (nrow(df) > max.points.shown) {
        df <- df[sample.int(nrow(df), max.points.shown), ]
    }

    g <- ggplot2::ggplot(df, ggplot2::aes(x = rank,
                                          y = rep.rank,
                                          color = idr)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::labs(x = xlab, y = ylab, color = "IDR", title = title)

    if (color.gradient == "rainbow") {
        g <- g + ggplot2::scale_color_gradientn(colours =
                                                    grDevices::rainbow(10),
                                                limits = c(0, 1.0),
                                                breaks = c(0.0, 0.25, 0.5,
                                                           0.75, 1.0))
    }
    return(g)
}

#' @title Create scatterplot of IDR values
#'
#' @description
#' Creates diagnostic plots to visualize the results of \code{estimateIDR}.
#'
#' @param df part of output of \code{estimateIDR}, a data frame with at
#' least the following named columns:
#' \tabular{rl}{
#'   \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the peaks\cr
#'   \code{rep.value} \tab numeric; value of corresponding
#'   replicate peak\cr
#'   \code{idr} \tab IDR of the peak and the
#'   corresponding peak in the other replicate.
#' }
#' @param remove.na logical; should NA values be removed?
#' @param remove.outliers logical; removes extreme data points
#' @param xlab character; x axis label
#' @param ylab character; y axis label
#' @param title character; plot title
#' @param color.gradient character; either "rainbow" or "default"
#' @param max.points.shown integer; default is 2500
#'
#' @return ggplot2 object; IDR value scatterplot
#'
#' @examples
#' idr.results <- estimateIDR1D(idr2d:::chipseq$rep1.df,
#'                              idr2d:::chipseq$rep2.df,
#'                              value.transformation = "log")
#' valueIDRScatterplot(idr.results$rep1.df)
#'
#' @importFrom dplyr filter
#' @importFrom stats complete.cases
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom grDevices rainbow
#' @export
valueIDRScatterplot <- function(df, remove.na = TRUE, remove.outliers = TRUE,
                                xlab = "transformed value in replicate 1",
                                ylab = "transformed value in replicate 2",
                                title = "value - IDR dependence",
                                color.gradient = c("rainbow", "default"),
                                max.points.shown = 2500) {
    # avoid CRAN warnings
    value <- rep.value <- idr <- NULL

    # argument handling
    color.gradient <- match.arg(color.gradient,
                                choices = c("rainbow", "default"))

    if (remove.na) {
        df <- dplyr::filter(df, stats::complete.cases(idr))
    }

    if (remove.outliers) {
        # remove outliers on upper tail (value)
        data.range <- range(df$value)[2] - range(df$value)[1]
        if (max(df$value) > 0) {
            truncated.df <- dplyr::filter(df, value < max(df$value) * 0.99)
        } else {
            truncated.df <- dplyr::filter(df, value < max(df$value) * 1.01)
        }

        truncated.data.range <- range(truncated.df$value)[2] -
            range(truncated.df$value)[1]

        if (truncated.data.range < 0.7 * data.range) {
            df <- truncated.df
        }

        # remove outliers on lower tail (value)
        data.range <- range(df$value)[2] - range(df$value)[1]
        if (min(df$value) > 0) {
            truncated.df <- dplyr::filter(df, value > min(df$value) * 1.01)
        } else {
            truncated.df <- dplyr::filter(df, value > min(df$value) * 0.99)
        }

        truncated.data.range <- range(truncated.df$value)[2] -
            range(truncated.df$value)[1]

        if (truncated.data.range < 0.7 * data.range) {
            df <- truncated.df
        }

        # remove outliers on upper tail (rep.value)
        data.range <- range(df$rep.value)[2] - range(df$rep.value)[1]
        if (max(df$rep.value) > 0) {
            truncated.df <- dplyr::filter(df, rep.value < max(df$rep.value) *
                                              0.99)
        } else {
            truncated.df <- dplyr::filter(df, rep.value < max(rep.value) * 1.01)
        }

        truncated.data.range <- range(truncated.df$rep.value)[2] -
            range(truncated.df$rep.value)[1]

        if (truncated.data.range < 0.7 * data.range) {
            df <- truncated.df
        }

        # remove outliers on lower tail (rep.value)
        data.range <- range(df$rep.value)[2] - range(df$rep.value)[1]
        if (min(df$rep.value) > 0) {
            truncated.df <- dplyr::filter(df, rep.value > min(df$rep.value) *
                                              1.01)
        } else {
            truncated.df <- dplyr::filter(df, rep.value > min(df$rep.value) *
                                              0.99)
        }

        truncated.data.range <- range(truncated.df$rep.value)[2] -
            range(truncated.df$rep.value)[1]

        if (truncated.data.range < 0.7 * data.range) {
            df <- truncated.df
        }
    }

    if (nrow(df) > max.points.shown) {
        df <- df[sample.int(nrow(df), max.points.shown), ]
    }

    g <- ggplot2::ggplot(df, ggplot2::aes(x = value,
                                          y = rep.value,
                                          color = idr)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::labs(x = xlab, y = ylab, color = "IDR", title = title)

    if (color.gradient == "rainbow") {
        g <- g + ggplot2::scale_color_gradientn(colours =
                                                    grDevices::rainbow(10),
                                                limits = c(0, 1.0),
                                                breaks = c(0.0, 0.25, 0.5,
                                                           0.75, 1.0))
    }
    return(g)
}