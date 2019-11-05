
#' @title Create histogram of IDR values
#'
#' @description
#' Creates diagnostic plots to visualize the results of
#' \code{\link{estimate_idr}}.
#'
#' @param df part of output of \code{\link{estimate_idr}}, a data frame with at
#' least the following named columns:
#' \tabular{rl}{
#'   \code{idr} \tab IDR of the peak and the
#'   corresponding peak in the other replicate.
#' }
#' @param remove_na logical; should NA values be removed?
#' @param xlab character; x axis label
#' @param ylab character; y axis label
#' @param title character; plot title
#'
#' @return ggplot2 object; IDR distribution histogram
#'
#' @examples
#' idr_results <- estimate_idr1d(idr2d:::chipseq$rep1_df,
#'                               idr2d:::chipseq$rep2_df,
#'                               value_transformation = "log")
#' draw_idr_distribution_histogram(idr_results$rep1_df)
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
draw_idr_distribution_histogram <- function(df, remove_na = TRUE,
                                            xlab = "IDR",
                                            ylab = "density",
                                            title = "IDR value distribution") {
    # avoid CRAN warnings
    idr <- NULL

    if (remove_na) {
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
#' Creates diagnostic plots to visualize the results of
#' \code{\link{estimate_idr}}.
#'
#' @param df part of output of \code{\link{estimate_idr}}, a data frame with at
#' least the following named columns:
#' \tabular{rl}{
#'   \code{rank} \tab integer; rank of the peak, established by
#'   value column, ascending order\cr
#'   \code{rep_rank} \tab integer; rank of corresponding
#'   replicate peak.\cr
#'   \code{idr} \tab IDR of the peak and the
#'   corresponding peak in the other replicate.
#' }
#' @param remove_na logical; should NA values be removed?
#' @param xlab character; x axis label
#' @param ylab character; y axis label
#' @param log_idr logical; use logarithmized IDRs for colors to better
#' distinguish highly significant IDRs
#' @param title character; plot title
#' @param color_gradient character; either "rainbow" or "default"
#' @param alpha numeric; transparency of dots, from 0.0 - 1.0, where 1.0 is
#' completely opaque; default is 1.0
#' @param max_points_shown integer; default is 2500
#'
#' @return ggplot2 object; IDR rank scatterplot
#'
#' @examples
#' idr_results <- estimate_idr1d(idr2d:::chipseq$rep1_df,
#'                               idr2d:::chipseq$rep2_df,
#'                               value_transformation = "log")
#' draw_rank_idr_scatterplot(idr_results$rep1_df)
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
draw_rank_idr_scatterplot <- function(
    df, remove_na = TRUE,
    xlab = "rank in replicate 1",
    ylab = "rank in replicate 2",
    log_idr = FALSE,
    title = "rank - IDR dependence",
    color_gradient = c("rainbow", "default"),
    alpha = 1.0,
    max_points_shown = 2500) {
    # avoid CRAN warnings
    rank <- rep_rank <- idr <- NULL

    # argument handling
    color_gradient <- match.arg(color_gradient,
                                choices = c("rainbow", "default"))

    if (remove_na) {
        df <- dplyr::filter(df, stats::complete.cases(idr))
    }

    if (nrow(df) > max_points_shown) {
        df <- df[sample.int(nrow(df), max_points_shown), ]
    }

    if (log_idr) {
        df$idr <- log(df$idr)
        idr_legend_label <- "log(IDR)"
    } else {
        idr_legend_label <- "IDR"
    }

    g <- ggplot2::ggplot(df, ggplot2::aes(x = rank,
                                          y = rep_rank,
                                          color = idr),
                         alpha = alpha) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::labs(x = xlab, y = ylab,
                      color = idr_legend_label, title = title)

    if (color_gradient == "rainbow") {
        if (log_idr) {
            g <- g + ggplot2::scale_color_gradientn(colors =
                                                        grDevices::rainbow(10))
        } else {
            g <- g + ggplot2::scale_color_gradientn(colours =
                                                        grDevices::rainbow(10),
                                                    limits = c(0, 1.0),
                                                    breaks = c(0.0, 0.25, 0.5,
                                                               0.75, 1.0))
        }
    }
    return(g)
}

#' @title Create scatterplot of IDR values
#'
#' @description
#' Creates diagnostic plots to visualize the results of
#' \code{\link{estimate_idr}}.
#'
#' @param df part of output of \code{\link{estimate_idr}}, a data frame with at
#' least the following named columns:
#' \tabular{rl}{
#'   \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the peaks\cr
#'   \code{rep_value} \tab numeric; value of corresponding
#'   replicate peak\cr
#'   \code{idr} \tab IDR of the peak and the
#'   corresponding peak in the other replicate.
#' }
#' @param remove_na logical; should NA values be removed?
#' @param remove_outliers logical; removes extreme data points
#' @param xlab character; x axis label
#' @param ylab character; y axis label
#' @param log_axes logical; show logarithmized values from replicate 1 and 2
#' (default value is FALSE)
#' @param log_idr logical; use logarithmized IDRs for colors to better
#' distinguish highly significant IDRs
#' (default value is FALSE)
#' @param title character; plot title
#' @param color_gradient character; either "rainbow" or "default"
#' @param alpha numeric; transparency of dots, from 0.0 - 1.0, where 1.0 is
#' completely opaque; default is 1.0
#' @param max_points_shown integer; default is 2500
#'
#' @return ggplot2 object; IDR value scatterplot
#'
#' @examples
#' idr_results <- estimate_idr1d(idr2d:::chipseq$rep1_df,
#'                               idr2d:::chipseq$rep2_df,
#'                               value_transformation = "log")
#' draw_value_idr_scatterplot(idr_results$rep1_df)
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
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom scales log10_trans
#' @importFrom scales trans_breaks
#' @importFrom scales trans_format
#' @importFrom scales math_format
#' @export
draw_value_idr_scatterplot <- function(
    df, remove_na = TRUE, remove_outliers = TRUE,
    xlab = "transformed value in replicate 1",
    ylab = "transformed value in replicate 2",
    log_axes = FALSE,
    log_idr = FALSE,
    title = "value - IDR dependence",
    color_gradient = c("rainbow", "default"),
    alpha = 1.0,
    max_points_shown = 2500) {
    # avoid CRAN warnings
    value <- rep_value <- idr <- .x <- NULL

    # argument handling
    color_gradient <- match.arg(color_gradient,
                                choices = c("rainbow", "default"))

    if (remove_na) {
        df <- dplyr::filter(df, stats::complete.cases(idr))
    }

    if (remove_outliers) {
        # remove outliers on upper tail (value)
        data_range <- range(df$value)[2] - range(df$value)[1]
        if (max(df$value) > 0) {
            truncated_df <- dplyr::filter(df, value < max(df$value) * 0.99)
        } else {
            truncated_df <- dplyr::filter(df, value < max(df$value) * 1.01)
        }

        truncated_data_range <- range(truncated_df$value)[2] -
            range(truncated_df$value)[1]

        if (truncated_data_range < 0.7 * data_range) {
            df <- truncated_df
        }

        # remove outliers on lower tail (value)
        data_range <- range(df$value)[2] - range(df$value)[1]
        if (min(df$value) > 0) {
            truncated_df <- dplyr::filter(df, value > min(df$value) * 1.01)
        } else {
            truncated_df <- dplyr::filter(df, value > min(df$value) * 0.99)
        }

        truncated_data_range <- range(truncated_df$value)[2] -
            range(truncated_df$value)[1]

        if (truncated_data_range < 0.7 * data_range) {
            df <- truncated_df
        }

        # remove outliers on upper tail (rep_value)
        data_range <- range(df$rep_value)[2] - range(df$rep_value)[1]
        if (max(df$rep_value) > 0) {
            truncated_df <- dplyr::filter(df, rep_value < max(df$rep_value) *
                                              0.99)
        } else {
            truncated_df <- dplyr::filter(df, rep_value < max(rep_value) * 1.01)
        }

        truncated_data_range <- range(truncated_df$rep_value)[2] -
            range(truncated_df$rep_value)[1]

        if (truncated_data_range < 0.7 * data_range) {
            df <- truncated_df
        }

        # remove outliers on lower tail (rep_value)
        data_range <- range(df$rep_value)[2] - range(df$rep_value)[1]
        if (min(df$rep_value) > 0) {
            truncated_df <- dplyr::filter(df, rep_value > min(df$rep_value) *
                                              1.01)
        } else {
            truncated_df <- dplyr::filter(df, rep_value > min(df$rep_value) *
                                              0.99)
        }

        truncated_data_range <- range(truncated_df$rep_value)[2] -
            range(truncated_df$rep_value)[1]

        if (truncated_data_range < 0.7 * data_range) {
            df <- truncated_df
        }
    }

    if (nrow(df) > max_points_shown) {
        df <- df[sample.int(nrow(df), max_points_shown), ]
    }

    if (log_axes) {
        df$value <- df$value + 1
        df$rep_value <- df$rep_value + 1
    }

    if (log_idr) {
        df$idr <- log(df$idr)
        idr_legend_label <- "log(IDR)"
    } else {
        idr_legend_label <- "IDR"
    }

    g <- ggplot2::ggplot(df, ggplot2::aes(x = value,
                                          y = rep_value,
                                          color = idr),
                         alpha = alpha) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::labs(x = xlab, y = ylab,
                      color = idr_legend_label, title = title)

    if (log_axes) {
        g <- g + ggplot2::scale_x_continuous(
            trans = scales::log10_trans(),
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) + ggplot2::scale_y_continuous(
            trans = scales::log10_trans(),
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
    }

    if (color_gradient == "rainbow") {
        if (log_idr) {
            g <- g + ggplot2::scale_color_gradientn(colors =
                                                        grDevices::rainbow(10))
        } else {
            g <- g + ggplot2::scale_color_gradientn(colours =
                                                        grDevices::rainbow(10),
                                                    limits = c(0, 1.0),
                                                    breaks = c(0.0, 0.25, 0.5,
                                                               0.75, 1.0))
        }
    }
    return(g)
}
