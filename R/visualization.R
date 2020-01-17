
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
#' @importFrom ggplot2 element_text
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
        ggplot2::theme(panel.border = ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(face = "bold")) +
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
#' @importFrom ggplot2 element_text
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
        ggplot2::theme(panel.border = ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(face = "bold"),
                       legend.title = ggplot2::element_text(face = "bold")) +
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
#' @importFrom ggplot2 element_text
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
        ggplot2::theme(panel.border = ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(face = "bold"),
                       legend.title = ggplot2::element_text(face = "bold")) +
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

pretty_print_block_size <- function(x) {
    if (x < 10000) {
        return(paste0(x, " bp"))
    } else if (x < 1000000) {
        return(paste0(round(x / 1000, digits = 1), " kbp"))
    } else if (x < 1000000000) {
        return(paste0(round(x / 1000000, digits = 1), " Mbp"))
    } else {
        return(paste0(round(x / 1000000000, digits = 1), " Gbp"))
    }
}

pretty_print_reads <- function(x, values_normalized) {
    if (values_normalized) {
        return(paste0(round(x, digits = 2), "\nnormalized reads"))
    } else {
        return(vapply(x, function(reads) {
            if (reads < 1000000) {
                return(paste0(round(reads), " reads"))
            } else if (reads < 1000000000) {
                return(paste0(round(reads / 1000000, digits = 1),
                              " million reads"))
            } else {
                return(paste0(round(reads / 1000000000, digits = 1),
                              " billion reads"))
            }
        }, FUN.VALUE = character(1)))
    }
}

determine_block_size <- function(coords) {
    coords <- sort(unique(coords))
    return(min(vapply(seq_len(length(coords) - 1), function(i) {
        return(as.integer(abs(coords[i] - coords[i + 1])))
    }, FUN.VALUE = integer(1))))
}

#' @title Create Hi-C contact map
#'
#' @description
#' Creates Hi-C contact maps to visualize the results of
#' \code{\link{estimate_idr2d_hic}}.
#'
#' @param df output of \code{\link{estimate_idr2d_hic}}, a data frame with
#' the following columns:
#' \tabular{rll}{
#'   column 1: \tab \code{interaction} \tab character; genomic location of
#'   interaction block (e.g., \code{"chr1:204940000-204940000"})\cr
#'   column 2: \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the interactions\cr
#'   column 3: \tab \code{"rep_value"} \tab numeric; value of corresponding
#'   replicate interaction\cr
#'   column 4: \tab \code{"rank"} \tab integer; rank of the interaction,
#'   established by value column, ascending order\cr
#'   column 5: \tab \code{"rep_rank"} \tab integer; rank of corresponding
#'   replicate interaction\cr
#'   column 6: \tab \code{"idr"} \tab integer; IDR of the block and the
#'   corresponding block in the other replicate
#' }
#' @param idr_cutoff numeric; only show blocks with IDR < \code{idr_cutoff},
#' shows all blocks by default
#' @param chromosome character; chromsome name or list of chromosome names to
#' be analyzed, e.g., UCSC chromosome 1, \code{"chr1"}, defaults to all
#' chromosomes (\code{chromosome = NULL})
#' @param start_coordinate integer; only show contact map window between
#' \code{"start_coordinate"} and \code{"end_coordinate"}, by default shows
#' entire chromosome
#' @param end_coordinate integer; only show contact map window between
#' \code{"start_coordinate"} and \code{"end_coordinate"}, by default shows
#' entire chromosome
#' @param title character; plot title
#' @param values_normalized logical; are read counts in value column
#' raw or normalized? Defaults to \code{FALSE}
#' @param log_values logical; log-transform value column? Defaults to
#' \code{TRUE}
#'
#' @return ggplot2 object; Hi-C contact map
#'
#' @examples
#' idr_results_df <- estimate_idr2d_hic(idr2d:::hic$rep1_df,
#'                                      idr2d:::hic$rep2_df)
#' draw_hic_contact_map(idr_results_df, idr_cutoff = 0.05, chromosome = "chr1")
#'
#' @importFrom dplyr filter
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 vars
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_colorbar
#' @importFrom grid unit
#' @export
draw_hic_contact_map <- function(df, idr_cutoff = NULL,
                                 chromosome = NULL,
                                 start_coordinate = NULL,
                                 end_coordinate = NULL,
                                 title = NULL,
                                 values_normalized = FALSE,
                                 log_values = TRUE) {
    # avoid CRAN warnings
    idr <- start <- end <- chr <- log_value <- NULL

    if (nrow(df) == 0) {
        stop("data frame is empty")
    }

    if (is.factor(df$interaction)) {
        df$interaction <- as.character(df$interaction)
    }
    location <- strsplit(df$interaction, ":|-")
    df$chr <- vapply(location, function(r) {return(as.character(r[1]))},
                     FUN.VALUE = character(1))
    df$start <- vapply(location, function(r) {return(as.integer(r[2]))},
                       FUN.VALUE = integer(1))
    df$end <- vapply(location, function(r) {return(as.integer(r[3]))},
                     FUN.VALUE = integer(1))

    block_size <- determine_block_size(df$start)
    block_size_label <- pretty_print_block_size(block_size)

    if (!is.null(chromosome)) {
        df <- dplyr::filter(df, chr == chromosome)
        if (nrow(df) == 0) {
            stop(paste0("no blocks with chromosome name ", chromosome))
        }
    }
    if (!is.null(start_coordinate) && !is.null(end_coordinate)) {
        df <- dplyr::filter(df,
                            start >= start_coordinate & end <= end_coordinate)
        if (nrow(df) == 0) {
            stop(paste0("no blocks in window ",
                        start_coordinate, " - ", end_coordinate))
        }
    }

    df$start <- df$start / block_size
    df$end <- df$end / block_size

    if (!is.null(idr_cutoff)) {
        df <- dplyr::filter(df, idr < idr_cutoff)
        if (nrow(df) == 0) {
            stop(paste0("no blocks with IDR < ", idr_cutoff))
        }
    }

    if (log_values) {
        if (min(df$value) < 0) {
            futile.logger::flog.warn("log transform negative numbers")
        }
        df$log_value <- log(df$value + 1)
    } else {
        df$log_value <- df$value
    }

    lower_triangle_df <- df
    lower_triangle_df$start <- df$end
    lower_triangle_df$end <- df$start
    lower_triangle_df <- dplyr::filter(lower_triangle_df, start != end)
    df <- dplyr::bind_rows(df, lower_triangle_df)

    g <- ggplot2::ggplot(df, ggplot2::aes(x = start,
                                          y = end,
                                          fill = log_value)) +
        ggplot2::scale_fill_gradient(
            low = "white", high = "red",
            breaks = c(min(df$log_value), max(df$log_value)),
            labels = pretty_print_reads(c(min(df$value), max(df$value)),
                                        values_normalized)) +
        ggplot2::geom_tile() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank(),
                       axis.title = ggplot2::element_text(face = "bold"),
                       legend.title = ggplot2::element_text(face = "bold"),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.position = "top",
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank()) +
        ggplot2::labs(x = paste0("block (block size = ", block_size_label, ")"),
                      y = NULL,
                      fill = NULL, title = title) +
        ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = grid::unit(3, "inch")))

    if (length(unique(df$chr)) > 1) {
        g <- g + ggplot2::facet_grid(rows = ggplot2::vars(chr))
    }

    return(g)
}
