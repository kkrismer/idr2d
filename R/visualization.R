
#' @export
IDRDistributionHistogram <- function(df, remove.na = TRUE,
                               xlab = "IDR",
                               ylab = "density",
                               title = "IDR value distribution") {
    # avoid CRAN warnings
    rank <- rep.rank <- idr <- NULL

    if (remove.na) {
        df <- dplyr::filter(df, complete.cases(idr))
    }

    g <- ggplot2::ggplot(df, ggplot2::aes(x = idr)) +
        ggplot2::geom_density() +
        ggplot2::scale_x_continuous(limits = c(0, 1.0)) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab(ylab) +
        ggplot2::ggtitle(title)

    return(g)
}

#' @export
rankIDRScatterplot <- function(df, remove.na = TRUE,
                               xlab = "rank in replicate 1",
                               ylab = "rank in replicate 2",
                               title = "rank - IDR dependence",
                               color.gradient = c("rainbow", "default")) {
    # avoid CRAN warnings
    rank <- rep.rank <- idr <- NULL

    # argument handling
    color.gradient <- match.arg(color.gradient,
                                choices = c("rainbow", "default"))

    if (remove.na) {
        df <- dplyr::filter(df, complete.cases(idr))
    }

    g <- ggplot2::ggplot(df, ggplot2::aes(x = rank,
                                          y = rep.rank,
                                          color = idr)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab(ylab) +
        ggplot2::labs(color = "IDR") +
        ggplot2::ggtitle(title)

    if (color.gradient == "rainbow") {
        g <- g + ggplot2::scale_color_gradientn(colours = grDevices::rainbow(10),
                                                limits = c(0, 1.0),
                                                breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))
    }
    return(g)
}

#' @export
valueIDRScatterplot <- function(df, remove.na = TRUE, remove.outliers = TRUE,
                               xlab = "transformed value in replicate 1",
                               ylab = "transformed value in replicate 2",
                               title = "value - IDR dependence",
                               color.gradient = c("rainbow", "default")) {
    # avoid CRAN warnings
    value <- rep.value <- idr <- NULL

    # argument handling
    color.gradient <- match.arg(color.gradient,
                                choices = c("rainbow", "default"))

    if (remove.na) {
        df <- dplyr::filter(df, complete.cases(idr))
    }

    if (remove.outliers) {
        # remove outliers on upper tail (value)
        data.range <- range(df$value)[2] - range(df$value)[1]
        if (max(df$value) > 0) {
            truncated.df <- dplyr::filter(df, value < max(df$value) * 0.99)
        } else {
            truncated.df <- dplyr::filter(df, value < max(df$value) * 1.01)
        }

        truncated.data.range <- range(truncated.df$value)[2] - range(truncated.df$value)[1]

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

        truncated.data.range <- range(truncated.df$value)[2] - range(truncated.df$value)[1]

        if (truncated.data.range < 0.7 * data.range) {
            df <- truncated.df
        }

        # remove outliers on upper tail (rep.value)
        data.range <- range(df$rep.value)[2] - range(df$rep.value)[1]
        if (max(df$rep.value) > 0) {
            truncated.df <- dplyr::filter(df, rep.value < max(df$rep.value) * 0.99)
        } else {
            truncated.df <- dplyr::filter(df, rep.value < max(rep.value) * 1.01)
        }

        truncated.data.range <- range(truncated.df$rep.value)[2] - range(truncated.df$rep.value)[1]

        if (truncated.data.range < 0.7 * data.range) {
            df <- truncated.df
        }

        # remove outliers on lower tail (rep.value)
        data.range <- range(df$rep.value)[2] - range(df$rep.value)[1]
        if (min(df$rep.value) > 0) {
            truncated.df <- dplyr::filter(df, rep.value > min(df$rep.value) * 1.01)
        } else {
            truncated.df <- dplyr::filter(df, rep.value > min(df$rep.value) * 0.99)
        }

        truncated.data.range <- range(truncated.df$rep.value)[2] - range(truncated.df$rep.value)[1]

        if (truncated.data.range < 0.7 * data.range) {
            df <- truncated.df
        }
    }

    g <- ggplot2::ggplot(df, ggplot2::aes(x = value,
                                          y = rep.value,
                                          color = idr)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank()) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab(ylab) +
        ggplot2::labs(color = "IDR") +
        ggplot2::ggtitle(title)

    if (color.gradient == "rainbow") {
        g <- g + ggplot2::scale_color_gradientn(colours = grDevices::rainbow(10),
                                                limits = c(0, 1.0),
                                                breaks = c(0.0, 0.25, 0.5, 0.75, 1.0))
    }
    return(g)
}
