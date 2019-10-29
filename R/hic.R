
#' @title Obtain Read Counts per Block in HiC Contact Matrices
#'
#' @description
#' \code{parse_hic_file} uses the Python package \code{hic-straw} internally
#'  to read .hic contact matrix files (see
#' \href{https://pypi.org/project/hic-straw/}{hic-straw on PyPI} or the
#' \href{https://pypi.org/project/hic-straw/}{Aiden lab GitHub repository}
#' for more information).
#'
#' The contact matrix is subdivided into blocks, where the block size is
#' determined by \code{resolution}. The reads per block are used to rank blocks
#' and replicate blocks are easily matched by genomic location.
#'
#' @references
#' Neva C. Durand, James T. Robinson, Muhammad S. Shamim, Ido Machol, Jill
#' P. Mesirov, Eric S. Lander, and Erez Lieberman Aiden. "Juicebox provides
#' a visualization system for Hi-C contact maps with unlimited zoom."
#' Cell Systems 3(1), 2016.
#'
#' @param hic_file path to .hic file (either local file path or URL).
#' @param resolution block resolution of HiC contact matrix in base pairs,
#' defaults to 1,000,000 bp (usually one of the following:
#' 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000)
#' @param normalization normalization step performed by Python package
#' \code{hic-straw}, one of the following: \code{"NONE"}, \code{"VC"},
#' \code{"VC_SQRT"}, \code{"KR"}.
#' @param chromosomes chromsome name or list of chromosome names in HiC file
#' to be analyzed, defaults to UCSC chromosome 1 (\code{"chr1"})
#' @param use_python if Python is not on PATH, specify path to Python binary
#' here (see \code{\link[reticulate:use_python]{use_python}})
#' @param use_virtualenv if Python package \code{hic-straw} is not in base
#' virtualenv environment, specify environment here (see
#' \code{\link[reticulate:use_python]{use_virtualenv}})
#' @param use_condaenv if Python package \code{hic-straw} is not in base
#' conda environment, specify environment here (see
#' \code{\link[reticulate:use_python]{use_condaenv}})
#'
#' @return Data frame with the following columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chr} \tab character; chromosome of block
#'   (e.g., \code{"chr3"})\cr
#'   column 2: \tab \code{region1} \tab integer; genomic location of side A of
#'   block in HiC contact matrix\cr
#'   column 3: \tab \code{region2} \tab integer; genomic location of side B of
#'   block in HiC contact matrix\cr
#'   column 4: \tab \code{value} \tab numeric; (normalized) read count in block
#' }
#'
#' @importFrom reticulate use_python
#' @importFrom reticulate use_virtualenv
#' @importFrom reticulate use_condaenv
#' @importFrom reticulate import
#' @importFrom futile.logger flog.info
#' @importFrom dplyr bind_rows
#' @export
parse_hic_file <- function(hic_file,
                           resolution = 1000000,
                           normalization = c("NONE", "VC", "VC_SQRT", "KR"),
                           chromosomes = NULL,
                           use_python = NULL, use_virtualenv = NULL,
                           use_condaenv = NULL) {
    normalization <- match.arg(normalization,
                               choices = c("NONE", "VC", "VC_SQRT", "KR"))
    if (is.null(chromosomes)) {
        chromosomes <- "chr1"
    }

    if (!is.null(use_python)) {
        reticulate::use_python(use_python)
    }

    if (!is.null(use_virtualenv)) {
        reticulate::use_virtualenv(use_virtualenv)
    }

    if (!is.null(use_condaenv)) {
        reticulate::use_condaenv(use_condaenv)
    }

    straw <- reticulate::import("straw")

    dfs <- lapply(chromosomes, function(chromosome) {
        counts <- tryCatch({
            straw$straw(normalization, hic_file,
                        chromosome, chromosome,
                        "BP", resolution)
        },
        error = function(e) {
            stop(paste0("error occurred while reading ", chromosome,
                        " of file '", hic_file, "': ", e,
                        "\npossible reasons:",
                        "\n(1) resolution might be too large or too small",
                        "\n(2) chromosome names are invalid, check species ",
                        "and chromosome name style - ",
                        "NCBI (e.g., \"1\"), UCSC (\"chr1\"), ",
                        "dbSNP (\"ch1\"), Ensembl (\"1\")",
                        "\n(3) Python package hic-straw is not installed",
                        "\n(4) .hic file does not exist"))
        })

        futile.logger::flog.info(paste0("read ", chromosome,
                                        ", replicate 1 (resolution = ",
                                        resolution, " bp)"))

        return(data.frame(chr = chromosome,
                          region1 = counts[[1]],
                          region2 = counts[[2]],
                          value = counts[[3]], stringsAsFactors = FALSE))
    })
    df <- dplyr::bind_rows(dfs)

    return(df)
}

#' @title Estimates IDR for Genomic Interactions measured by HiC experiments
#'
#' @description
#' This method estimates Irreproducible Discovery Rates (IDR) of
#' genomic interactions between two replicates of HiC experiments.
#'
#' The contact matrix is subdivided into blocks, where the block size is
#' determined by \code{resolution}. The reads per block are used to rank blocks
#' and replicate blocks are easily matched by genomic location.
#'
#' @references
#' Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011) Measuring
#' reproducibility of high-throughput experiments. Annals of Applied
#' Statistics, Vol. 5, No. 3, 1752-1779.
#'
#' @param rep1_df data frame of parsed .hic file for replicate 1 (output of
#' \code{\link{parse_hic_file}})
#' @param rep2_df data frame of parsed .hic file for replicate 2 (output of
#' \code{\link{parse_hic_file}})
#' @param combined_min_value exclude blocks with a combined (replicate 1 +
#' replicate 2) read count or normalized read count of less than
#'  \code{combined_min_value} (default is 20 reads)
#' @param combined_max_value exclude blocks with a combined (replicate 1 +
#' replicate 2) read count or normalized read count of more than
#'  \code{combined_max_value} (default is infinity)
#' @inheritParams estimate_idr2d
#'
#' @return Data frame with the following columns:
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
#'
#' @importFrom futile.logger flog.info
#' @importFrom stringr str_trim
#' @importFrom dplyr full_join
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr select
#' @importFrom idr est.IDR
#' @importFrom futile.logger flog.warn
#' @export
estimate_idr2d_hic <- function(rep1_df, rep2_df,
                               combined_min_value = 30,
                               combined_max_value = Inf,
                               max_factor = 1.5, jitter_factor = 0.0001,
                               mu = 0.1, sigma = 1.0, rho = 0.2, p = 0.5,
                               eps = 0.001, max_iteration = 30,
                               local_idr = TRUE) {
    value <- rep_value <- idr <- combined_value <- NULL

    rep1_df <- data.frame(interaction = paste0(
        rep1_df[, 1], ":",
        stringr::str_trim(format(rep1_df[, 2], scientific = FALSE)), "-",
        stringr::str_trim(format(rep1_df[, 3], scientific = FALSE))),
        value = rep1_df[, 4], stringsAsFactors = FALSE)

    rep2_df <- data.frame(interaction = paste0(
        rep2_df[, 1], ":",
        stringr::str_trim(format(rep2_df[, 2], scientific = FALSE)), "-",
        stringr::str_trim(format(rep2_df[, 3], scientific = FALSE))),
        rep_value = rep2_df[, 4], stringsAsFactors = FALSE)

    df <- dplyr::full_join(rep1_df, rep2_df, by = "interaction")
    df$value[is.na(df$value)] <- 0
    df$rep_value[is.na(df$rep_value)] <- 0

    if (is.infinite(combined_min_value) && is.infinite(combined_max_value)) {
        futile.logger::flog.info(paste0("number of interaction blocks: ",
                                        nrow(df)))
    } else {
        futile.logger::flog.info(paste0("number of interaction blocks ",
                                        "(before filtering): ", nrow(df)))

        df$combined_value <- df$value + df$rep_value
        if (!is.infinite(combined_min_value)) {
            df <- dplyr::filter(df, combined_value >= combined_min_value)
        }

        if (!is.infinite(combined_max_value)) {
            df <- dplyr::filter(df, combined_value <= combined_max_value)
        }
        df$combined_value <- NULL

        futile.logger::flog.info(paste0("number of interaction blocks ",
                                        "(after filtering): ", nrow(df)))
    }

    df <- dplyr::arrange(df, dplyr::desc(value))
    df$rank <- seq_len(nrow(df))
    df <- dplyr::arrange(df, dplyr::desc(rep_value))
    df$rep_rank <- seq_len(nrow(df))

    df$value <- preprocess(df$value, "identity",
                           max_factor = max_factor,
                           jitter_factor = jitter_factor)

    df$rep_value <- preprocess(df$rep_value, "identity",
                               max_factor = max_factor,
                               jitter_factor = jitter_factor)

    idr_matrix <- as.matrix(dplyr::select(df, value, rep_value))

    invisible(tryCatch({
        idr_results <- idr::est.IDR(idr_matrix, mu, sigma, rho, p,
                                    eps = eps,
                                    max.ite = max_iteration)
        if (local_idr) {
            df$idr <- idr_results$idr
        } else {
            df$idr <- idr_results$IDR
        }
    },
    error = function(e) {
        df$idr <- as.numeric(NA)
        futile.logger::flog.warn(stringr::str_trim(e))
    }))

    df <- dplyr::arrange(df, idr)

    class(df) <- c("idr2d_hic_result", class(df))
    return(df)
}


#' @importFrom methods is
#' @importFrom dplyr filter
#' @export
summary.idr2d_hic_result <- function(object, ...) {
    idr <- NULL

    stopifnot(methods::is(object, "idr2d_hic_result"))

    num_sig_df <- dplyr::filter(object, idr < 0.05)
    num_high_sig_df <- dplyr::filter(object, idr < 0.01)

    res <- list(
        analysis_type = "IDR2D HiC",
        num_blocks = nrow(object),
        num_significant_blocks = nrow(num_sig_df),
        num_highly_significant_blocks = nrow(num_high_sig_df)
    )

    class(res) <- c("idr2d_result_summary", class(res))
    return(res)
}
