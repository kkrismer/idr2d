
get_standard_chromosomes <- function(species, style) {
    system.file(package = "GenomeInfoDb", "extdata", "dataFiles")
}

#' @title Estimates IDR for Genomic Interactions measured by HiC experiments
#'
#' @description
#' This method estimates Irreproducible Discovery Rates (IDR) of
#' genomic interactions between two replicates of HiC experiments.
#'
#' \code{estimate_idr2d_hic} uses the Python package \code{hic-straw} internally
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
#' Q. Li, J. B. Brown, H. Huang and P. J. Bickel. (2011) Measuring
#' reproducibility of high-throughput experiments. Annals of Applied
#' Statistics, Vol. 5, No. 3, 1752-1779.
#'
#' @param rep1_hic_file path to .hic file for replicate 1
#' (either local file path or URL)
#' @param rep2_hic_file path to .hic file for replicate 2
#' (either local file path or URL)
#' @param resolution block resolution of HiC contact matrix in base pairs,
#' defaults to 10,000 bp
#' @param normalization normalization step performed by Python package
#' \code{hic-straw}, one of the following: \code{"NONE"}, \code{"VC"},
#' \code{"VC_SQRT"}, \code{"KR"}.
#' @param chromosomes list of chromosome names in HiC files, defaults to UCSC
#' human chromosome names (chr1, ..., chr22, chrX, chrY, chrM)
#' @param use_python if Python is not on PATH, specify path to Python binary
#' here (see \code{\link[reticulate:use_python]{use_python}})
#' @param use_virtualenv if Python package \code{hic-straw} is not in base
#' virtualenv environment, specify environment here (see
#' \code{\link[reticulate:use_python]{use_virtualenv}})
#' @param use_condaenv if Python package \code{hic-straw} is not in base
#' conda environment, specify environment here (see
#' \code{\link[reticulate:use_python]{use_condaenv}})
#' @inheritParams estimate_idr2d
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
#' @importFrom reticulate use_python
#' @importFrom reticulate use_virtualenv
#' @importFrom reticulate use_condaenv
#' @importFrom reticulate import
#' @importFrom futile.logger flog.info
#' @importFrom stringr str_trim
#' @importFrom dplyr bind_rows
#' @importFrom dplyr full_join
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom dplyr select
#' @importFrom idr est.IDR
#' @importFrom futile.logger flog.warn
#' @export
estimate_idr2d_hic <- function(rep1_hic_file, rep2_hic_file, resolution = 10000,
                               normalization = c("NONE", "VC", "VC_SQRT", "KR"),
                               chromosomes = NULL,
                               max_factor = 1.5, jitter_factor = 0.0001,
                               mu = 0.1, sigma = 1.0, rho = 0.2, p = 0.5,
                               eps = 0.001, max_iteration = 30,
                               local_idr = TRUE,
                               use_python = NULL, use_virtualenv = NULL,
                               use_condaenv = NULL) {
    value <- rep_value <- idr <- NULL

    normalization <- match.arg(normalization,
                               choices = c("NONE", "VC", "VC_SQRT", "KR"))
    if (is.null(chromosomes)) {
        chromosomes <- paste0("chr", c(seq_len(22), "X", "Y", "M"))
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

    rep1 <- lapply(chromosomes, function(chromosome) {
        counts <- tryCatch({
            straw$straw(normalization, rep1_hic_file,
                        chromosome, chromosome,
                        "BP", resolution)
        },
        error = function(e) {
            stop(paste0("error occurred while reading ", chromosome,
                        " of replicate 1: ", e,
                        "\npossible reasons:",
                        "\n(1) resolution might be too large or too small",
                        "\n(2) chromosome names are invalid, check species ",
                        "and chromosome name style - ",
                        "NCBI (e.g., \"1\"), UCSC (\"chr1\"), ",
                        "dbSNP (\"ch1\"), Ensembl (\"1\")",
                        "\n(3) Python package hic-straw is not installed"))
        })

        futile.logger::flog.info(paste0("read ", chromosome,
                                        ", replicate 1 (resolution = ",
                                        resolution, " bp)"))
        return(data.frame(interaction = paste0(
            chromosome, ":",
            stringr::str_trim(format(counts[[1]], scientific = FALSE)), "-",
            stringr::str_trim(format(counts[[2]], scientific = FALSE))),
            value = counts[[3]], stringsAsFactors = FALSE))
    })
    rep1_df <- dplyr::bind_rows(rep1)

    rep2 <- lapply(chromosomes, function(chromosome) {
        counts <- tryCatch({
            straw$straw(normalization, rep2_hic_file,
                        chromosome, chromosome,
                        "BP", resolution)
        },
        error = function(e) {
            stop(paste0("error occurred while reading ", chromosome,
                        " of replicate 2: ", e,
                        "\nresolution might be too large or too small"))
        })

        futile.logger::flog.info(paste0("read ", chromosome,
                                        ", replicate 2 (resolution = ",
                                        resolution, " bp)"))
        return(data.frame(interaction = paste0(
            chromosome, ":",
            stringr::str_trim(format(counts[[1]], scientific = FALSE)), "-",
            stringr::str_trim(format(counts[[2]], scientific = FALSE))),
            rep_value = counts[[3]], stringsAsFactors = FALSE))
    })
    rep2_df <- dplyr::bind_rows(rep2)

    df <- dplyr::full_join(rep1_df, rep2_df, by = "interaction")

    futile.logger::flog.info(paste0("number of interaction blocks: ", nrow(df)))

    df$value[is.na(df$value)] <- 0
    df$rep_value[is.na(df$rep_value)] <- 0

    df <- dplyr::arrange(df, dplyr::desc(value))
    df$rank <- seq_len(nrow(df))
    df <- dplyr::arrange(df, dplyr::desc(rep_value))
    df$rep_rank <- seq_len(nrow(df))

    idr_matrix <- as.matrix(dplyr::select(df,
                                          value,
                                          rep_value))

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

    return(df)
}
