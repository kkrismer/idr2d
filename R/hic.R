
getStandardChromosomes <- function(species, style) {
    system.file(package = "GenomeInfoDb", "extdata", "dataFiles")
}

#' @title TODO
#'
#' @description
#' TODO
#'
#' @param rep1.hic.file path to .hic file for replicate 1
#' @param rep1.hic.file path to .hic file for replicate 2
#' @param resolution block resolution of HiC contact matrix in base pairs,
#' defaults to 10,000 bp
#' @param chromosomes list of chromosome names in HiC files, defaults to UCSC
#' human chromosome names (chr1, ..., chr22, chrX, chrY, chrM)
#' @param use.python if Python is not on PATH, specify path to Python binary
#' here (see \code{\link[reticulate:use_python]{use_python}})
#' @param use.virtualenv if Python package \code{hic-straw} is not in base
#' virtualenv environment, specify environment here (see
#' \code{\link[reticulate:use_python]{use_virtualenv}})
#' @param use.condaenv if Python package \code{hic-straw} is not in base
#' conda environment, specify environment here (see
#' \code{\link[reticulate:use_python]{use_condaenv}})
#' @inheritParams estimateIDR2D
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
#' rep1.hic.file <- "D:/MIT/Research/PhD/Gifford/projects/chia-pet/data/GSE71831/GSE71831_Patski_paternal.hic"
#' rep2.hic.file <- "D:/MIT/Research/PhD/Gifford/projects/chia-pet/data/GSE71831/GSE71831_Patski_maternal.hic"
#' resolution <- 1000000
#' mus musculus
#' chromosomes <- paste0("chr", c(1:19, "X", "Y"))
#' df <- estimateIDR2DHiC(rep1.hic.file, rep2.hic.file, resolution, chromosomes = chromosomes)
#'
#' @export
estimateIDR2DHiC <- function(rep1.hic.file, rep2.hic.file, resolution = 10000,
                             normalization = c("NONE", "VC", "VC_SQRT", "KR"),
                             chromosomes = NULL,
                             max.factor = 1.5, jitter.factor = 0.0001,
                             mu = 0.1, sigma = 1.0, rho = 0.2, p = 0.5,
                             eps = 0.001, max.iteration = 30, local.idr = TRUE,
                             use.python = NULL, use.virtualenv = NULL,
                             use.condaenv = NULL) {
    value <- rep.value <- idr <- NULL

    normalization <- match.arg(normalization,
                               choices = c("NONE", "VC", "VC_SQRT", "KR"))
    if (is.null(chromosomes)) {
        chromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))
    }

    if (!is.null(use.python)) {
        reticulate::use_python(use.python)
    }

    if (!is.null(use.virtualenv)) {
        reticulate::use_virtualenv(use.virtualenv)
    }

    if (!is.null(use.condaenv)) {
        reticulate::use_condaenv(use.condaenv)
    }

    straw <- reticulate::import("straw")

    rep1 <- lapply(chromosomes, function(chromosome) {
        print(chromosome)
        counts <- straw$straw(normalization, rep1.hic.file,
                              chromosome, chromosome,
                              'BP', resolution)
        return(data.frame(interaction = paste0(
            chromosome, ":",
            format(counts[[1]], scientific = FALSE), "-",
            format(counts[[2]], scientific = FALSE)),
            value = counts[[3]], stringsAsFactors = FALSE))
    })
    rep1.df <- dplyr::bind_rows(rep1)

    rep2 <- lapply(chromosomes, function(chromosome) {
        counts <- straw$straw(normalization, rep2.hic.file,
                              chromosome, chromosome,
                              'BP', resolution)
        futile.logger::flog.info(paste0("completed ", chromosome))
        return(data.frame(interaction = paste0(
            chromosome, ":",
            format(counts[[1]], scientific = FALSE), "-",
            format(counts[[2]], scientific = FALSE)),
            rep.value = counts[[3]], stringsAsFactors = FALSE))
    })
    rep2.df <- dplyr::bind_rows(rep2)

    df <- dplyr::full_join(rep1.df, rep2.df, by = "interaction")

    df$value[is.na(df$value)] <- 0
    df$rep.value[is.na(df$rep.value)] <- 0

    df <- dplyr::arrange(df, dplyr::desc(value))
    df$rank <- seq_len(nrow(df))
    df <- dplyr::arrange(df, dplyr::desc(rep.value))
    df$rep.rank <- seq_len(nrow(df))

    idr.matrix <- as.matrix(dplyr::select(df,
                                          value,
                                          rep.value))

    invisible(tryCatch({
        idr.results <- idr::est.IDR(idr.matrix, mu, sigma, rho, p,
                                    eps = eps,
                                    max.ite = max.iteration)
        if (local.idr) {
            df$idr <- idr.results$idr
        } else {
            df$idr <- idr.results$IDR
        }
    }, error = function(e) {
        df$idr <- as.numeric(NA)
        futile.logger::flog.warn(stringr::str_trim(e))
    }))

    df <- dplyr::arrange(df, idr)

    return(df)
}
