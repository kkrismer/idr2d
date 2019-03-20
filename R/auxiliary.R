
#' @title Identifies Overlapping Anchors
#'
#' @description
#' Identifies all overlapping anchor pairs (m:n mapping).
#'
#' @param rep1.anchor data frame with the following columns:
#' \tabular{rl}{
#'   \code{chr} \tab character; genomic location of anchor in replicate 1 -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   \code{start} \tab integer; genomic location of anchor in replicate 1 -
#'   start coordinate\cr
#'   \code{end} \tab integer; genomic location of anchor in replicate 1 -
#'   end coordinate
#' }
#' @param rep2.anchor data frame with the following columns:
#' \tabular{rl}{
#'   \code{chr} \tab character; genomic location of anchor in replicate 2 -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   \code{start} \tab integer; genomic location of anchor in replicate 2 -
#'   start coordinate\cr
#'   \code{end} \tab integer; genomic location of anchor in replicate 2 -
#'   end coordinate
#' }
#' @param max.gap integer; maximum gap in nucleotides allowed between two
#' anchors for
#' them to be considered as overlapping (defaults to \code{1000L})
#'
#' @return A data frame where each row is an overlapping anchor pair.
#' There are two columns:
#' \tabular{rl}{
#'   \code{"rep1.idx"} \tab anchor index in data frame \code{rep1.anchor} \cr
#'   \code{"rep2.idx"} \tab anchor index in data frame \code{rep2.anchor}
#' }
#'
#' @examples
#' rep1.df <- idr2d:::chiapet$rep1.df
#' rep2.df <- idr2d:::chiapet$rep2.df
#'
#' rep1.anchor.a <- data.frame(chr = rep1.df[, 1],
#'                             start = rep1.df[, 2],
#'                             end = rep1.df[, 3])
#' rep2.anchor.a <- data.frame(chr = rep2.df[, 1],
#'                             start = rep2.df[, 2],
#'                             end = rep2.df[, 3])
#'
#' anchor.a.overlap <- anchorOverlap(rep1.anchor.a, rep2.anchor.a)
#'
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom stringr str_sort
#' @importFrom GenomicRanges findOverlaps
#' @export
anchorOverlap <- function(rep1.anchor, rep2.anchor, max.gap = 1000L) {
    rep1.ranges <- GenomicRanges::makeGRangesFromDataFrame(rep1.anchor)
    rep2.ranges <- GenomicRanges::makeGRangesFromDataFrame(rep2.anchor)

    # adjust seq levels
    seq.levels.r1 <- GenomeInfoDb::seqlevels(rep1.ranges)
    seq.levels.r2 <- GenomeInfoDb::seqlevels(rep2.ranges)
    combined.seq.levels <- stringr::str_sort(union(seq.levels.r1,
                                                   seq.levels.r2))
    GenomeInfoDb::seqlevels(rep1.ranges) <- combined.seq.levels
    GenomeInfoDb::seqlevels(rep2.ranges) <- combined.seq.levels

    # get overlap between replicates, accept 1000 bp gap
    overlap.df <- data.frame(GenomicRanges::findOverlaps(rep1.ranges,
                                                         rep2.ranges,
                                                         maxgap = max.gap))
    colnames(overlap.df) <- c("rep1.idx", "rep2.idx")
    return(overlap.df)
}

#' @title Distance between Midpoints of two Peaks
#'
#' @description
#' Calculates the distance in nucleotides between the midpoints of
#' two peaks.
#'
#' Note: peaks must be on the same
#' chromosome; start coordinate is always less than end coordinate
#'
#' @inheritParams calculateRelativeOverlap1D
#'
#' @return positive integer vector; distances between peak pairs
#'
#' @examples
#' # identical, zero distance
#' calculateMidpointDistance1D(100, 120,
#'                           100, 120)
#'
#' # centered, zero distance
#' calculateMidpointDistance1D(100, 120,
#'                           90, 130)
#'
#' # off by 10 per anchor
#' calculateMidpointDistance1D(100, 120,
#'                          110, 130)
#'
#' # vectorized example
#' calculateMidpointDistance1D(c(100, 100, 100),
#'                           c(120, 120, 120),
#'                           c(100, 90, 110),
#'                           c(120, 130, 130))
#'
#' @export
calculateMidpointDistance1D <- function(peak1.start, peak1.end,
                                        peak2.start, peak2.end) {
      midpoint.peak1 <- abs(peak1.start + (peak1.end - peak1.start) / 2)
      midpoint.peak2 <- abs(peak2.start + (peak2.end - peak2.start) / 2)
      return(as.integer(abs(midpoint.peak1 - midpoint.peak2)))
}

#' @title Distance between Anchor Midpoints of two Interactions
#'
#' @description
#' Calculates the distance in nucleotides between the anchor midpoints of
#' two interactions, which is the sum of the distance between midpoints of
#' anchor A in interaction 1 and anchor A in interaction 2, and the distance
#' between midpoints of
#' anchor B in interaction 1 and anchor B in interaction 2.
#'
#' Note: all anchors must be on the same
#' chromosome; start coordinate is always less than end coordinate
#'
#' @inheritParams calculateRelativeOverlap2D
#'
#' @return positive integer vector; distances between interaction pairs
#'
#' @examples
#' # identical, zero distance
#' calculateMidpointDistance2D(100, 120, 240, 260,
#'                           100, 120, 240, 260)
#'
#' # centered, zero distance
#' calculateMidpointDistance2D(100, 120, 240, 260,
#'                           90, 130, 230, 270)
#'
#' # off by 10 per anchor
#' calculateMidpointDistance2D(100, 120, 240, 250,
#'                          110, 130, 230, 240)
#'
#' # off by 10 (anchor B only)
#' calculateMidpointDistance2D(100, 120, 240, 250,
#'                          90, 130, 250, 260)
#'
#' # vectorized example
#' calculateMidpointDistance2D(c(100, 100, 100, 100),
#'                           c(120, 120, 120, 120),
#'                           c(240, 240, 240, 240),
#'                           c(260, 260, 250, 250),
#'                           c(100, 90, 110, 90),
#'                           c(120, 130, 130, 130),
#'                           c(240, 230, 230, 250),
#'                           c(260, 270, 240, 260))
#'
#' @export
calculateMidpointDistance2D <- function(int1.anchor.a.start, int1.anchor.a.end,
                                      int1.anchor.b.start, int1.anchor.b.end,
                                      int2.anchor.a.start, int2.anchor.a.end,
                                      int2.anchor.b.start, int2.anchor.b.end) {
    midpoint.int1.anchor.a <- abs(int1.anchor.a.start +
                                      (int1.anchor.a.end - int1.anchor.a.start) / 2)
    midpoint.int1.anchor.b <- abs(int1.anchor.b.start +
                                      (int1.anchor.b.end - int1.anchor.b.start) / 2)
    midpoint.int2.anchor.a <- abs(int2.anchor.a.start +
                                      (int2.anchor.a.end - int2.anchor.a.start) / 2)
    midpoint.int2.anchor.b <- abs(int2.anchor.b.start +
                                      (int2.anchor.b.end - int2.anchor.b.start) / 2)
    return(as.integer(abs(midpoint.int1.anchor.a - midpoint.int2.anchor.a) +
                          abs(midpoint.int1.anchor.b - midpoint.int2.anchor.b)))
}

#' @title Relative Anchor Overlap of two Peaks
#'
#' @description
#' Calculates the overlap between anchor A of interaction 1 and anchor
#' A of interaction 2, as well as anchor B of interaction 1 and anchor B of
#' interaction 2. The overlap (in nucleotides) is then normalized by the length
#' of the anchors.
#'
#' @param peak1.start integer vector; genomic start coordinate(s)
#' of peak in replicate 1
#' @param peak1.end integer vector; genomic end coordinate(s)
#' of peak in replicate 1
#' @param peak2.start integer vector; genomic start coordinate(s)
#' of peak in replicate 2
#' @param peak2.end integer vector; genomic end coordinate(s)
#' of peak in replicate 2
#'
#' @return numeric vector; relative overlaps between peak pairs
#'
#' @examples
#' # 100% overlap
#' calculateRelativeOverlap1D(100, 120,
#'                          100, 120)
#'
#' # 50% overlap
#' calculateRelativeOverlap1D(100, 120,
#'                          100, 110)
#'
#' # negative overlap
#' calculateRelativeOverlap1D(100, 120,
#'                          130, 140)
#'
#' # larger negative overlap
#' calculateRelativeOverlap1D(100, 120,
#'                          200, 220)
#'
#' # vectorized example
#' calculateRelativeOverlap1D(c(100, 100, 100, 100),
#'                          c(120, 120, 120, 120),
#'                          c(100, 100, 130, 200),
#'                          c(120, 110, 140, 220))
#' @export
calculateRelativeOverlap1D <- function(peak1.start, peak1.end,
                                       peak2.start, peak2.end) {
    peak.overlap <- pmin(peak1.end, peak2.end) -
        pmax(peak1.start, peak2.start)

    peak.combined.length <- pmax(peak1.end, peak2.end) -
        pmin(peak1.start, peak2.start)

    return(peak.overlap / peak.combined.length)
}

#' @title Relative Anchor Overlap of two Interactions
#'
#' @description
#' Calculates the overlap between anchor A of interaction 1 and anchor
#' A of interaction 2, as well as anchor B of interaction 1 and anchor B of
#' interaction 2. The overlap (in nucleotides) is then normalized by the length
#' of the anchors.
#'
#' Note: anchors A and B of the same interaction have to be on the same
#' chromosome; start coordinate is always less than end coordinate
#'
#' @param int1.anchor.a.start integer vector; genomic start coordinate(s)
#' of anchor A in replicate 1 interaction
#' @param int1.anchor.a.end integer vector; genomic end coordinate(s)
#' of anchor A in replicate 1 interaction
#' @param int1.anchor.b.start integer vector; genomic start coordinate(s)
#' of anchor B in replicate 1 interaction
#' @param int1.anchor.b.end integer vector; genomic end coordinate(s)
#' of anchor B in replicate 1 interaction
#' @param int2.anchor.a.start integer vector; genomic start coordinate(s)
#' of anchor A in replicate 2 interaction
#' @param int2.anchor.a.end integer vector; genomic end coordinate(s)
#' of anchor A in replicate 2 interaction
#' @param int2.anchor.b.start integer vector; genomic start coordinate(s)
#' of anchor B in replicate 2 interaction
#' @param int2.anchor.b.end integer vector; genomic end coordinate(s)
#' of anchor B in replicate 2 interaction
#'
#' @return numeric vector; relative overlaps between interaction pairs
#'
#' @examples
#' # 100% overlap
#' calculateRelativeOverlap2D(100, 120, 240, 260,
#'                          100, 120, 240, 260)
#'
#' # 50% overlap
#' calculateRelativeOverlap2D(100, 120, 240, 250,
#'                          100, 110, 240, 260)
#'
#' # negative overlap
#' calculateRelativeOverlap2D(100, 120, 240, 250,
#'                          130, 140, 260, 280)
#'
#' # larger negative overlap
#' calculateRelativeOverlap2D(100, 120, 240, 250,
#'                          200, 220, 340, 350)
#'
#' # vectorized example
#' calculateRelativeOverlap2D(c(100, 100, 100, 100),
#'                          c(120, 120, 120, 120),
#'                          c(240, 240, 240, 240),
#'                          c(260, 250, 250, 250),
#'                          c(100, 100, 130, 200),
#'                          c(120, 110, 140, 220),
#'                          c(240, 240, 260, 340),
#'                          c(260, 260, 280, 350))
#' @export
calculateRelativeOverlap2D <- function(int1.anchor.a.start, int1.anchor.a.end,
                                     int1.anchor.b.start, int1.anchor.b.end,
                                     int2.anchor.a.start, int2.anchor.a.end,
                                     int2.anchor.b.start, int2.anchor.b.end) {
    anchor.a.overlap <- pmin(int1.anchor.a.end, int2.anchor.a.end) -
        pmax(int1.anchor.a.start, int2.anchor.a.start)
    anchor.b.overlap <- pmin(int1.anchor.b.end, int2.anchor.b.end) -
        pmax(int1.anchor.b.start, int2.anchor.b.start)

    anchor.a.combined.length <- pmax(int1.anchor.a.end, int2.anchor.a.end) -
        pmin(int1.anchor.a.start, int2.anchor.a.start)
    anchor.b.combined.length <- pmax(int1.anchor.b.end, int2.anchor.b.end) -
        pmin(int1.anchor.b.start, int2.anchor.b.start)

    return((anchor.a.overlap + anchor.b.overlap) /
               (anchor.a.combined.length + anchor.b.combined.length))
}

#' @title Establish m:n Mapping Between Peaks from Replicate 1 and 2
#'
#' @description
#' This method returns all overlapping interactions between two replicates.
#' For each pair of overlapping interactions, the
#' \emph{ambiguity resolution value} (ARV) is calculated, which helps to reduce
#' the m:n mapping to a 1:1 mapping. The semantics of the ARV depend on the
#' specified \code{ambiguity.resolution.method}, but in general interaction
#' pairs with lower ARVs have priority over interaction pairs with higher ARVs
#' when the bijective mapping is established.
#'
#' @param rep1.df data frame of observations (i.e., genomic peaks) of
#' replicate 1, with at least the following columns (position of columns
#' matter, column names are irrelevant):
#' \tabular{rl}{
#'   column 1 (\code{chr}) \tab character; genomic location of peak -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2 \code{start}) \tab integer; genomic location of peak -
#'   start coordinate\cr
#'   column 3 (\code{end}) \tab integer; genomic location of peak -
#'   end coordinate\cr
#'   column 7 (\code{value}) \tab numeric; p-value, FDR, or heuristic used to
#'   rank the interactions
#' }
#' @param rep2.df data frame of observations (i.e., genomic peaks) of
#' replicate 2, with the following columns (position of columns
#' matter, column names are irrelevant):
#' \tabular{rl}{
#'   column 1 (\code{chr}) \tab character; genomic location of peak -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2 \code{start}) \tab integer; genomic location of peak -
#'   start coordinate\cr
#'   column 3 (\code{end}) \tab integer; genomic location of peak -
#'   end coordinate\cr
#'   column 7 (\code{value}) \tab numeric; p-value, FDR, or heuristic used to
#'   rank the interactions
#' }
#' @param ambiguity.resolution.method defines how ambiguous assignments
#' (when one interaction in replicate 1 overlaps with multiple interactions in
#' replicate 2 or vice versa)
#' are resolved. Available methods:
#' \tabular{rl}{
#'   \code{"value"} \tab interactions are prioritized by ascending or descending
#'   \code{value} column (see \code{sorting.direction}), e.g., if two
#'   interactions in replicate 1 overlap with one interaction in replicate 2,
#'   the interaction from replicate 1 is chosen which has a lower (if
#'   \code{sorting.direction} is \code{"ascending"}) or higher (if
#'   \code{"descending"}) value \cr
#'   \code{"overlap"} \tab the interaction pair is chosen which has the highest
#'   relative overlap, i.e., overlap in nucleotides of replicate 1 interaction
#'   anchor A and replicate 2 interaction anchor A,
#'   plus replicate 1 interaction anchor B and replicate 2 interaction anchor B,
#'   normalized by their lengths\cr
#'   \code{"midpoint"} \tab the interaction pair is chosen which has the
#'   smallest
#'   distance between their anchor midpoints, i.e., distance from midpoint of
#'   replicate 1 interaction anchor A to midpoint of
#'   replicate 2 interaction anchor A, plus distance from midpoint of
#'   replicate 1 interaction anchor B to midpoint of
#'   replicate 2 interaction anchor B
#' }
#' @inheritParams anchorOverlap
#'
#' @return data frame with the following columns:
#' \tabular{rl}{
#'   \code{rep1.idx} \tab index of interaction in replicate 1 (i.e., row
#'   index in \code{rep1.df})\cr
#'   \code{rep2.idx} \tab index of interaction in replicate 2 (i.e., row
#'   index in \code{rep2.df})\cr
#'   \code{arv} \tab ambiguity resolution value used turn m:n mapping into 1:1
#'   mapping. Interaction pairs with lower \code{arv} are prioritized.
#' }
#'
#' @examples
#' rep1.df <- idr2d:::chipseq$rep1.df
#' rep1.df$fdr <- preprocess(rep1.df$fdr, "log.additive.inverse")
#'
#' rep2.df <- idr2d:::chipseq$rep2.df
#' rep2.df$fdr <- preprocess(rep2.df$fdr, "log.additive.inverse")
#'
#' # shuffle to break preexisting order
#' rep1.df <- rep1.df[sample.int(nrow(rep1.df)), ]
#' rep2.df <- rep2.df[sample.int(nrow(rep2.df)), ]
#'
#' # sort by value column
#' rep1.df <- dplyr::arrange(rep1.df, rep1.df$fdr)
#' rep2.df <- dplyr::arrange(rep2.df, rep2.df$fdr)
#'
#' pairs.df <- overlap1D(rep1.df, rep2.df)
#'
#' @export
overlap1D <- function(rep1.df, rep2.df,
                    ambiguity.resolution.method = c("value", "overlap",
                                                    "midpoint"),
                    max.gap = 1000L) {
    # argument handling
    ambiguity.resolution.method <- match.arg(ambiguity.resolution.method,
                                             choices = c("value",
                                                         "overlap",
                                                         "midpoint"))

    if (nrow(rep1.df) > 0 && nrow(rep2.df) > 0) {
        overlaps.df <- anchorOverlap(
            data.frame(chr = rep1.df[, 1],
                       start = rep1.df[, 2],
                       end = rep1.df[, 3]),
            data.frame(chr = rep2.df[, 1],
                       start = rep2.df[, 2],
                       end = rep2.df[, 3]),
            max.gap = max.gap
        )

        idx.rep1 <- overlaps.df$rep1.idx
        idx.rep2 <- overlaps.df$rep2.idx

        # arv is ambiguity resolution value
        # the lower, the better
        if (ambiguity.resolution.method == "value") {
            arv <- rep1.df[idx.rep1, 4] + rep2.df[idx.rep2, 4]
            arv <- (-1) * arv
        } else if (ambiguity.resolution.method == "overlap") {
            arv <- calculateRelativeOverlap1D(
                rep1.df[idx.rep1, 2], rep1.df[idx.rep1, 3],
                rep2.df[idx.rep2, 2], rep2.df[idx.rep2, 3]
            )
            arv <- (-1) * arv
        } else if (ambiguity.resolution.method == "midpoint") {
            arv <- calculateMidpointDistance1D(
                rep1.df[idx.rep1, 2], rep1.df[idx.rep1, 3],
                rep2.df[idx.rep2, 2], rep2.df[idx.rep2, 3]
            )
        } else {
            stop(paste0("unknown ambiguity resolution method: ",
                        ambiguity.resolution.method))
        }

        pairs.df <- data.frame(rep1.idx = idx.rep1,
                               rep2.idx = idx.rep2,
                               arv = arv)
    } else {
        pairs.df <- data.frame()
    }

    return(pairs.df)
}


#' @title Establish m:n mapping between interactions from replicate 1 and 2
#'
#' @description
#' This method returns all overlapping interactions between two replicates.
#' For each pair of overlapping interactions, the
#' \emph{ambiguity resolution value} (ARV) is calculated, which helps to reduce
#' the m:n mapping to a 1:1 mapping. The semantics of the ARV depend on the
#' specified \code{ambiguity.resolution.method}, but in general interaction
#' pairs with lower ARVs have priority over interaction pairs with higher ARVs
#' when the bijective mapping is established.
#'
#' @param rep1.df data frame of observations (i.e., genomic interactions) of
#' replicate 1, with at least the following columns (position of columns
#' matter, column names are irrelevant):
#' \tabular{rl}{
#'   column 1 (\code{chrA}) \tab character; genomic location of anchor A -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2 \code{startA}) \tab integer; genomic location of anchor A -
#'   start coordinate\cr
#'   column 3 (\code{endA}) \tab integer; genomic location of anchor A -
#'   end coordinate\cr
#'   column 4 (\code{chrB}) \tab character; genomic location of anchor B -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 5 (\code{startB}) \tab integer; genomic location of anchor B -
#'   start coordinate\cr
#'   column 6 (\code{endB}) \tab integer; genomic location of anchor B -
#'   end coordinate\cr
#'   column 7 (\code{value}) \tab numeric; p-value, FDR, or heuristic used to
#'   rank the interactions
#' }
#' @param rep2.df data frame of observations (i.e., genomic interactions) of
#' replicate 2, with the following columns (position of columns
#' matter, column names are irrelevant):
#' \tabular{rl}{
#'   column 1 (\code{chrA}) \tab character; genomic location of anchor A -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2 \code{startA}) \tab integer; genomic location of anchor A -
#'   start coordinate\cr
#'   column 3 (\code{endA}) \tab integer; genomic location of anchor A -
#'   end coordinate\cr
#'   column 4 (\code{chrB}) \tab character; genomic location of anchor B -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 5 (\code{startB}) \tab integer; genomic location of anchor B -
#'   start coordinate\cr
#'   column 6 (\code{endB}) \tab integer; genomic location of anchor B -
#'   end coordinate\cr
#'   column 7 (\code{value}) \tab numeric; p-value, FDR, or heuristic used to
#'   rank the interactions
#' }
#' @param ambiguity.resolution.method defines how ambiguous assignments
#' (when one interaction in replicate 1 overlaps with multiple interactions in
#' replicate 2 or vice versa)
#' are resolved. Available methods:
#' \tabular{rl}{
#'   \code{"value"} \tab interactions are prioritized by ascending or descending
#'   \code{value} column (see \code{sorting.direction}), e.g., if two
#'   interactions in replicate 1 overlap with one interaction in replicate 2,
#'   the interaction from replicate 1 is chosen which has a lower (if
#'   \code{sorting.direction} is \code{"ascending"}) or higher (if
#'   \code{"descending"}) value \cr
#'   \code{"overlap"} \tab the interaction pair is chosen which has the highest
#'   relative overlap, i.e., overlap in nucleotides of replicate 1 interaction
#'   anchor A and replicate 2 interaction anchor A,
#'   plus replicate 1 interaction anchor B and replicate 2 interaction anchor B,
#'   normalized by their lengths\cr
#'   \code{"midpoint"} \tab the interaction pair is chosen which has the
#'   smallest
#'   distance between their anchor midpoints, i.e., distance from midpoint of
#'   replicate 1 interaction anchor A to midpoint of
#'   replicate 2 interaction anchor A, plus distance from midpoint of
#'   replicate 1 interaction anchor B to midpoint of
#'   replicate 2 interaction anchor B
#' }
#' @inheritParams anchorOverlap
#'
#' @return data frame with the following columns:
#' \tabular{rl}{
#'   \code{rep1.idx} \tab index of interaction in replicate 1 (i.e., row
#'   index in \code{rep1.df})\cr
#'   \code{rep2.idx} \tab index of interaction in replicate 2 (i.e., row
#'   index in \code{rep2.df})\cr
#'   \code{arv} \tab ambiguity resolution value used turn m:n mapping into 1:1
#'   mapping. Interaction pairs with lower \code{arv} are prioritized.
#' }
#'
#' @examples
#' rep1.df <- idr2d:::chiapet$rep1.df
#' rep1.df$fdr <- preprocess(rep1.df$fdr, "log.additive.inverse")
#'
#' rep2.df <- idr2d:::chiapet$rep2.df
#' rep2.df$fdr <- preprocess(rep2.df$fdr, "log.additive.inverse")
#'
#' # shuffle to break preexisting order
#' rep1.df <- rep1.df[sample.int(nrow(rep1.df)), ]
#' rep2.df <- rep2.df[sample.int(nrow(rep2.df)), ]
#'
#' # sort by value column
#' rep1.df <- dplyr::arrange(rep1.df, rep1.df$fdr)
#' rep2.df <- dplyr::arrange(rep2.df, rep2.df$fdr)
#'
#' pairs.df <- overlap2D(rep1.df, rep2.df)
#'
#' @export
overlap2D <- function(rep1.df, rep2.df,
                    ambiguity.resolution.method = c("value", "overlap",
                                                    "midpoint"),
                    max.gap = 1000L) {
    # argument handling
    ambiguity.resolution.method <- match.arg(ambiguity.resolution.method,
                                             choices = c("value",
                                                         "overlap",
                                                         "midpoint"))
    if (nrow(rep1.df) > 0 && nrow(rep2.df) > 0) {
        overlaps.anchorsA.df <- anchorOverlap(
            data.frame(chr = rep1.df[, 1],
                       start = rep1.df[, 2],
                       end = rep1.df[, 3]),
            data.frame(chr = rep2.df[, 1],
                       start = rep2.df[, 2],
                       end = rep2.df[, 3]),
            max.gap = max.gap
        )
        overlaps.anchorsB.df <- anchorOverlap(
            data.frame(chr = rep1.df[, 4],
                       start = rep1.df[, 5],
                       end = rep1.df[, 6]),
            data.frame(chr = rep2.df[, 4],
                       start = rep2.df[, 5],
                       end = rep2.df[, 6]),
            max.gap = max.gap
        )

        a <- paste0(overlaps.anchorsA.df$rep1.idx, "-",
                    overlaps.anchorsA.df$rep2.idx)
        b <- paste0(overlaps.anchorsB.df$rep1.idx, "-",
                    overlaps.anchorsB.df$rep2.idx)

        replicates <- intersect(a, b)

        idx.rep1 <- unlist(lapply(strsplit(replicates, "-", fixed = TRUE),
                                  function(interaction) {
                                      return(as.integer(interaction[1]))
                                  }))

        idx.rep2 <- unlist(lapply(strsplit(replicates, "-", fixed = TRUE),
                                  function(interaction) {
                                      return(as.integer(interaction[2]))
                                  }))

        # arv is ambiguity resolution value
        # the lower, the better
        if (ambiguity.resolution.method == "value") {
            arv <- rep1.df[idx.rep1, 7] + rep2.df[idx.rep2, 7]
            arv <- (-1) * arv
        } else if (ambiguity.resolution.method == "overlap") {
            arv <- calculateRelativeOverlap2D(
                rep1.df[idx.rep1, 2], rep1.df[idx.rep1, 3],
                rep1.df[idx.rep1, 5], rep1.df[idx.rep1, 6],
                rep2.df[idx.rep2, 2], rep2.df[idx.rep2, 3],
                rep2.df[idx.rep2, 5], rep2.df[idx.rep2, 6]
            )
            arv <- (-1) * arv
        } else if (ambiguity.resolution.method == "midpoint") {
            arv <- calculateMidpointDistance2D(
                rep1.df[idx.rep1, 2], rep1.df[idx.rep1, 3],
                rep1.df[idx.rep1, 5], rep1.df[idx.rep1, 6],
                rep2.df[idx.rep2, 2], rep2.df[idx.rep2, 3],
                rep2.df[idx.rep2, 5], rep2.df[idx.rep2, 6]
            )
        } else {
            stop(paste0("unknown ambiguity resolution method: ",
                        ambiguity.resolution.method))
        }

        pairs.df <- data.frame(rep1.idx = idx.rep1,
                               rep2.idx = idx.rep2,
                               arv = arv)
    } else {
        pairs.df <- data.frame()
    }

    return(pairs.df)
}


#' @title Removes Peaks on Non-standard Chromosomes
#' @param x data frame of genomic peaks, with the following columns
#' (position of columns matter, column names are irrelevant):
#' \tabular{rl}{
#'   column 1 (\code{chr}) \tab character; genomic location of peak -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2 \code{start}) \tab integer; genomic location of peak  -
#'   start coordinate\cr
#'   column 3 (\code{end}) \tab integer; genomic location of peak -
#'   end coordinate\cr
#'   column 4 (\code{value}) \tab numeric; p-value, FDR, or heuristic used to
#'   rank the peaks
#' }
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @export
removeNonstandardChromosomes1D <- function(x) {
    x <- GenomicRanges::GRanges(x$chr,
                                IRanges::IRanges(x$start, x$end),
                                value = x$value)
    x <- GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = "coarse")
    x.df <- as.data.frame(x)
    x.df$width <- NULL
    x.df$strand <- NULL
    colnames(x.df) <- c("chr", "start", "end", "value")
    return(x.df)
}

#' @title Removes Interactions on Non-standard Chromosomes
#' @param x data frame of genomic interactions, with the following columns
#' (position of columns matter, column names are irrelevant):
#' \tabular{rl}{
#'   column 1 (\code{chr.a}) \tab character; genomic location of anchor A -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2 \code{start.a}) \tab integer; genomic location of anchor A -
#'   start coordinate\cr
#'   column 3 (\code{end.a}) \tab integer; genomic location of anchor A -
#'   end coordinate\cr
#'   column 4 (\code{chr.b}) \tab character; genomic location of anchor B -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 5 (\code{start.b}) \tab integer; genomic location of anchor B -
#'   start coordinate\cr
#'   column 6 (\code{end.c}) \tab integer; genomic location of anchor B -
#'   end coordinate\cr
#'   column 7 (\code{value}) \tab numeric; p-value, FDR, or heuristic used to
#'   rank the interactions
#' }
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom dplyr inner_join
#' @importFrom dplyr select
#' @export
removeNonstandardChromosomes2D <- function(x) {
    # avoid CRAN warnings
    chr.a <- start.a <- end.a <- chr.b <- start.b <- end.b <- value <- NULL

    x$idx <- seq_len(nrow(x))
    anchor.a <- GenomicRanges::GRanges(x$chr.a,
                                       IRanges::IRanges(x$start.a, x$end.a),
                                       value = x$value,
                                       idx = x$idx)
    anchor.a <- GenomeInfoDb::keepStandardChromosomes(anchor.a,
                                                      pruning.mode = "coarse")
    anchor.a.df <- as.data.frame(anchor.a)
    anchor.a.df$width <- NULL
    anchor.a.df$strand <- NULL
    colnames(anchor.a.df) <- c("chr.a", "start.a", "end.a", "value", "idx")

    anchor.b <- GenomicRanges::GRanges(x$chr.b,
                                       IRanges::IRanges(x$start.b, x$end.b),
                                       idx = x$idx)
    anchor.b <- GenomeInfoDb::keepStandardChromosomes(anchor.b,
                                                      pruning.mode = "coarse")
    anchor.b.df <- as.data.frame(anchor.b)
    anchor.b.df$width <- NULL
    anchor.b.df$strand <- NULL
    colnames(anchor.b.df) <- c("chr.b", "start.b", "end.b", "idx")

    x.df <- dplyr::inner_join(anchor.a.df, anchor.b.df, by = "idx")
    x.df <- dplyr::select(x.df, chr.a, start.a, end.a,
                          chr.b, start.b, end.b, value)

    return(x.df)
}
