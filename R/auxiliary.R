
#' @title Identifies Overlapping Anchors
#'
#' @description
#' Identifies all overlapping anchor pairs (m:n mapping).
#'
#' @param rep1_anchor data frame with the following columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chr} \tab character; genomic location of anchor in
#'   replicate 1 - chromosome (e.g., \code{"chr3"})\cr
#'   column 2: \tab \code{start} \tab integer; genomic location of anchor in
#'   replicate 1 - start coordinate\cr
#'   column 3: \tab \code{end} \tab integer; genomic location of anchor in
#'   replicate 1 - end coordinate
#' }
#' @param rep2_anchor data frame with the following columns:
#' \tabular{rll}{
#'   column 1: \tab \code{chr} \tab character; genomic location of anchor in
#'   replicate 2 - chromosome (e.g., \code{"chr3"})\cr
#'   column 2: \tab \code{start} \tab integer; genomic location of anchor in
#'   replicate 2 - start coordinate\cr
#'   column 3: \tab \code{end} \tab integer; genomic location of anchor in
#'   replicate 2 - end coordinate
#' }
#' @param max_gap integer; maximum gap in nucleotides allowed between two
#' anchors for
#' them to be considered as overlapping
#' (defaults to zero, no gap between anchors)
#'
#' @return A data frame containing overlapping
#'  anchor pairs with the following columns:
#' \tabular{rll}{
#'   column 1: \tab \code{rep1_idx} \tab anchor index in data frame
#'   \code{rep1_anchor} \cr
#'   column 2: \tab \code{rep2_idx} \tab anchor index in data frame
#'   \code{rep2_anchor}
#' }
#'
#' @examples
#' rep1_df <- idr2d:::chiapet$rep1_df
#' rep2_df <- idr2d:::chiapet$rep2_df
#'
#' rep1_anchor_a <- data.frame(chr = rep1_df[, 1],
#'                             start = rep1_df[, 2],
#'                             end = rep1_df[, 3])
#' rep2_anchor_a <- data.frame(chr = rep2_df[, 1],
#'                             start = rep2_df[, 2],
#'                             end = rep2_df[, 3])
#'
#' anchor_a_overlap <- anchorOverlap(rep1_anchor_a, rep2_anchor_a)
#'
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom stringr str_sort
#' @importFrom GenomicRanges findOverlaps
#' @export
anchorOverlap <- function(rep1_anchor, rep2_anchor, max_gap = 0L) {
    rep1_ranges <- GenomicRanges::makeGRangesFromDataFrame(rep1_anchor)
    rep2_ranges <- GenomicRanges::makeGRangesFromDataFrame(rep2_anchor)

    # adjust seq levels
    seq_levels_r1 <- GenomeInfoDb::seqlevels(rep1_ranges)
    seq_levels_r2 <- GenomeInfoDb::seqlevels(rep2_ranges)
    combined_seq_levels <- stringr::str_sort(union(seq_levels_r1,
                                                   seq_levels_r2))
    GenomeInfoDb::seqlevels(rep1_ranges) <- combined_seq_levels
    GenomeInfoDb::seqlevels(rep2_ranges) <- combined_seq_levels

    # get overlap between replicates, accept 1000 bp gap
    overlap_df <- data.frame(GenomicRanges::findOverlaps(rep1_ranges,
                                                         rep2_ranges,
                                                         maxgap = max_gap))
    colnames(overlap_df) <- c("rep1_idx", "rep2_idx")
    return(overlap_df)
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
calculateMidpointDistance1D <- function(peak1_start, peak1_end,
                                        peak2_start, peak2_end) {
    midpoint_peak1 <- abs(peak1_start + (peak1_end - peak1_start) / 2)
    midpoint_peak2 <- abs(peak2_start + (peak2_end - peak2_start) / 2)
    return(as.integer(abs(midpoint_peak1 - midpoint_peak2)))
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
calculateMidpointDistance2D <- function(int1_anchor_a_start,
                                        int1_anchor_a_end,
                                        int1_anchor_b_start,
                                        int1_anchor_b_end,
                                        int2_anchor_a_start,
                                        int2_anchor_a_end,
                                        int2_anchor_b_start,
                                        int2_anchor_b_end) {
    midpoint_int1_anchor_a <- abs(int1_anchor_a_start +
                                      (int1_anchor_a_end -
                                           int1_anchor_a_start) / 2)
    midpoint_int1_anchor_b <- abs(int1_anchor_b_start +
                                      (int1_anchor_b_end -
                                           int1_anchor_b_start) / 2)
    midpoint_int2_anchor_a <- abs(int2_anchor_a_start +
                                      (int2_anchor_a_end -
                                           int2_anchor_a_start) / 2)
    midpoint_int2_anchor_b <- abs(int2_anchor_b_start +
                                      (int2_anchor_b_end -
                                           int2_anchor_b_start) / 2)
    return(as.integer(abs(midpoint_int1_anchor_a - midpoint_int2_anchor_a) +
                          abs(midpoint_int1_anchor_b - midpoint_int2_anchor_b)))
}

#' @title Relative Anchor Overlap of two Peaks
#'
#' @description
#' Calculates the overlap between anchor A of interaction 1 and anchor
#' A of interaction 2, as well as anchor B of interaction 1 and anchor B of
#' interaction 2. The overlap (in nucleotides) is then normalized by the length
#' of the anchors.
#'
#' @param peak1_start integer vector; genomic start coordinate(s)
#' of peak in replicate 1
#' @param peak1_end integer vector; genomic end coordinate(s)
#' of peak in replicate 1
#' @param peak2_start integer vector; genomic start coordinate(s)
#' of peak in replicate 2
#' @param peak2_end integer vector; genomic end coordinate(s)
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
calculateRelativeOverlap1D <- function(peak1_start, peak1_end,
                                       peak2_start, peak2_end) {
    peak_overlap <- pmin(peak1_end, peak2_end) -
        pmax(peak1_start, peak2_start)

    peak_combined_length <- pmax(peak1_end, peak2_end) -
        pmin(peak1_start, peak2_start)

    return(peak_overlap / peak_combined_length)
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
#' @param int1_anchor_a_start integer vector; genomic start coordinate(s)
#' of anchor A in replicate 1 interaction
#' @param int1_anchor_a_end integer vector; genomic end coordinate(s)
#' of anchor A in replicate 1 interaction
#' @param int1_anchor_b_start integer vector; genomic start coordinate(s)
#' of anchor B in replicate 1 interaction
#' @param int1_anchor_b_end integer vector; genomic end coordinate(s)
#' of anchor B in replicate 1 interaction
#' @param int2_anchor_a_start integer vector; genomic start coordinate(s)
#' of anchor A in replicate 2 interaction
#' @param int2_anchor_a_end integer vector; genomic end coordinate(s)
#' of anchor A in replicate 2 interaction
#' @param int2_anchor_b_start integer vector; genomic start coordinate(s)
#' of anchor B in replicate 2 interaction
#' @param int2_anchor_b_end integer vector; genomic end coordinate(s)
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
calculateRelativeOverlap2D <- function(int1_anchor_a_start, int1_anchor_a_end,
                                       int1_anchor_b_start, int1_anchor_b_end,
                                       int2_anchor_a_start, int2_anchor_a_end,
                                       int2_anchor_b_start, int2_anchor_b_end) {
    anchor_a_overlap <- pmin(int1_anchor_a_end, int2_anchor_a_end) -
        pmax(int1_anchor_a_start, int2_anchor_a_start)
    anchor_b_overlap <- pmin(int1_anchor_b_end, int2_anchor_b_end) -
        pmax(int1_anchor_b_start, int2_anchor_b_start)

    anchor_a_combined_length <- pmax(int1_anchor_a_end, int2_anchor_a_end) -
        pmin(int1_anchor_a_start, int2_anchor_a_start)
    anchor_b_combined_length <- pmax(int1_anchor_b_end, int2_anchor_b_end) -
        pmin(int1_anchor_b_start, int2_anchor_b_start)

    return((anchor_a_overlap + anchor_b_overlap) /
               (anchor_a_combined_length + anchor_b_combined_length))
}

#' @title Establish m:n Mapping Between Peaks from Replicate 1 and 2
#'
#' @description
#' This method returns all overlapping interactions between two replicates.
#' For each pair of overlapping interactions, the
#' \emph{ambiguity resolution value} (ARV) is calculated, which helps to reduce
#' the m:n mapping to a 1:1 mapping. The semantics of the ARV depend on the
#' specified \code{ambiguity_resolution_method}, but in general interaction
#' pairs with lower ARVs have priority over interaction pairs with higher ARVs
#' when the bijective mapping is established.
#'
#' @param rep1_df data frame of observations (i.e., genomic peaks) of
#' replicate 1, with at least the following columns (position of columns
#' matter, column names are irrelevant):
#' \tabular{rll}{
#'   column 1:  \tab \code{chr} \tab character; genomic location of peak -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2:  \tab \code{start} \tab integer; genomic location of peak -
#'   start coordinate\cr
#'   column 3:  \tab \code{end} \tab integer; genomic location of peak -
#'   end coordinate\cr
#'   column 4:  \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the interactions
#' }
#' @param rep2_df data frame of observations (i.e., genomic peaks) of
#' replicate 2, with the following columns (position of columns
#' matter, column names are irrelevant):
#' \tabular{rll}{
#'   column 1:  \tab \code{chr} \tab character; genomic location of peak -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2:  \tab \code{start} \tab integer; genomic location of peak -
#'   start coordinate\cr
#'   column 3:  \tab \code{end} \tab integer; genomic location of peak -
#'   end coordinate\cr
#'   column 4:  \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the interactions
#' }
#' @param ambiguity_resolution_method defines how ambiguous assignments
#' (when one interaction in replicate 1 overlaps with multiple interactions in
#' replicate 2 or vice versa)
#' are resolved. Available methods:
#' \tabular{rl}{
#'   \code{"value"} \tab interactions are prioritized by ascending or descending
#'   \code{value} column (see \code{sorting_direction}), e.g., if two
#'   interactions in replicate 1 overlap with one interaction in replicate 2,
#'   the interaction from replicate 1 is chosen which has a lower (if
#'   \code{sorting_direction} is \code{"ascending"}) or higher (if
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
#' \tabular{rll}{
#'   column 1:  \tab \code{rep1_idx} \tab index of interaction in replicate 1
#'   (i.e., row index in \code{rep1_df})\cr
#'   column 2:  \tab \code{rep2_idx} \tab index of interaction in replicate 2
#'   (i.e., row index in \code{rep2_df})\cr
#'   column 3:  \tab \code{arv} \tab ambiguity resolution value used turn
#'   m:n mapping into 1:1 mapping. Interaction pairs with lower \code{arv}
#'   are prioritized.
#' }
#'
#' @examples
#' rep1_df <- idr2d:::chipseq$rep1_df
#' rep1_df$value <- preprocess(rep1_df$value, "log_additive_inverse")
#'
#' rep2_df <- idr2d:::chipseq$rep2_df
#' rep2_df$value <- preprocess(rep2_df$value, "log_additive_inverse")
#'
#' # shuffle to break preexisting order
#' rep1_df <- rep1_df[sample.int(nrow(rep1_df)), ]
#' rep2_df <- rep2_df[sample.int(nrow(rep2_df)), ]
#'
#' # sort by value column
#' rep1_df <- dplyr::arrange(rep1_df, value)
#' rep2_df <- dplyr::arrange(rep2_df, value)
#'
#' pairs_df <- overlap1D(rep1_df, rep2_df)
#'
#' @export
overlap1D <- function(rep1_df, rep2_df,
                      ambiguity_resolution_method = c("overlap",
                                                      "midpoint",
                                                      "value"),
                      max_gap = 0L) {
    # argument handling
    ambiguity_resolution_method <- match.arg(ambiguity_resolution_method,
                                             choices = c("overlap",
                                                         "midpoint",
                                                         "value"))

    if (nrow(rep1_df) > 0 && nrow(rep2_df) > 0) {
        overlaps_df <- anchorOverlap(
            data.frame(chr = rep1_df[, 1],
                       start = rep1_df[, 2],
                       end = rep1_df[, 3]),
            data.frame(chr = rep2_df[, 1],
                       start = rep2_df[, 2],
                       end = rep2_df[, 3]),
            max_gap = max_gap
        )

        idx_rep1 <- overlaps_df$rep1_idx
        idx_rep2 <- overlaps_df$rep2_idx

        # arv is ambiguity resolution value
        # the lower, the better
        if (ambiguity_resolution_method == "value") {
            arv <- rep1_df[idx_rep1, 4] + rep2_df[idx_rep2, 4]
            arv <- (-1) * arv
        } else if (ambiguity_resolution_method == "overlap") {
            arv <- calculateRelativeOverlap1D(
                rep1_df[idx_rep1, 2], rep1_df[idx_rep1, 3],
                rep2_df[idx_rep2, 2], rep2_df[idx_rep2, 3]
            )
            arv <- (-1) * arv
        } else if (ambiguity_resolution_method == "midpoint") {
            arv <- calculateMidpointDistance1D(
                rep1_df[idx_rep1, 2], rep1_df[idx_rep1, 3],
                rep2_df[idx_rep2, 2], rep2_df[idx_rep2, 3]
            )
        } else {
            stop(paste0("unknown ambiguity resolution method: ",
                        ambiguity_resolution_method))
        }

        pairs_df <- data.frame(rep1_idx = idx_rep1,
                               rep2_idx = idx_rep2,
                               arv = arv)
    } else {
        pairs_df <- data.frame()
    }

    return(pairs_df)
}


#' @title Establish m:n mapping between interactions from replicate 1 and 2
#'
#' @description
#' This method returns all overlapping interactions between two replicates.
#' For each pair of overlapping interactions, the
#' \emph{ambiguity resolution value} (ARV) is calculated, which helps to reduce
#' the m:n mapping to a 1:1 mapping. The semantics of the ARV depend on the
#' specified \code{ambiguity_resolution_method}, but in general interaction
#' pairs with lower ARVs have priority over interaction pairs with higher ARVs
#' when the bijective mapping is established.
#'
#' @param rep1_df data frame of observations (i.e., genomic interactions) of
#' replicate 1, with at least the following columns (position of columns
#' matter, column names are irrelevant):
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
#'   column 7: \tab \code{value} \tab numeric; p-value, FDR, or heuristic
#'   used to rank the interactions
#' }
#' @param rep2_df data frame of observations (i.e., genomic interactions) of
#' replicate 2, with the following columns (position of columns
#' matter, column names are irrelevant):
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
#'   to rank the interactions
#' }
#' @param ambiguity_resolution_method defines how ambiguous assignments
#' (when one interaction in replicate 1 overlaps with multiple interactions in
#' replicate 2 or vice versa)
#' are resolved. Available methods:
#' \tabular{rl}{
#'   \code{"value"} \tab interactions are prioritized by ascending or descending
#'   \code{value} column (see \code{sorting_direction}), e.g., if two
#'   interactions in replicate 1 overlap with one interaction in replicate 2,
#'   the interaction from replicate 1 is chosen which has a lower (if
#'   \code{sorting_direction} is \code{"ascending"}) or higher (if
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
#' \tabular{rll}{
#'   column 1: \tab \code{rep1_idx} \tab index of interaction in replicate 1
#'   (i.e., row index in \code{rep1_df})\cr
#'   column 2: \tab \code{rep2_idx} \tab index of interaction in replicate 2
#'   (i.e., row index in \code{rep2_df})\cr
#'   column 3: \tab \code{arv} \tab ambiguity resolution value used turn
#'   m:n mapping into 1:1 mapping. Interaction pairs with lower \code{arv}
#'   are prioritized.
#' }
#'
#' @examples
#' rep1_df <- idr2d:::chiapet$rep1_df
#' rep1_df$fdr <- preprocess(rep1_df$fdr, "log_additive_inverse")
#'
#' rep2_df <- idr2d:::chiapet$rep2_df
#' rep2_df$fdr <- preprocess(rep2_df$fdr, "log_additive_inverse")
#'
#' # shuffle to break preexisting order
#' rep1_df <- rep1_df[sample.int(nrow(rep1_df)), ]
#' rep2_df <- rep2_df[sample.int(nrow(rep2_df)), ]
#'
#' # sort by value column
#' rep1_df <- dplyr::arrange(rep1_df, rep1_df$fdr)
#' rep2_df <- dplyr::arrange(rep2_df, rep2_df$fdr)
#'
#' pairs_df <- overlap2D(rep1_df, rep2_df)
#'
#' @export
overlap2D <- function(rep1_df, rep2_df,
                      ambiguity_resolution_method = c("overlap",
                                                      "midpoint",
                                                      "value"),
                      max_gap = 0L) {
    # argument handling
    ambiguity_resolution_method <- match.arg(ambiguity_resolution_method,
                                             choices = c("overlap",
                                                         "midpoint",
                                                         "value"))
    if (nrow(rep1_df) > 0 && nrow(rep2_df) > 0) {
        overlaps_anchors_a_df <- anchorOverlap(
            data.frame(chr = rep1_df[, 1],
                       start = rep1_df[, 2],
                       end = rep1_df[, 3]),
            data.frame(chr = rep2_df[, 1],
                       start = rep2_df[, 2],
                       end = rep2_df[, 3]),
            max_gap = max_gap
        )
        overlaps_anchors_b_df <- anchorOverlap(
            data.frame(chr = rep1_df[, 4],
                       start = rep1_df[, 5],
                       end = rep1_df[, 6]),
            data.frame(chr = rep2_df[, 4],
                       start = rep2_df[, 5],
                       end = rep2_df[, 6]),
            max_gap = max_gap
        )

        a <- paste0(overlaps_anchors_a_df$rep1_idx, "-",
                    overlaps_anchors_a_df$rep2_idx)
        b <- paste0(overlaps_anchors_b_df$rep1_idx, "-",
                    overlaps_anchors_b_df$rep2_idx)

        replicates <- intersect(a, b)

        idx_rep1 <- unlist(lapply(strsplit(replicates, "-", fixed = TRUE),
                                  function(interaction) {
                                      return(as.integer(interaction[1]))
                                  }))

        idx_rep2 <- unlist(lapply(strsplit(replicates, "-", fixed = TRUE),
                                  function(interaction) {
                                      return(as.integer(interaction[2]))
                                  }))

        # arv is ambiguity resolution value
        # the lower, the better
        if (ambiguity_resolution_method == "value") {
            arv <- rep1_df[idx_rep1, 7] + rep2_df[idx_rep2, 7]
            arv <- (-1) * arv
        } else if (ambiguity_resolution_method == "overlap") {
            arv <- calculateRelativeOverlap2D(
                rep1_df[idx_rep1, 2], rep1_df[idx_rep1, 3],
                rep1_df[idx_rep1, 5], rep1_df[idx_rep1, 6],
                rep2_df[idx_rep2, 2], rep2_df[idx_rep2, 3],
                rep2_df[idx_rep2, 5], rep2_df[idx_rep2, 6]
            )
            arv <- (-1) * arv
        } else if (ambiguity_resolution_method == "midpoint") {
            arv <- calculateMidpointDistance2D(
                rep1_df[idx_rep1, 2], rep1_df[idx_rep1, 3],
                rep1_df[idx_rep1, 5], rep1_df[idx_rep1, 6],
                rep2_df[idx_rep2, 2], rep2_df[idx_rep2, 3],
                rep2_df[idx_rep2, 5], rep2_df[idx_rep2, 6]
            )
        } else {
            stop(paste0("unknown ambiguity resolution method: ",
                        ambiguity_resolution_method))
        }

        pairs_df <- data.frame(rep1_idx = idx_rep1,
                               rep2_idx = idx_rep2,
                               arv = arv)
    } else {
        pairs_df <- data.frame()
    }

    return(pairs_df)
}


#' @title Removes Peaks on Non-standard Chromosomes
#' @param x data frame of genomic peaks, with the following columns
#' (position of columns matter, column names are irrelevant):
#' \tabular{rll}{
#'   column 1: \tab \code{chr} \tab character; genomic location of peak -
#'   chromosome (e.g., \code{"chr3"})\cr
#'   column 2: \tab \code{start} \tab integer; genomic location of peak  -
#'   start coordinate\cr
#'   column 3: \tab \code{end} \tab integer; genomic location of peak -
#'   end coordinate\cr
#'   column 4: \tab \code{value} \tab numeric; p-value, FDR, or heuristic used
#'   to rank the peaks
#' }
#'
#' @return \code{x} without non-standard chromosomes.
#'
#' @examples
#' rep1_df <- removeNonstandardChromosomes1D(idr2d:::chipseq$rep1_df)
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @export
removeNonstandardChromosomes1D <- function(x) {
    x <- GenomicRanges::GRanges(x[, 1],
                                IRanges::IRanges(x[, 2], x[, 3]),
                                value = x[, 4])
    x <- GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = "coarse")
    x_df <- as.data.frame(x)
    x_df$width <- NULL
    x_df$strand <- NULL
    colnames(x_df) <- c("chr", "start", "end", "value")
    return(x_df)
}

#' @title Removes Interactions on Non-standard Chromosomes
#' @param x data frame of genomic interactions, with the following columns
#' (position of columns matter, column names are irrelevant):
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
#'   to rank the interactions
#' }
#'
#' @examples
#' rep1_df <- removeNonstandardChromosomes2D(idr2d:::chiapet$rep1_df)
#'
#' @return \code{x} without non-standard chromosomes.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom dplyr inner_join
#' @importFrom dplyr select
#' @export
removeNonstandardChromosomes2D <- function(x) {
    # avoid CRAN warnings
    chr_a <- start_a <- end_a <- chr_b <- start_b <- end_b <- value <- NULL

    x$idx <- seq_len(nrow(x))
    anchor_a <- GenomicRanges::GRanges(x[, 1],
                                       IRanges::IRanges(x[, 2], x[, 3]),
                                       value = x[, 7],
                                       idx = x$idx)
    anchor_a <- GenomeInfoDb::keepStandardChromosomes(anchor_a,
                                                      pruning.mode = "coarse")
    anchor_a_df <- as.data.frame(anchor_a)
    anchor_a_df$width <- NULL
    anchor_a_df$strand <- NULL
    colnames(anchor_a_df) <- c("chr_a", "start_a", "end_a", "value", "idx")

    anchor_b <- GenomicRanges::GRanges(x[, 4],
                                       IRanges::IRanges(x[, 5], x[, 6]),
                                       idx = x$idx)
    anchor_b <- GenomeInfoDb::keepStandardChromosomes(anchor_b,
                                                      pruning.mode = "coarse")
    anchor_b_df <- as.data.frame(anchor_b)
    anchor_b_df$width <- NULL
    anchor_b_df$strand <- NULL
    colnames(anchor_b_df) <- c("chr_b", "start_b", "end_b", "idx")

    x_df <- dplyr::inner_join(anchor_a_df, anchor_b_df, by = "idx")
    x_df <- dplyr::select(x_df, chr_a, start_a, end_a,
                          chr_b, start_b, end_b, value)

    return(x_df)
}
