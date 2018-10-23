
#' @importFrom dplyr arrange_
#' @importFrom futile.logger flog.warn
#' @importFrom futile.logger flog.info
sortBy <- function(df, column.name) {
    if (nrow(df) > 0) {
        # shuffle to break preexisting order
        df <- df[sample.int(nrow(df)), ]

        df <- dplyr::arrange_(df, column.name)

        futile.logger::flog.info("data sorted")
    } else {
        futile.logger::flog.warn("empty data frame")
    }
    return(df)
}

#' @title Identifies Overlapping Anchors
#'
#' @description
#' TODO
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
#' @param max.gap integer; maximum gap in nucleotides allowed between two anchors for
#' them to be considered as overlapping (defaults to \code{1000L})
#'
#' @return A data frame with the following columns:
#' \tabular{rl}{
#'   \code{"queryHits"} \tab TODO\cr
#'   \code{"subjectHits"} \tab TODO
#' }
#'
#' @examples
#' TODO
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
  combined.seq.levels <- stringr::str_sort(union(seq.levels.r1, seq.levels.r2))
  GenomeInfoDb::seqlevels(rep1.ranges) <- combined.seq.levels
  GenomeInfoDb::seqlevels(rep2.ranges) <- combined.seq.levels

  # get overlap between replicates, accept 1000 bp gap
  return(data.frame(GenomicRanges::findOverlaps(rep1.ranges, rep2.ranges,
    maxgap = max.gap
  )))
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
#' @inheritParams calculateRelativeOverlap
#'
#' @return positive integer vector; distances between interaction pairs
#'
#' @examples
#' # identical, zero distance
#' calculateMidpointDistance(100, 120, 240, 260,
#'                           100, 120, 240, 260)
#'
#' # centered, zero distance
#' calculateMidpointDistance(100, 120, 240, 260,
#'                           90, 130, 230, 270)
#'
#' # off by 10 per anchor
#' calculateMidpointDistance(100, 120, 240, 250,
#'                          110, 130, 230, 240)
#'
#' # off by 10 (anchor B only)
#' calculateMidpointDistance(100, 120, 240, 250,
#'                          90, 130, 250, 260)
#'
#' # vectorized example
#' calculateMidpointDistance(c(100, 100, 100, 100),
#'                           c(120, 120, 120, 120),
#'                           c(240, 240, 240, 240),
#'                           c(260, 260, 250, 250),
#'                           c(100, 90, 110, 90),
#'                           c(120, 130, 130, 130),
#'                           c(240, 230, 230, 250),
#'                           c(260, 270, 240, 260))
#'
#' @export
calculateMidpointDistance <- function(int1.anchor.a.start, int1.anchor.a.end,
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
#' calculateRelativeOverlap(100, 120, 240, 260,
#'                          100, 120, 240, 260)
#'
#' # 50% overlap
#' calculateRelativeOverlap(100, 120, 240, 250,
#'                          100, 110, 240, 260)
#'
#' # negative overlap
#' calculateRelativeOverlap(100, 120, 240, 250,
#'                          130, 140, 260, 280)
#'
#' # larger negative overlap
#' calculateRelativeOverlap(100, 120, 240, 250,
#'                          200, 220, 340, 350)
#'
#' # vectorized example
#' calculateRelativeOverlap(c(100, 100, 100, 100),
#'                          c(120, 120, 120, 120),
#'                          c(240, 240, 240, 240),
#'                          c(260, 250, 250, 250),
#'                          c(100, 100, 130, 200),
#'                          c(120, 110, 140, 220),
#'                          c(240, 240, 260, 340),
#'                          c(260, 260, 280, 350))
#' @export
calculateRelativeOverlap <- function(int1.anchor.a.start, int1.anchor.a.end,
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
#' @title TODO
#'
#' @description
#' TODO
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
#'   \code{"midpoint"} \tab the interaction pair is chosen which has the smallest
#'   distance between their anchor midpoints, i.e., distance from midpoint of
#'   replicate 1 interaction anchor A to midpoint of
#'   replicate 2 interaction anchor A, plus distance from midpoint of
#'   replicate 1 interaction anchor B to midpoint of
#'   replicate 2 interaction anchor B \cr
#'   \code{"expansion"} \tab if one interaction in replicate 1 overlaps with
#'   \emph{n} interactions in replicate 2, the interaction of replicate 1 is
#'   copied \emph{n} times and each copy is assigned to one of the interactions
#'   in replicate 2
#' }
#' @inheritParams anchorOverlap
#'
#' @return data frame with the following columns:
#' \tabular{rl}{
#'   \code{rep1.idx} \tab index of interaction in replicate 1 (i.e., row
#'   index in \code{rep1.df})\cr
#'   \code{rep2.idx} \tab index of interaction in replicate 2 (i.e., row
#'   index in \code{rep2.df})\cr
#'   \code{arv} ambiguity resolution value used turn m:n mapping into 1:1
#'   mapping. Interaction pairs with lower \code{arv} are prioritized.\tab
#' }
#'
#'
#' @examples
#' TODO
#'
#' @importFrom futile.logger flog.warn
#' @importFrom futile.logger flog.info
#' @export
overlap <- function(rep1.df, rep2.df,
                    ambiguity.resolution.method = c("value", "overlap",
                                                    "midpoint", "expansion"),
                    max.gap = 1000L) {
    ambiguity.resolution.method <- match.arg(ambiguity.resolution.method,
                                             choices = c("value",
                                                         "overlap",
                                                         "midpoint",
                                                         "expansion"))
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

        a <- paste0(overlaps.anchorsA.df$queryHits, "-",
                    overlaps.anchorsA.df$subjectHits)
        b <- paste0(overlaps.anchorsB.df$queryHits, "-",
                    overlaps.anchorsB.df$subjectHits)

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
            arv <- calculateRelativeOverlap(
                rep1.df[idx.rep1, 2], rep1.df[idx.rep1, 3],
                rep1.df[idx.rep1, 5], rep1.df[idx.rep1, 6],
                rep2.df[idx.rep2, 2], rep2.df[idx.rep2, 3],
                rep2.df[idx.rep2, 5], rep2.df[idx.rep2, 6]
            )
            arv <- (-1) * arv
        } else if (ambiguity.resolution.method == "midpoint") {
            arv <- calculateMidpointDistance(
                rep1.df[idx.rep1, 2], rep1.df[idx.rep1, 3],
                rep1.df[idx.rep1, 5], rep1.df[idx.rep1, 6],
                rep2.df[idx.rep2, 2], rep2.df[idx.rep2, 3],
                rep2.df[idx.rep2, 5], rep2.df[idx.rep2, 6]
            )
        } else if (ambiguity.resolution.method == "expansion") {
            # no action required
            arv <- rep(0, length(idx.rep1))
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
