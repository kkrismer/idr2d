context("auxiliary functions")
library(idr2d)

test_that("determine_anchor_overlap", {
    rep1_df <- idr2d:::chiapet$rep1_df
    rep2_df <- idr2d:::chiapet$rep2_df

    rep1_anchor_a <- data.frame(chr = rep1_df[, 1],
                                start = rep1_df[, 2],
                                end = rep1_df[, 3])
    rep2_anchor_a <- data.frame(chr = rep2_df[, 1],
                                start = rep2_df[, 2],
                                end = rep2_df[, 3])

    anchor_a_overlap <- determine_anchor_overlap(rep1_anchor_a, rep2_anchor_a)

    expect_equal(nrow(anchor_a_overlap), 42322)
    expect_equal(colnames(anchor_a_overlap), c("rep1_idx", "rep2_idx"))
})

test_that("establish_overlap1d", {
    rep1_df <- idr2d:::chipseq$rep1_df
    rep1_df$value <- preprocess(rep1_df$value, "log_additive_inverse")

    rep2_df <- idr2d:::chipseq$rep2_df
    rep2_df$value <- preprocess(rep2_df$value, "log_additive_inverse")

    rep1_df <- rep1_df[sample.int(nrow(rep1_df)), ]
    rep2_df <- rep2_df[sample.int(nrow(rep2_df)), ]

    rep1_df <- dplyr::arrange(rep1_df, value)
    rep2_df <- dplyr::arrange(rep2_df, value)

    pairs_df <- establish_overlap1d(rep1_df, rep2_df)

    expect_equal(nrow(pairs_df), 505)
    expect_equal(colnames(pairs_df), c("rep1_idx", "rep2_idx", "arv"))
    expect_equal(sum(pairs_df$arv), -298.2899,
                 tolerance = 0.00002)

    pairs_df <- establish_overlap1d(rep1_df, rep2_df,
                                    ambiguity_resolution_method = "value")
    expect_equal(sum(pairs_df$arv), 3204.912,
                 tolerance = 0.00002)

    pairs_df <- establish_overlap1d(rep1_df, rep2_df,
                                    ambiguity_resolution_method = "midpoint")
    expect_equal(sum(pairs_df$arv), 55044,
                 tolerance = 0.00002)

    expect_error(establish_overlap1d(rep1_df, rep2_df,
                                     ambiguity_resolution_method = "xx"))


    pairs_df <- establish_overlap1d(rep1_df, data.frame(chr = character(0),
                                                        start = integer(0),
                                                        end = integer(0),
                                                        value = numeric(0)))
    expect_equal(nrow(pairs_df), 0)
})

test_that("establish_overlap2d", {
    rep1_df <- idr2d:::chiapet$rep1_df
    rep1_df$fdr <- preprocess(rep1_df$fdr, "log_additive_inverse")

    rep2_df <- idr2d:::chiapet$rep2_df
    rep2_df$fdr <- preprocess(rep2_df$fdr, "log_additive_inverse")

    # shuffle to break preexisting order
    rep1_df <- rep1_df[sample.int(nrow(rep1_df)), ]
    rep2_df <- rep2_df[sample.int(nrow(rep2_df)), ]

    # sort by value column
    rep1_df <- dplyr::arrange(rep1_df, rep1_df$fdr)
    rep2_df <- dplyr::arrange(rep2_df, rep2_df$fdr)

    pairs_df <- establish_overlap2d(rep1_df, rep2_df)

    expect_equal(nrow(pairs_df), 5979)
    expect_equal(colnames(pairs_df), c("rep1_idx", "rep2_idx", "arv"))
    expect_equal(sum(pairs_df$arv), -5298.153,
                 tolerance = 0.00002)


    pairs_df <- establish_overlap2d(rep1_df, rep2_df,
                                    ambiguity_resolution_method = "value")
    expect_equal(sum(pairs_df$arv), -34394.35,
                 tolerance = 0.00002)

    pairs_df <- establish_overlap2d(rep1_df, rep2_df,
                                    ambiguity_resolution_method = "midpoint")
    expect_equal(sum(pairs_df$arv), 1464651,
                 tolerance = 0.00002)

    expect_error(establish_overlap2d(rep1_df, rep2_df,
                                     ambiguity_resolution_method = "xx"))


    pairs_df <- establish_overlap2d(rep1_df, data.frame())
    expect_equal(nrow(pairs_df), 0)
})

test_that("remove_nonstandard_chromosomes1d", {
    df <- data.frame(chr = c("chr1", "chr2", "chr1_gl456210_random"),
                     start = c(1, 2, 3),
                     end = c(3, 4, 5),
                     value = c(10, 20, 30))
    df <- remove_nonstandard_chromosomes1d(df)

    expect_equal(nrow(df), 2)
    expect_equal(df$chr, factor(c("chr1", "chr2")))
})

test_that("remove_nonstandard_chromosomes2d", {
    df <- data.frame(chr_a = c("chr1", "chr2", "chr1_gl456210_random"),
                     start_a = c(1, 2, 3),
                     end_a = c(3, 4, 5),
                     chr_b = c("chr1", "chr2", "chr1"),
                     start_b = c(1, 2, 3),
                     end_b = c(3, 4, 5),
                     value = c(10, 20, 30))
    df <- remove_nonstandard_chromosomes2d(df)

    expect_equal(nrow(df), 2)
    expect_equal(df$chr_a, factor(c("chr1", "chr2")))
})
