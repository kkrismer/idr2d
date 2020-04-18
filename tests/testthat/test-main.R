context("main functions")
library(idr2d)

test_that("establish_bijection1d", {
    rep1_df <- idr2d:::chipseq$rep1_df
    rep1_df$value <- preprocess(rep1_df$value, "log")

    rep2_df <- idr2d:::chipseq$rep2_df
    rep2_df$value <- preprocess(rep2_df$value, "log")

    mapping <- establish_bijection1d(rep1_df, rep2_df)
    expect_equal(length(mapping), 2)
    expect_equal(names(mapping), c("rep1_df","rep2_df"))

    rep1_df <- mapping$rep1_df
    rep2_df <- mapping$rep2_df
    expect_equal(nrow(rep1_df), 20978)
    expect_equal(colnames(rep1_df), c("chr", "start", "end", "value",
                                      "rep_value", "rank", "rep_rank",
                                      "idx", "rep_idx"))
    expect_equal(nrow(rep2_df), 20979)
})

test_that("establish_bijection2d", {
    rep1_df <- idr2d:::chiapet$rep1_df
    rep1_df$fdr <- preprocess(rep1_df$fdr, "log_additive_inverse")

    rep2_df <- idr2d:::chiapet$rep2_df
    rep2_df$fdr <- preprocess(rep2_df$fdr, "log_additive_inverse")

    mapping <- establish_bijection2d(rep1_df, rep2_df)

    expect_equal(length(mapping), 2)
    expect_equal(names(mapping), c("rep1_df","rep2_df"))

    rep1_df <- mapping$rep1_df
    rep2_df <- mapping$rep2_df
    expect_equal(nrow(rep1_df), 9928)
    expect_equal(colnames(rep1_df), c("chr_a", "start_a", "end_a",
                                      "chr_b", "start_b", "end_b",
                                      "value", "rep_value", "rank",
                                      "rep_rank", "idx", "rep_idx"))
    expect_equal(nrow(rep2_df), 10326)
})

test_that("establish_bijection", {
    rep1_df <- idr2d:::chipseq$rep1_df
    rep1_df$value <- preprocess(rep1_df$value, "log")

    rep2_df <- idr2d:::chipseq$rep2_df
    rep2_df$value <- preprocess(rep2_df$value, "log")

    mapping <- establish_bijection(rep1_df, rep2_df, analysis_type = "IDR1D")

    expect_equal(length(mapping), 2)
    expect_equal(names(mapping), c("rep1_df","rep2_df"))

    rep1_df <- mapping$rep1_df
    rep2_df <- mapping$rep2_df
    expect_equal(nrow(rep1_df), 20978)
    expect_equal(colnames(rep1_df), c("chr", "start", "end", "value",
                                      "rep_value", "rank", "rep_rank",
                                      "idx", "rep_idx"))
    expect_equal(nrow(rep2_df), 20979)

    expect_error(establish_bijection(rep1_df, rep2_df, analysis_type = "xyz"))

    mapping <- establish_bijection(rep1_df, data.frame(), analysis_type = "IDR1D")
    expect_equal(length(mapping), 2)
    expect_equal(names(mapping), c("rep1_df","rep2_df"))

    rep2_df <- mapping$rep2_df
    expect_equal(nrow(rep2_df), 0)
})

test_that("preprocess", {
    rep1_df <- idr2d:::chiapet$rep1_df

    expect_equal(sum(preprocess(rep1_df$fdr, "log_additive_inverse")), 14836.63,
                 tolerance = 0.00002)
    expect_equal(sum(preprocess(rep1_df$fdr, "identity")), 9339.174,
                 tolerance = 0.00002)
    expect_equal(sum(preprocess(rep1_df$fdr, "additive_inverse")), -9339.174,
                 tolerance = 0.00002)
    expect_equal(sum(preprocess(c(0.1, 0.54, 2.1, 0.1, 0.3, 0.4, 0.1,
                                  0.4, 22.21, 2.111, 0.331, 0.21, 21.2),
                                "multiplicative_inverse")), 49.01033,
                 tolerance = 0.00002)
    expect_equal(sum(preprocess(rep1_df$fdr, "log")), -14836.63,
                 tolerance = 0.00002)
    expect_equal(sum(preprocess(rep1_df$fdr, "log_additive_inverse")), 14836.63,
                 tolerance = 0.00002)
})

test_that("estimate_idr1d", {
    set.seed(3)
    idr_results <- estimate_idr1d(idr2d:::chipseq$rep1_df,
                                  idr2d:::chipseq$rep2_df,
                                  value_transformation = "log")

    expect_equal(length(idr_results), 3)
    expect_equal(names(idr_results), c("rep1_df", "rep2_df", "analysis_type"))


    rep1_df <- idr_results$rep1_df
    expect_equal(nrow(rep1_df), 20978)
    expect_equal(sum(rep1_df$idr, na.rm = TRUE), 350.6883,
                 tolerance = 0.1)
})

test_that("estimate_idr2d", {
    set.seed(3)
    idr_results <- estimate_idr2d(idr2d:::chiapet$rep1_df,
                                  idr2d:::chiapet$rep2_df,
                                  value_transformation = "log_additive_inverse")

    expect_equal(length(idr_results), 3)
    expect_equal(names(idr_results), c("rep1_df", "rep2_df", "analysis_type"))

    rep1_df <- idr_results$rep1_df
    expect_equal(nrow(rep1_df), 9928)
    expect_equal(sum(rep1_df$idr, na.rm = TRUE), 5252.755,
                 tolerance = 0.1)
})

test_that("estimate_idr", {
    set.seed(3)
    idr_results <- estimate_idr(idr2d:::chiapet$rep1_df,
                                idr2d:::chiapet$rep2_df,
                                analysis_type = "IDR2D",
                                value_transformation = "identity")

    expect_equal(length(idr_results), 3)
    expect_equal(names(idr_results), c("rep1_df", "rep2_df", "analysis_type"))

    rep1_df <- idr_results$rep1_df
    expect_equal(nrow(rep1_df), 9928)
    expect_equal(sum(rep1_df$idr, na.rm = TRUE), 5287.196,
                 tolerance = 0.1)

    pairs_df <- estimate_idr(idr2d:::chiapet$rep1_df,
                             data.frame(chr = character(0),
                                        start = integer(0),
                                        end = integer(0),
                                        value = numeric(0)),
                             analysis_type = "IDR1D")
    expect_equal(nrow(pairs_df$rep2_df), 0)
})
