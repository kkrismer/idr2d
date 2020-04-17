context("Hi-C functions")
library(idr2d)

test_that("estimate_idr2d_hic", {
    futile.logger::flog.threshold(futile.logger::WARN)
    idr_results_df <- estimate_idr2d_hic(idr2d:::hic$rep1_df,
                                         idr2d:::hic$rep2_df)

    expect_equal(nrow(idr_results_df), 3081)
    expect_equal(colnames(idr_results_df), c("interaction", "value",
                                             "rep_value", "rank",
                                             "rep_rank", "idr"))
    expect_equal(sum(idr_results_df$value), 5025075)
})
