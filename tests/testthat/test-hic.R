context("Hi-C functions")

test_that("estimate_idr2d_hic", {
    futile.logger::flog.threshold(futile.logger::WARN)
    set.seed(3)
    idr_results_df <- estimate_idr2d_hic(idr2d:::hic$rep1_df,
                                         idr2d:::hic$rep2_df)

    expect_equal(nrow(idr_results_df), 3081)
    expect_equal(colnames(idr_results_df), c("interaction", "value",
                                             "rep_value", "rank",
                                             "rep_rank", "idr"))
    expect_equal(sum(idr_results_df$value), 5025075)
    expect_equal(sum(idr_results_df$idr), 7.262495,
                 tolerance = 0.2)

    summarized_results <- summary(idr_results_df)
    expect_equal(summarized_results$analysis_type, "IDR2D Hi-C")
    expect_equal(summarized_results$num_blocks, 3081)
    expect_equal(summarized_results$num_significant_blocks, 3063)
    expect_equal(summarized_results$num_highly_significant_blocks, 3022)
    expect_output(print(summary(idr_results_df)))

    idr_results_df <- estimate_idr2d_hic(idr2d:::hic$rep1_df,
                                         idr2d:::hic$rep2_df, local_idr = FALSE,
                                         combined_max_value = 10^7,
                                         min_value = -20,
                                         max_value = 10^7)
    expect_equal(sum(idr_results_df$idr), 1.091252,
                 tolerance = 0.2)
})
