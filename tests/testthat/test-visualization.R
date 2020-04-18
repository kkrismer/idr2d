context("visualization functions")

test_that("pretty_print_block_size", {
    expect_equal(idr2d:::pretty_print_block_size(1), "1 bp")
    expect_equal(idr2d:::pretty_print_block_size(20000), "20 kbp")
    expect_equal(idr2d:::pretty_print_block_size(2000000), "2 Mbp")
    expect_equal(idr2d:::pretty_print_block_size(2234332301), "2.2 Gbp")
    expect_equal(idr2d:::pretty_print_block_size(2934332301), "2.9 Gbp")
    expect_equal(idr2d:::pretty_print_block_size(2994332301), "3 Gbp")
})

test_that("draw_idr_distribution_histogram", {
    idr_results <- estimate_idr1d(idr2d:::chipseq$rep1_df,
                                  idr2d:::chipseq$rep2_df,
                                  value_transformation = "log")
    p <- draw_idr_distribution_histogram(idr_results$rep1_df)
    expect_equal(class(p), c("gg", "ggplot"))
})

test_that("draw_rank_idr_scatterplot", {
    idr_results <- estimate_idr1d(idr2d:::chipseq$rep1_df,
                                  idr2d:::chipseq$rep2_df,
                                  value_transformation = "log")
    p <- draw_rank_idr_scatterplot(idr_results$rep1_df)
    expect_equal(class(p), c("gg", "ggplot"))
    p <- draw_rank_idr_scatterplot(idr_results$rep1_df, log_idr = TRUE)
    expect_equal(class(p), c("gg", "ggplot"))
})

test_that("draw_value_idr_scatterplot", {
    idr_results <- estimate_idr1d(idr2d:::chipseq$rep1_df,
                                  idr2d:::chipseq$rep2_df,
                                  value_transformation = "log")
    p <- draw_value_idr_scatterplot(idr_results$rep1_df)
    expect_equal(class(p), c("gg", "ggplot"))
    p <- draw_value_idr_scatterplot(idr_results$rep1_df, log_axes = TRUE)
    expect_equal(class(p), c("gg", "ggplot"))
    p <- draw_value_idr_scatterplot(idr_results$rep1_df, log_idr = TRUE)
    expect_equal(class(p), c("gg", "ggplot"))
})

test_that("draw_hic_contact_map", {
    futile.logger::flog.threshold(futile.logger::WARN)
    idr_results_df <- estimate_idr2d_hic(idr2d:::hic$rep1_df,
                                         idr2d:::hic$rep2_df)
    p <- draw_hic_contact_map(idr_results_df, idr_cutoff = 0.05,
                              chromosome = "chr1")
    expect_equal(class(p), c("gg", "ggplot"))

    expect_error(draw_hic_contact_map(data.frame(), idr_cutoff = 0.05,
                                      chromosome = "chr1"))
    expect_error(draw_hic_contact_map(idr_results_df, idr_cutoff = 0.05,
                                      chromosome = "chrWRONG"))
    p <- draw_hic_contact_map(idr_results_df, idr_cutoff = 0.05,
                              chromosome = "chr1", log_values = FALSE)
    expect_equal(class(p), c("gg", "ggplot"))
})

