context("visualization functions")

test_that("pretty_print_block_size", {
    expect_equal(idr2d:::pretty_print_block_size(1), "1 bp")
    expect_equal(idr2d:::pretty_print_block_size(20000), "20 kbp")
    expect_equal(idr2d:::pretty_print_block_size(2000000), "2 Mbp")
    expect_equal(idr2d:::pretty_print_block_size(2234332301), "2.2 Gbp")
    expect_equal(idr2d:::pretty_print_block_size(2934332301), "2.9 Gbp")
    expect_equal(idr2d:::pretty_print_block_size(2994332301), "3 Gbp")
})

test_that("pretty_print_reads", {
    expect_equal(idr2d:::pretty_print_reads(1, TRUE), "1\nnormalized reads")
    expect_equal(idr2d:::pretty_print_reads(20000, FALSE), "20000 reads")
    expect_equal(idr2d:::pretty_print_reads(2000000, FALSE), "2 million reads")
    expect_equal(idr2d:::pretty_print_reads(2234332301, FALSE),
                 "2.2 billion reads")
    expect_equal(idr2d:::pretty_print_reads(2934332301, FALSE),
                 "2.9 billion reads")
    expect_equal(idr2d:::pretty_print_reads(2994332301, FALSE),
                 "3 billion reads")
    expect_equal(idr2d:::pretty_print_reads(29943323010, FALSE),
                 "29.9 billion reads")
    expect_equal(idr2d:::pretty_print_reads(299743323010, FALSE),
                 "299.7 billion reads")
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
    p <- draw_rank_idr_scatterplot(idr_results$rep1_df, log_idr = TRUE,
                                   max_points_shown = 250)
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
    p <- draw_value_idr_scatterplot(idr_results$rep1_df, log_idr = TRUE,
                                    max_points_shown = 250)
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

