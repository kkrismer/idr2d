context("visualization functions")
library(idr2d)

test_that("pretty_print_block_size", {
    expect_equal(idr2d:::pretty_print_block_size(1), "1 bp")
    expect_equal(idr2d:::pretty_print_block_size(20000), "20 kbp")
    expect_equal(idr2d:::pretty_print_block_size(2000000), "2 Mbp")
    expect_equal(idr2d:::pretty_print_block_size(2234332301), "2.2 Gbp")
    expect_equal(idr2d:::pretty_print_block_size(2934332301), "2.9 Gbp")
    expect_equal(idr2d:::pretty_print_block_size(2994332301), "3 Gbp")
})
