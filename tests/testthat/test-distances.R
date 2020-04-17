context("Overlaps and distances between interactions")
library(idr2d)

test_that("calculate_relative_overlap2d", {
    expect_equal(calculate_relative_overlap2d(100, 120, 240, 260,
                                              100, 120, 240, 260), 1.0)
    expect_equal(calculate_relative_overlap2d(100, 120, 240, 250,
                                              100, 110, 240, 260), 0.5)
    expect_equal(calculate_relative_overlap2d(100, 120, 240, 250,
                                              130, 140, 260, 280), -0.25)
    expect_equal(calculate_relative_overlap2d(100, 120, 240, 250,
                                              200, 220, 340, 350), -0.7391304,
                 tolerance = .00002)

    expect_equal(calculate_relative_overlap2d(c(100, 100, 100, 100),
                                              c(120, 120, 120, 120),
                                              c(240, 240, 240, 240),
                                              c(260, 250, 250, 250),
                                              c(100, 100, 130, 200),
                                              c(120, 110, 140, 220),
                                              c(240, 240, 260, 340),
                                              c(260, 260, 280, 350)),
                 c(1.0, 0.5, -0.25, -0.7391304),
                 tolerance = .00002)
})

test_that("calculate_midpoint_distance2d", {
    expect_identical(calculate_midpoint_distance2d(100, 120, 240, 260,
                                                   100, 120, 240, 260), 0L)
    expect_identical(calculate_midpoint_distance2d(100, 120, 240, 260,
                                                   90, 130, 230, 270), 0L)
    expect_identical(calculate_midpoint_distance2d(100, 120, 240, 250,
                                                   110, 130, 230, 240), 20L)
    expect_identical(calculate_midpoint_distance2d(100, 120, 240, 250,
                                                   90, 130, 250, 260), 10L)
    expect_identical(calculate_midpoint_distance2d(c(100, 100, 100, 100),
                                                   c(120, 120, 120, 120),
                                                   c(240, 240, 240, 240),
                                                   c(260, 260, 250, 250),
                                                   c(100, 90, 110, 90),
                                                   c(120, 130, 130, 130),
                                                   c(240, 230, 230, 250),
                                                   c(260, 270, 240, 260)),
                     c(0L, 0L, 20L, 10L))
})

test_that("calculate_relative_overlap1d", {
    expect_equal(calculate_relative_overlap1d(100, 120, 100, 120), 1.0)
    expect_equal(calculate_relative_overlap1d(100, 120, 100, 110), 0.5)
    expect_equal(calculate_relative_overlap1d(100, 120, 130, 140), -0.25)
    expect_equal(calculate_relative_overlap1d(100, 120, 200, 220), -0.666666,
                 tolerance = .00002)

    expect_equal(calculate_relative_overlap1d(c(100, 100, 100, 100),
                                              c(120, 120, 120, 120),
                                              c(100, 100, 130, 200),
                                              c(120, 110, 140, 220)),
                 c(1.0, 0.5, -0.25, -0.666666),
                 tolerance = .00002)
})

test_that("calculate_midpoint_distance1d", {
    expect_identical(calculate_midpoint_distance1d(100, 120, 100, 120), 0L)
    expect_identical(calculate_midpoint_distance1d(0, 0, 100, 120), 110L)
    expect_identical(calculate_midpoint_distance1d(100, 120, 90, 130), 0L)
    expect_identical(calculate_midpoint_distance1d(100, 120, 110, 130), 10L)
    expect_identical(calculate_midpoint_distance1d(100, 120, 90, 200), 35L)
    expect_identical(calculate_midpoint_distance1d(c(100, 100, 100, 100),
                                                   c(120, 120, 120, 120),
                                                   c(100, 90, 110, 90),
                                                   c(120, 130, 130, 200)),
                     c(0L, 0L, 10L, 35L))
})
