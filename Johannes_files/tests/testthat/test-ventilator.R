test_that("add volume works", {
    sample_vscapture_rt$RT
    w_volume <- add_volume_signal(sample_vscapture_rt$RT)

    expect_equal(mean(tail(w_volume$volume, 1e5)), 171.7, tolerance = 0.01)
})
