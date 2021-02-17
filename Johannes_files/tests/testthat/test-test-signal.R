test_that("Test filter", {
  cvp_filt <- filter_signal(sample_record$vital$tracks$Intellivue$CVP, cutoff_frequency = 0.5)
  cvp_filt2 <- filter_signal(cvp_filt, signal_col = 'CVP_filt', cutoff_frequency = 0.1, type = 'high',
                            postfix = '_high', trim_ends = 1)

  expect_equal(nrow(cvp_filt2), 246375)
  expect_equal(names(cvp_filt2)[4], 'CVP_filt_high')

})

test_that("Test filter - band pass", {
    cvp_filt <- filter_signal(sample_record$vital$tracks$Intellivue$CVP[1e5:2e5,], cutoff_frequency = c(0.5, 1), type = 'pass')

    expect_lt(mean(cvp_filt$CVP_filt), 0.01)
    expect_lt(var(cvp_filt$CVP_filt), 1)
    expect_gt(var(cvp_filt$CVP_filt), 0.7)

})
