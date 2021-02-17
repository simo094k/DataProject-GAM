context('Test functions in fix_waveforms')

test_that("fix_time_shifts work", {
  signals <- fix_time_shifts(sample_record$vital$tracks, show_diagnostics = FALSE)
  expect_lt(max(diff(signals$Intellivue$ECG_II$time)), 0.0021)
  expect_lt(max(diff(signals$Intellivue$ABP$time)), 0.0081)
})

test_that("fix_time_shifts can skip", {
  signals_break <- fix_time_shifts(sample_record$vital$tracks, skip_shift = 1, show_diagnostics = FALSE)
  expect_gt(max(diff(signals_break$Intellivue$ECG_II$time)), 1.6)
})

test_that("fix_repeating_values runs", {
    sv_short <- fix_repeating_values(sample_record$vital$tracks$EV1000$SV)
    expect_lt(max(diff(sv_short$time)), 25)
    expect_gt(min(diff(sv_short$time)), 5)
    expect_equal(nrow(sv_short), 100)
})
