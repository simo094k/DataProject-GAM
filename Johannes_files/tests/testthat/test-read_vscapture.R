context('Load and test VSCapture')

vscapture_path <- system.file('extdata','vs_capture_test', package = 'waveformtools')


test_that("slow data loads", {
  slow_test <- read_vscapture_slow(vscapture_path)
  expect_length(slow_test$TextMessages, 12)
  expect_equal(nrow(slow_test$MeasuredCP1), 26)

})

test_that("rt data loads", {
  expect_warning(rt_test <- read_vscapture_rt(vscapture_path), 'calculated time')
  expect_length(rt_test$RT$flow, 208581)
  expect_length(rt_test, 2)

})

test_that("organize_slow_data works", {
  simple_settings <- gen_device_settings_annotation(sample_record$vscapture$DeviceSettings,
                                                    min_setting_len = 5)
  expect_length(simple_settings$time, 12)

})

test_that("volume calculation works", {
  rt_vol <- add_volume_signal(sample_vscapture_rt$RT)
  expect_length(rt_vol, 6)
  expect_lt(quantile(rt_vol$volume, .90), 450)
  expect_gt(quantile(rt_vol$volume, .90), 440)
})
