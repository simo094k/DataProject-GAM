context('Misc Functions')

test_that("join nearest works", {
  sample_slow <- sample_record$vital$tracks$Intellivue$ECG_II[seq(1, 10000, by = 1000), ]
  sample_abp <- sample_record$vital$tracks$Intellivue$ABP[1:5000,]
  res <- join_nearest(sample_slow, sample_abp)

  # Expect unchanged
  expect_known_hash(res, '164acef222')
  expect_lt(abs(mean(with(res, time - time.y))), 0.001)

  slow_df <- data.frame(key_slow = seq(0, 100, by = 10), slow_val = "a")
  fast_df <- data.frame(key_fast = seq(1, 100, by = 3), fast_val = seq(1, 100, by = 3))

  res2 <- join_nearest(slow_df, fast_df, xkey = 'key_slow', ykey = 'key_fast')

})

test_that('Subset record works', {
  interval <- c(lubridate::ymd_hms('2020-05-27 10:24:05', tz = 'CEST'),
                lubridate::ymd_hms('2020-05-27 10:25:05', tz = 'CEST'))
  sample_rec_sub <- subset_record(sample_record, interval)

  expect_equal(nrow(sample_rec_sub$vital$tracks$Intellivue$ABP), 7500)
})
