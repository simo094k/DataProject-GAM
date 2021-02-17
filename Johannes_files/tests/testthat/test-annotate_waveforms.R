context('Add annotations to waveforms')

abp_subset <- head(sample_record$vital$tracks$Intellivue$ABP, 30000)

test_that('find_abp_beats finds all peaks', {
    abp_beats <- find_abp_beats(abp)
    expect_equal(nrow(abp_beats), 52)

    abp_beats2 <- find_abp_beats(abp_subset)
    expect_equal(nrow(abp_beats2), 199)

    expect_equal(head(abp_beats2$PP), c(44, 44.125, 45.6250, 43.9375, 43.3125, 43.5))
    expect_equal(tail(abp_beats2$PP), c(39.750, 33.625, 34.8125, 34.6875, 36.1875, 35.3750))

})

test_that('gen_annotation_index creats reasonable values', {
    abp_beats2 <- find_abp_beats(abp_subset)
    resp_start <- abp_subset$time[seq(20, length(abp_subset$time), by = 1e3)]
    resp_index <- gen_annotation_index(abp_beats2, time_annotation = resp_start)

    # Annotation index between
    expect_gt(max(resp_index$ann_index, na.rm = TRUE), 7.9)
    expect_lte(max(resp_index$ann_index, na.rm = TRUE), 8)
    expect_gte(min(resp_index$ann_index, na.rm = TRUE), 0)
    expect_lt(min(resp_index$ann_index, na.rm = TRUE), 0.1)
})

test_that('is_extra_systole() finds extrasystoles', {
    abp_beats2 <- find_abp_beats(abp_subset)
    abp_beats2 <- dplyr::mutate(abp_beats2, extra_systole = is_extra_systole(time))
    expect_equal(sum(abp_beats2$extra_systole), 17)
})

test_that('flag_beats works', {
    abp_beats2 <- find_abp_beats(abp_subset)
    flagged <- flag_beats(abp_beats2, max_pos_after_sys = 2)
    expect_equal(which(flagged),
                 c(62L, 85L, 86L, 87L, 88L, 95L, 96L, 97L, 98L, 99L, 100L, 101L,
                   102L, 103L, 104L, 105L, 109L, 110L, 117L, 118L, 119L, 120L, 121L,
                   147L, 175L, 176L, 177L, 178L, 179L, 180L, 181L, 182L, 183L, 184L,
                   185L, 186L, 187L, 188L, 189L, 190L, 191L, 192L, 193L, 194L, 195L,
                   196L, 197L))
})
