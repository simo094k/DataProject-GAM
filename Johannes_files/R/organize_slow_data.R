#' Generate simple data frame of device settings
#'
#' @param dev_settings Data frame of device settings, from vscapture
#' @param min_setting_len Settings must be changed one at a time.
#' if a setting is kept for shorter than this duration (in seconds) it is considered
#' to be part of a change.
#'
#' @export
gen_device_settings_annotation <- function(dev_settings, min_setting_len = 0) {

    # Mark each change in settings, to ensure that distinct does not remove previously used settings
    res <- dplyr::mutate(dev_settings,
                         .rr_change = c(diff(tidyr::replace_na(Respiratory_rate, 0)), 0) != 0,
                         .vt_change = c(diff(tidyr::replace_na(Inspiratory_tidal_volume_L, 0)), 0) != 0,
                         .setting_num = cumsum(.rr_change + .vt_change)
                         )

    res <- dplyr::distinct(res, Inspiratory_tidal_volume_L, Respiratory_rate, .setting_num,
                               .keep_all = TRUE)



    res <- dplyr::mutate(res,
                         # Add end time to each setting
                         time_end = dplyr::lead(time, default = dev_settings$time[length(dev_settings$time)]),
                         label = glue::glue("RR: {Respiratory_rate}, Vt: {Inspiratory_tidal_volume_L}"))

    duration <- c(as.numeric(diff(res$time),units = 'secs'), Inf)

    dplyr::select(res[duration > min_setting_len, ], !dplyr::starts_with('.'))

}
