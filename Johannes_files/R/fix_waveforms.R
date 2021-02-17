#' Fix Time shifts in vital waveforms
#'
#' Vital Recorder may introduce short shifts in time in waveform recordings.
#' This is probably an error, and can be removed.
#' This function looks through a waveform and for every diff higher than 1/samplerate,
#' subtracts this diff - samplerate from all signal-times that are later than the culprit.
#'
#' @param signal_list list of vital signals
#' @param skip_shift Index (vector) of shift to skip (allow in output)
#' @param show_diagnostics Show diagnostic prints and plots
#' @param diagnostic_signal Signal to use in diagnostic plots
#'
#' @export
fix_time_shifts <- function(signal_list,
                            skip_shift = NULL,
                            show_diagnostics = TRUE,
                            diagnostic_signal = 'CVP') {
    assertthat::assert_that(all(c("Intellivue", "VITAL") %in% names(signal_list)))

    base_wave_times <- signal_list$Intellivue$ECG_II$time
    base_wave_times_diff_s <- diff(as_seconds(base_wave_times))
    base_wave_samplen <- 1 / attr(signal_list$Intellivue$ECG_II, 'signal.samplerate')
    long_base_wave_diff <- base_wave_times_diff_s > base_wave_samplen * 1.1

    if (!any(long_base_wave_diff)) {
        if (show_diagnostics) message('No long time differences')
        return(signal_list)
    }

    # Find shifts to remove
    ## The break length is the difference minus the sample length
    break_len <- base_wave_times_diff_s[long_base_wave_diff] - base_wave_samplen
    break_time <- base_wave_times[which(long_base_wave_diff)]

    if (show_diagnostics) {
        if (requireNamespace('gridExtra', quietly = TRUE)) {
        message(glue::glue("{length(break_len)} breaks found in record: {paste(round(break_len, 2), collapse = ', ')} second(s)"))

        pre_plots <- lapply(break_time, function(t) {
            assertthat::assert_that(diagnostic_signal %in% names(signal_list$Intellivue))

            plot_signal(dplyr::filter(signal_list$Intellivue[[diagnostic_signal]],
                                            dplyr::between(time, t - 2, t + 2))) +
                      ggplot2::ggtitle('Before break is removed') +
                      ggplot2::geom_vline(xintercept = t)
        })
        } else {
            stop('install "gridExtra" to use show_diagnostics')
        }
    }


    trim_breaks <- function(signal) {
        # Work backwards, to ensure that a timeshift does not invalidate the position of the next break
        for (i in length(break_time): 1) {

            if (i %in% skip_shift) {
                next
            }

            signal$time[signal$time > (break_time[i] + break_len[i] / 2)] <-
                signal$time[signal$time > (break_time[i] + break_len[i] / 2)] - break_len[i]

        }

        signal
    }

    signal_list$Intellivue <- lapply(signal_list$Intellivue, trim_breaks)

    if ('EV1000' %in% names(signal_list)) {
        signal_list$EV1000 <- lapply(signal_list$EV1000, trim_breaks)
    }

    if (show_diagnostics) {
        for (i in 1:length(break_time)) {

            # Break time must be corrected if more than on break is fixed
            break_time_corrected <- break_time[i] - sum(c(0, break_len)[1:i])

            if (i %in% skip_shift) post_title <- ggplot2::ggtitle('SKIPPED: break not removed')
            else post_title <- ggplot2::ggtitle('After break is removed')

            post_plot <- plot_signal(dplyr::filter(signal_list$Intellivue[[diagnostic_signal]],
                                      dplyr::between(time, break_time_corrected - 2, break_time_corrected + 2))) +
                post_title +
                ggplot2::geom_vline(xintercept = break_time_corrected)

            gridExtra::grid.arrange(pre_plots[[i]], post_plot, nrow = 1)
            if (i < length(break_time)) readline(prompt="Press [enter] to show next diagnostic plot")
        }

    }

    signal_list
}

#' Remove repeating instances of the same sample
#'
#' A device may only return a value every 20 seconds. If data is represented with 2 sec intervals, 9/10 are just
#' repeats. This function takes a dataframe with repeating values, and tries to only return 1 observation per
#' measured value.
#'
#' eg. Vital recorder saves EV1000 samples every other second. Data is only calculated every 20 seconds.
#' This function tries to match these 20 second samples.
#'
#' @param df Signal to fix (data frame)
#' @param expected_interval The interval expected from the device.
#' @param max_difftime If the interval between two changes is higher than max_timediff,
#' this is interpreted as representing multiple true measurements.
#'
#' @return Dataframe with only relevant values (1 every 20 seconds)
#' @export
#'
fix_repeating_values <- function(df, expected_interval = 20, max_difftime = 25) {
    short_df <- dplyr::filter(df,
                              c(TRUE, # First value
                                diff(df[[2]][1:(nrow(df)-1)]) != 0, # All changes
                                TRUE)) #Last value

    while(any(as_seconds(diff(short_df$time)) > max_difftime)) {
        gaps <- which(as_seconds(diff(short_df$time)) > max_difftime)

        for (i in rev(gaps)) {
            short_df <- tibble::add_row(short_df,
                                       dplyr::mutate(short_df[i, ], time = time + expected_interval),
                                       .after = i)
        }
    }

    short_df
}
