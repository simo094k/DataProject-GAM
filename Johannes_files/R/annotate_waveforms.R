#' Create DataFrame of arterial pressure waveform annotations.
#'
#' Position and value of systolic and diastolic pressures.
#' Each row represents one beat with first a diastolic and then a systolic value.
#'
#' .noise_pos_after_sys is the sum of positive changes in ABP excluding the
#' systole and the dicrotic notch, divided by the beat length in seconds.
#' 2 mmhg/s seems like a good cutoff.
#'
#' @param data Vector of arterial blood pressure.
#' @param abp_col Index or name of column with abp data
#' @param min_PP Minimum ratio of PP/mean pressure.
#' If a beat has a lower PP, it is considered a fluctuation around the threshold.
#' The dicrotic notch is a common cause of such fluctuations.
#'
#' @param min_beat_width_s Minimum beat width (in seconds)
#' The default is 0.3 seconds
#'
#' @param win_size_avg Window size for moving avg calculation.
#' This is used to calculate the cutoff value for when a peak represents a new beat and not just noise.
#' Lower values increase flexibility. Use visualize_abp_peak_detection() to dial in.
#'
#' @param time_col Vector with the same length as `abp`, used to indicate the timing of each sample.
#' If this variable is provided, sample indexes will be substituted for this vector.
#'
#' @param show.plot Show diagnostic plot
#'
#' @param include_waveform Include waveform for each beat in output
#'
#' @param sample_rate Used if data is a vector of samples.
#'
#' @export
find_abp_beats <- function(data,
                           abp_col = 2,
                           time_col = 1,
                           min_PP = 0.20,
                           min_beat_width_s = 0.3,
                           win_size_avg = 2000,
                           show.plot = FALSE,
                           include_waveform = FALSE,
                           sample_rate = NULL) {

    assertthat::not_empty(data)


    if (is.vector(data)) {

        abp <- data
    } else {
        abp <- data[[abp_col]]
        time_vector <- data[[time_col]]
    }

    if(is.na(abp[1])) {
        return(NA)
    }

    if (is.null(sample_rate)) {
        sample_rate <- attr(data, 'signal.samplerate')
        min_beat_width <- as.integer(min_beat_width_s * 125)
        if (is.null(sample_rate)) message('Sample rate not set. Returning beat length in samples ')
    }

    #create cuttoff pressure from average pressure (movingaves is fast)
    moving_mean_p <- accelerometry::movingaves(
        #repeat avg of last measurements (* movingavg_win) to get equal lengths
        c(abp,rep(mean(tail(abp, win_size_avg)), win_size_avg-1)),
        window = win_size_avg)

    # group cycles
    # finds every crossing of cutoff pressure, and keeps only changes from high to low
    # abp < cutoff gives a series of T/F:  TTTTFFFFTTTTFFFFTT
    # diff find only changes (+ = +1):     000-000+000-000+00
    # == 1 keeps only +1 change:           000000010000000100
    # (cumsum creates groups:               000000011111111222)
    cross_index <- which(c(diff(abp < (moving_mean_p * 1.1)), 0) == 1)

    cross_groups <- splitAt(abp, cross_index, trim_ends = FALSE)

    # Split abp by diastoles
    dia_index <- purrr::map_int(cross_groups, which.min) + c(0, cross_index - 1)

    # Check that no beats are shorter than mean beat width.
    while (any(diff(dia_index) < min_beat_width)) {
        short_i <- which(diff(dia_index) < min_beat_width)[1]

        # If a shorter interval is found. Keep the index with the lowest abp (the true diastole)
        # If the shorter index is the last, just remove it.
        if (short_i == length(dia_index) || abp[dia_index[short_i+1]] < abp[dia_index[short_i]]) {
            dia_index <- dia_index[-short_i]
        }
        else dia_index <- dia_index[-(short_i + 1)]

    }

    beat_groups <- splitAt(abp, dia_index, trim_ends = FALSE) # Do not trim, since the start is needed for timing.

    res <- dplyr::tibble(beat_wave = beat_groups)
    res <- dplyr::mutate(res,
                         dia = purrr::map_dbl(beat_wave, ~.x[1]),
                         sys = purrr::map_dbl(beat_wave, ~max(.x)),
                         which.sys = purrr::map_int(beat_wave, ~which.max(.x)),
                         segment_len = purrr::map_int(beat_wave, ~length(.x)),
                         dia_pos = dplyr::lag(cumsum(segment_len), default = 1),
                         sys_pos = which.sys + dia_pos - 1, # Both are 1 indexed, so to make the global position 1 indexed subtract 1
                         PP = sys - dia)

    if(!is.null(sample_rate)) res$segment_len <- res$segment_len * 1/sample_rate

    res <- dplyr::filter(res,
                         #PP/moving_mean_p[segment_start] > min_PP,
                         dia_pos > 20, # First min cannot be among the first 20 samples
                         dia_pos < length(abp) - 100) # and last min cannot be among the last 60
                                                     # (to avoid detecting a dicrotic notch as a diastole)


    # Detect Noise
    pos_after_sys <- function(x) {
        if (is.null(sample_rate)) return(NA)
        diff_x <- diff(x)

        slopes <- rle(diff_x > 0)

        # Two positive slopes longer than 5 samples are accepted
        # (systole and dichrotic notch)
        pos_slopes_long <- which(slopes$lengths > 5 & slopes$values == TRUE)

        slopes$values[head(pos_slopes_long, 2)] <- FALSE

        sum(diff_x[inverse.rle(slopes)]) / (length(diff_x) * (1/sample_rate))
    }

    res <- dplyr::mutate(res,
                         # Wiggliness is calculated similarly to gam smooth funtion wiggliness
                         # Sum of squared 2nd derivative.
                         .noise_wiggliness = purrr::map_dbl(beat_wave, ~sum(diff(diff(.x))^2)),
                         .noise_pos_after_sys = purrr::map_dbl(beat_wave, pos_after_sys))

    if(!include_waveform) {
        res$beat_wave <- NULL
    }

    res <- dplyr::select(res,
                         dia_pos,
                         dia,
                         sys_pos,
                         sys,
                         PP,
                         beat_len = segment_len,
                         contains('.noise'),
                         contains('beat_wave'))

    if (show.plot) {
        plot(abp, type = 'l')
        points(res$dia_pos, res$dia)
        points(res$sys_pos, res$sys)
        lines(moving_mean_p * 1.1, col = 'red')
    }
    if (is.vector(data)) {
        res
    } else {

        res <- dplyr::mutate(res, time = time_vector[dia_pos],
                      time_systole = time_vector[sys_pos])
        dplyr::select(res,
               time,
               dia,
               sys,
               PP,
               beat_len,
               time_systole,
               dplyr::everything())
    }

}

#' Find extra systoles from beat times
#'
#' @param beat_times time of each beat (vector of datetimes)
#' @param threshold if a beat interval (RR) is shorter than the median of
#' the preceding 10 intervals by more than this threshold (relative), it is
#' considered and extra systole (default = 0.1, 10%).
#'
#' @export
is_extra_systole <- function(beat_times, threshold = 0.1) {
    assertthat::assert_that(length(beat_times) > 5)

    beat_interval <- diff(as.numeric(beat_times, units = 'secs'))

    es <- rep(FALSE, length(beat_times))

    for (i in 4:length(beat_times)) {
        es[i] <- beat_interval[i-1] < ((1-threshold) * median(beat_interval[max(i-11, 1): (i-2)]))
    }

    es
}

#' Flag Abnormal Beats
#'
#' @param beats_df Dataframe of beats (from `find_abp_beats()`)
#' @param max_PP_change Maximum allowed PP change (%) compared to median.
#' @param max_pos_after_sys Maximum allow positive deflection excluding 2 (systole and dicrotic notch)
#' mmHg/s
#'
#' @return a vector indicating whether a beat is noisy / abnormal
#'
#' @export
flag_beats <- function(beats_df, max_pos_after_sys = 2, max_PP_change = 0.15) {
    # Median pp of 9 beats (4 on each side).
    median_PP_9 <- dplyr::lag(runmed(beats_df$PP, k = 9), default = beats_df$PP[1])

    beats_df$PP > median_PP_9 * (1 + max_PP_change) |
        beats_df$PP < median_PP_9 * (1 - max_PP_change) |
        beats_df$.noise_pos_after_sys > max_pos_after_sys
}


#' Generate index vector for signal vector.
#'
#' Index is the time since the most recent occurrence of some annotation. e.g. QRS complex
#'
#' @param data A data frame with at least a time column
#' @param time_col Index or name of time column
#' @param time_annotation A vector of time stamps for some annotation
#' @param prefix Prefix to new columns
#' @param incl_prev_lengths Include length of previous 3 cycles as variables
#' @export
gen_annotation_index <- function(data, time_annotation, time_col = 1, prefix = 'ann',
                                 incl_prev_lengths = FALSE) {
    if (nrow(data) == 0) stop("No Data")

    if (!is.data.frame(data)) {
        time_vec <- data
        data <- dplyr::tibble(time = data)
    } else {
        time_vec <- data[[time_col]]
    }


    if ("POSIXct" %in% class(time_vec) & "POSIXct" %in% class(time_annotation)) {
        time_vec <- as.numeric(time_vec, units = 'secs')
        time_annotation <- as.numeric(time_annotation, units = 'secs')
    }


    time_annotation <- sort(time_annotation)

    # add cycle length of each annotation (e.g. resp cycle length)
    # The last cycle does not end, and therefore does not have a length
    cycle_length <- diff(time_annotation)

    if (incl_prev_lengths) {
        cycle_length_previous_1 <- dplyr::lag(cycle_length)
        cycle_length_previous_2 <- dplyr::lag(cycle_length, 2)
        cycle_length_previous_3 <- dplyr::lag(cycle_length, 3)
    }

    # Create for loop, that goes through time_vec, and finds latest annotation.
    # Time since last annotation
    ann_index  <- rep(NA, length(time_vec))
    cycle_len  <- rep(NA, length(time_vec))
    ann_n      <- rep(NA, length(time_vec))

    if (incl_prev_lengths) {
        cycle_len_prev_1 <- rep(NA, length(time_vec))
        cycle_len_prev_2 <- rep(NA, length(time_vec))
        cycle_len_prev_3 <- rep(NA, length(time_vec))
    }

    i_ann <- 1L

    for (i in 1:length(time_vec)) {
        # If annotation is after current time, check next time
        if (time_annotation[i_ann] > time_vec[i]) next
        # If the next annotation is also before the current time_vec, increment i_ann by one.
        while (i_ann < length(time_annotation) && time_annotation[i_ann + 1] < time_vec[i]) {
            i_ann <- i_ann + 1L
        }

        ann_index[i]      <- time_vec[i] - time_annotation[i_ann]
        cycle_len[i]      <- cycle_length[i_ann]
        ann_n[i]          <- i_ann

        if (incl_prev_lengths) {
            cycle_len_prev_1[i] <- cycle_length_previous_1[i_ann]
            cycle_len_prev_2[i] <- cycle_length_previous_2[i_ann]
            cycle_len_prev_3[i] <- cycle_length_previous_3[i_ann]
        }

    }

    res <- dplyr::tibble(index = ann_index,
                  n = ann_n,
                  cycle_len = cycle_len,
                  rel_index = index / cycle_len)

    if (incl_prev_lengths) {
        res <- dplyr::bind_cols(res,
                       data.frame(
                           cycle_len_prev_1,
                           cycle_len_prev_2,
                           cycle_len_prev_3
                       ))
    }

    names(res) <- paste(prefix, names(res), sep = '_')

    dplyr::bind_cols(data, res)
}

#' Find R peaks in ECG
#' @description  Wrapper for detect_rpeaks() (Pan-Tompkins algorithm from {rsleep}).
#' @param ecg_df Vector of ecg waveform with a
#' @param freq Sample frequency
#' @param ... Passed to detect_rpeaks()
#' @export
gen_r_peaks_df <- function(ecg_df, freq = NULL, ...) {
    assertthat::are_equal(length(ecg_df), 2)
    assertthat::assert_that(class(ecg_df[[1]])[1] == "POSIXct")

    if (is.null(freq)) {
        freq <- 1/mean(as.numeric(diff(ecg_df[[1]][1:100])))
    }

    r_index <- detect_rpeaks(ecg_df[[2]], freq, ...)
    dplyr::tibble(time = ecg_df$time[r_index], label = 'R peak')

}

## Slightly modified from Rsleep
detect_rpeaks <- function (signal, sRate, lowcut = 0, highcut = 15, filter_order = 1,
          integration_window = 15, refractory = 200)
{
    nyquist_freq = 0.5 * sRate
    low = lowcut/nyquist_freq
    high = highcut/nyquist_freq
    bandpass <- signal::butter(n = filter_order, W = c((lowcut/(0.5 *
                                                                    sRate)), high), type = "pass")
    signal_filt <- signal::filtfilt(bandpass, c(rep(signal[1],
                                                    sRate), signal, rep(signal[length(signal)], sRate)))
    signal_filt <- signal_filt[(sRate + 1):(length(signal_filt) -
                                                sRate)]
    signal_diff <- diff(signal_filt)
    signal_squared <- signal_diff^2
    signal_squared <- c(signal_squared[1], signal_squared)
    split_ecg <- splitAt(signal_squared, 1:(length(signal_squared) %/% 1e5) * 1e5)
    split_ecg2 <- lapply(split_ecg, function(x) {
        xc <- cladoRcpp::rcpp_convolve(x, rep(1, integration_window))
        difflen <- length(xc) - length(x)
        xc <- xc[(difflen/2 + 1):(length(xc) - (difflen/2))]
        xc
    })
    signal_conv <- unlist(split_ecg2, use.names = FALSE)
    peaks <- c(0)
    limit <- mean(signal_conv) * 3
    refractory <- sRate * (refractory/1000)
    x <- signal_conv
    for (i in c(2:(length(x) - 1))) {
        if ((x[i] > limit) && (x[i] > x[i - 1]) && (x[i] > x[i +
                                                           1]) && ((i - peaks[length(peaks)]) > refractory)) {
            peaks <- c(peaks, i)
        }
    }
    peaks <- peaks[2:length(peaks)]

    return(peaks)
}
