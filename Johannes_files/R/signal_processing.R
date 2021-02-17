#' High or low pass filter using 2nd order Butterworth filter
#' @param signal Vector or signal dataframe with time column
#' @param cutoff_frequency Critical frequency of filter
#' @param signal_col column with signal to filter
#' @param sample_rate Sample rate of signal. If NULL, this is calculated from the time column.
#' @param type passed to [signal::butter()]. One of "low", "high", "stop", "pass"
#' @param trim_ends Seconds to trim off each end after filtering. To remove tails towards 0
#' @param postfix Postfix to the new column
#'
#' @export
filter_signal <- function(signal, cutoff_frequency, signal_col = 2,
                          sample_rate = NULL, type = 'low', trim_ends = 0,
                          postfix = '_filt') {
    if (is.data.frame(signal)) {
        signal_vec <- signal[[signal_col]]
        signal_name <- ifelse(is.numeric(signal_col), names(signal)[signal_col], signal_col)
        if (is.null(sample_rate)) {
            sample_rate <- 1 / mean(as.numeric(diff(head(signal$time))))
        }
    }
    # To set the half amplitude cutoff, divide with sample_rate/2
    # According to https://stackoverflow.com/a/49905340/9665680
    bf <- signal::butter(2, cutoff_frequency/(sample_rate/2), type=type)
    filt <- signal::filtfilt(bf, signal_vec)



    if (is.data.frame(signal)) {
        signal[[paste0(signal_name, postfix)]]  <- filt

        # Cut of tails towards 0
        signal[round(sample_rate * trim_ends):(length(filt)- round(sample_rate * trim_ends)), ]
    } else {

        # Cut of tails towards 0
        filt[round(sample_rate * trim_ends):(length(filt)- round(sample_rate * trim_ends))]
    }
}

#' Peak detection algorithm from https://github.com/stas-g/findPeaks
#'
#' @description To find valleys, use -x
#'
#' @param x signal vector
#' @param m Number of points of both sides of the peak with a lower value.
#' 2m = width of peak
#' @export
find_peaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))
        if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
    pks <- unlist(pks)
    pks
}
