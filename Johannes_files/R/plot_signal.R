#' Plot a signal
#'
#' @param signal Waveform data frame with time and value
#' @param signal_col Column with y values
#' @param seconds Limit plot to first x seconds
#' @export
plot_signal <- function(signal, signal_col = 2, seconds = NULL) {
    if (!is.null(seconds)) {
        signal <- dplyr::filter(signal, time < (time[1] + lubridate::seconds(seconds)))
    }

    if (is.numeric(signal_col)) {
        signal_name <- names(signal)[signal_col]
    } else {
        signal_name <- signal_col
    }
    baseplot <- ggplot2::ggplot(signal, ggplot2::aes(time, .data[[signal_name]]))

    if (nrow(signal) < 1000) {
        baseplot + ggplot2::geom_step()
    } else{
        baseplot + ggplot2::geom_line()
    }
}

#' Geom to add point annotations to plot
#'
#' @param signal must be an EVENT signal
#' @param seconds Limit plot to first x seconds
#'
#' @export
geom_annotation <- function(signal, seconds = NULL) {
    if (!is.null(seconds)) {
        signal <- dplyr::filter(signal, time < (time[1] + lubridate::seconds(seconds)))
    }

    # Add space before labels to nudge from bottom of plot
    signal[[2]] <- paste0(" ", signal[[2]])

    eventcolname <- names(signal)[2]

    list(
        ggplot2::geom_vline(ggplot2::aes(xintercept = time), data = signal, linetype = 2),
        ggplot2::geom_text(ggplot2::aes(x = time, label = !!rlang::ensym(eventcolname)),
                           y = -Inf, nudge_y = 10, vjust = -0.3, hjust = 0, angle = 90,
                           data = signal)
    )
}
