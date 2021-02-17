#' Plot a signal with dygraphs
#'
#' @param signal Signal dataframe
#' @param ... Tidyselect argument to select signal columns to plot. If empty, plot everything.
#' @param time_col Column to use for x-axis
#' @param main Graph title
#' @param height height of dygraph HTML element
#' @param width Width of dygraph HTML element
#' @param maxlen maximum number of samples.
#' The signal is down sampled to this number.
#' @export
dygraph_signal <- function(signal, ..., time_col = 1, main = NULL, height = NULL, width = "100%", maxlen = 200e3) {
    if (is.null(main)) main <- deparse(substitute(signal)) # Convert argument name to string


    if (! is.data.frame(signal)){
        assertthat::assert_that(is.numeric(signal))
        assertthat::assert_that(rlang::dots_n(...) == 0, msg = 'For vector signals, `...` must be empty')
        signal <- data.frame(index = 1:length(signal), signal = signal)
    }

    while (nrow(signal) > maxlen) {
        message(glue::glue('Signal too long (more than {maxlen} rows): Sample rate halved'))
        signal <- signal[seq(2, nrow(signal), by = 2),]
    }

    if (rlang::dots_n(...) == 0) {
        sig <- dplyr::relocate(signal, {{time_col}})
    } else {
        sig <- dplyr::select(signal, {{time_col}}, ...)
    }

    if ('POSIXt' %in% class(sig[[1]])) {
        sig <- xts::xts(sig[,-1], order.by = sig[[1]])
    }

    dygraphs::dygraph(sig, group = '1', ylab = '', xlab = 'Time', main = main,
                      height = height, width = width) %>%
        dygraphs::dyRangeSelector() %>%
        dygraphs::dyAxis('x', valueFormatter = 'function(ms){    var d = new Date(ms);
                          return Dygraph.zeropad(d.getHours()) + ":" +
                          Dygraph.zeropad(d.getMinutes()) + ":" +
                          Dygraph.zeropad(d.getSeconds()) + "." +
                          Dygraph.zeropad(d.getMilliseconds());}')
}

#' Add events to dygraph
#'
#' @param dygraph Dygraph to add events to
#' @param events Events dataframe
#' @param label_col Column containing labels
#' @param time_col Column to use for x-axis
#' @export
dygraph_events <- function(dygraph, events, label_col = 2, time_col = 1) {

    dygraphs::dyEvent(dygraph, events[[time_col]], events[[label_col]], labelLoc = "bottom")
}

#' Combine Dygraphs
#' @param ... Dygraphs to combine.
#' @export
dygraph_c <- function(...) {
    graphs <- list(...)

    htmltools::browsable(htmltools::tagList(graphs))
}
