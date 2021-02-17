#' Load csv files with slow data created by VSCapture DrgVent
#'
#' @param path Path to .csv file
#' @param start,end start and end time for import
#'
#' @export
read_vscapture_slow_csv <- function(path, start = NULL, end = NULL) {
    # var_name = gsub(signal_name_pattern, '\\1', path)
    res <- readr::read_csv(path,
                           col_names = c('time', 'name', 'value'),
                           readr::cols(
                              time = readr::col_datetime(format = '%d-%m-%Y %H:%M:%OS'),
                              name = readr::col_character(),
                              value = readr::col_character()
                           ),
                    locale = readr::locale(tz = 'CET'))

    res <- tidyr::pivot_wider(res)

    if (!is.null(start)) {
        res <- dplyr::filter(res, time > start)
    }

    if (!is.null(end)) {
        res <- dplyr::filter(res, time < end)
    }

    res
}

#' Load all csv files with slow data created by VSCapture DrgVent
#'
#' @param path Path to folder containing vital capture .csv files
#' @param start,end start and end time for import
#'
#' @export
read_vscapture_slow <- function(path, start = NULL, end = NULL) {
    files_in_folder <- list.files(path)
    csv_in_folder <- grep('\\.csv', files_in_folder, value = TRUE)

    slow_tables <- grepl('DrgVentMedibusX', csv_in_folder)

    # Find table names
    table_names <- gsub('DrgVentMedibusX(.*)DataExport.csv', '\\1', csv_in_folder[slow_tables])

    files_to_load <- csv_in_folder[slow_tables]

    paths_to_load <- file.path(path, files_to_load)

    #Apply names, to make lapply create a named list
    names(paths_to_load) <- table_names
    lapply(paths_to_load, read_vscapture_slow_csv, start, end)
}

#' Load csv file with realtime data created by VSCapture DrgVent
#'
#' @param path Path to .csv file
#' @export
read_vscapture_rt_csv <- function(path) {
    readr::read_csv(path,
                           col_names = c('time', 'unix_ms', 'value', 'insp'),
                           col_types = list(readr::col_datetime(format = '%m/%d/%Y %H:%M:%S'),
                                            readr::col_double(),
                                            readr::col_double(),
                                            readr::col_character()),
                           na = character(),
                           locale = readr::locale(tz = 'CET'))
}

#' Load csv files with realtime data created by VSCapture DrgVent
#'
#' @param path Path to VSCapture csv folder
#' @param tz time zone
#' @param return_raw return original flow and pressure tables and new corrected times as separate entries in a list.
#' @param start,end start and end time for import
#' @param show_diagnostics show plot of expected vs recorded times (recorded times are wrong, but should be correct on average).
#' @export
read_vscapture_rt <- function(path, tz = 'CET', return_raw = FALSE, start = NULL, end = NULL,
                              show_diagnostics = TRUE) {
    path_pressure <- paste0(path, '/Airway_pressureDrgWaveExport.csv')
    path_flow <- paste0(path, '/Flow_inspiratory_expiratoryDrgWaveExport.csv')

    pressure <- read_vscapture_rt_csv(path_pressure)
    flow <- read_vscapture_rt_csv(path_flow)

    if (!is.null(start)) {
        stopifnot("POSIXct" %in% class(start))
        pressure <- dplyr::filter(pressure, time > start)
        flow <- dplyr::filter(flow, time > start)
    }

    if (!is.null(end)) {
        stopifnot("POSIXct" %in% class(end))
        pressure <- dplyr::filter(pressure, time < end)
        flow <- dplyr::filter(flow, time < end)
    }

    # Check that length and time are equal
    stopifnot(pressure$unix_ms[1] == flow$unix_ms[1])

    if (! (nrow(pressure) == nrow(flow))) {
        warning(paste('Unequal length of pressure and flow tables: nrow(pressure) - nrow(flow)) = ',
                      nrow(pressure) - nrow(flow)))

        # reduce the longer table.
        if (nrow(pressure) > nrow(flow)) pressure <- head(pressure, nrow(flow))
        else flow <- head(flow, nrow(pressure))
    }

    unix_ms <- pressure$unix_ms

    diff_unix_ms <- diff(unix_ms)

    # Each consecutive group is a series of timestamps with a mean of 20 ms between.
    # Optimally, there should only be one. More can occur if the recording is temporarily interrupted.
    consecutive_groups <- cumsum(c(0, diff_unix_ms > 1000))

    if(max(consecutive_groups) > 0) {
        message(glue::glue('Obs: {max(consecutive_groups)} break(s) in ventilator waveforms'))
    }

    gen_naive_ms <- function(time_ms) {
        if (length(time_ms) < 200) return(NULL)
        if (length(time_ms) < 500) base_n <- length(time_ms)
        else base_n <- 500


        # Set starttime to the 1% quantile of unix_ms - index * 20. This finds the true baseline,
        # when times fluctuate due to slow data.
        starttime <- quantile(time_ms[1:base_n] - (0:(base_n-1) * 20), 0.01)

        starttime + 0:(length(time_ms)-1) * 20
    }

    # Generate time vector as uniform 20 ms intervals. Restart for each consecutive group
    naive_ms <- unlist(lapply(split(unix_ms, consecutive_groups), gen_naive_ms))

    if (show_diagnostics) {

        # Check that there is no drift in the difference between naive ms and unix_ms
        offset <- unix_ms - naive_ms
        median_offset_10s <- sapply(split(offset, 1:length(offset) %/% 500), median)
        quant1_offset_10s <- sapply(split(offset, 1:length(offset) %/% 500), quantile, 0.01)

        #remove start and end
        median_offset_10s <- median_offset_10s[2:(length(median_offset_10s) - 1)]

        diagnostic_plots <- FALSE

        if (any(median_offset_10s > 100)) {
            warning('read_vscapture: The median difference between calculated time and recorded time is more than 100 ms')
            diagnostic_plots <- TRUE
        }

        if (sd(median_offset_10s) > 10) {
            warning('read_vscapture: The median difference between calculated time and recorded time is not constant')
            diagnostic_plots <- TRUE
        }

        if (diagnostic_plots) {
            plot(median_offset_10s, ylim = c(min(quant1_offset_10s), max(median_offset_10s)))
            points(quant1_offset_10s,  pch='x')
        }
    }

    time <- as.POSIXct(naive_ms/1000, origin="1970-01-01", tz = tz)

    insp <- cumsum(pressure$insp == 'InspStart')

    if (as.numeric(time[1] - pressure$time[1], units = 'secs') > 1) {
        warning('read_vscapture: The median difference between calculated time and recorded datetime off by more than a second.
                \nThis is probably due to an erroneously set timezone.')
    }

    if (return_raw) {
        return(list(pressure = pressure, flow = flow, time = time))
    }

    res <- list(
        RT   = dplyr::tibble(time, pressure = pressure$value, flow = flow$value, insp),
        insp_start = dplyr::tibble(time = time[pressure$insp == "InspStart"], label = 'Insp start'))

    attr(res$RT$pressure, 'signal.unit') <- 'mBar'
    attr(res$RT$flow, 'signal.unit') <- 'L/min'
    attr(res$RT$pressure, 'signal.samplerate') <- 50
    attr(res$RT$flow, 'signal.samplerate') <- 50

    res
}

#' Load all csv files created by VSCapture DrgVent.
#' Converts column types
#'
#' @param path Path to folder containing vital capture .csv files
#' @param start,end start and end time for import
#' @param show_diagnostics show plot of expected vs recorded times (recorded times are wrong, but should be correct on average).
#' @export
read_vscapture <- function(path, start = NULL, end = NULL, show_diagnostics = TRUE) {
    rt <- read_vscapture_rt(path, start = start, end = end, show_diagnostics = show_diagnostics)
    slow <- read_vscapture_slow(path, start = start, end = end)

    slow$DeviceSettings <- dplyr::mutate_at(slow$DeviceSettings, dplyr::vars(-time), as.numeric)
    slow$MeasuredCP1 <- dplyr::mutate_at(slow$MeasuredCP1, dplyr::vars(-time), as.numeric)
    slow$MeasuredCP2 <- dplyr::mutate_at(slow$MeasuredCP2, dplyr::vars(-time), as.numeric)

    c(rt, slow)
}
