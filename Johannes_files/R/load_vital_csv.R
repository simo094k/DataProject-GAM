signal_trkid_pattern <- '.*_(\\d+)\\.csv.*'

#' Read Header CSV (tracks.csv) from parse_vital
#'
#' @param path Path to folder containing the files exported by parse_vital
#' @importFrom rlang .data
#' @export
read_vital_header <- function(path) {
    stopifnot(length(path) == 1)

    res <- readr::read_csv(paste0(path, '/tracks.csv'),
                           col_types = readr::cols(
                               trkid = readr::col_integer(),
                               name = readr::col_character(),
                               unit = readr::col_character(),
                               srate = readr::col_double(),
                               devname = readr::col_character()
                           )
    )


    res
}

#' Load csv files created by custom python module.
#'
#' @param path Path to .csv file
#' @param var_name name given to the values column (e.g. CVP)
#' @param tz time zone
#' @export
read_vital_csv <- function(path, var_name = 'value', tz = 'UTC') {
    readr::read_csv(path,
                    col_names = c('time', var_name),
                    col_types = list(readr::col_datetime(),
                                     readr::col_guess()),
                    locale = readr::locale(tz = tz))
}

#' Read Vital
#'
#' Loads a folder of gz compressed csv files, exported by parse_vital
#'
#' @param path Path to folder containing the files exported by parse_vital
#' @param tz Time zone
#' @param nested_list Create a nested list of tracks inside a list of devices.
#' Necessary to deal with duplicate track names between devices.
#' @param tracks_only Do not include header in returned list.
#' If FALSE, an error is given if there are any duplicate track names.
#' @param recalc_time recalculate time using first value and sample rate.
#' (Should probably not be used).
#' @return A (nested) list of tracks.
#' @examples
#' \dontrun{
#' read_vital(folder, tz = 'CET')
#' }
#' @importFrom rlang .data
#' @export
read_vital <- function(path, tz = 'UTC', nested_list = TRUE, tracks_only = FALSE, recalc_time = FALSE) {
    header <- read_vital_header(path)

    files_in_folder <- list.files(path)
    csv_in_folder <- grep('\\.csv', files_in_folder, value = TRUE)

    # Find trackids
    trackids <- gsub(signal_trkid_pattern, '\\1', csv_in_folder)

    files_to_load <- csv_in_folder[csv_in_folder != 'tracks.csv']
    paths_to_load <- file.path(path, files_to_load)

    # Give a name to the EVENT 'device'
    header$devname[header$name == 'EVENT'] <- 'VITAL'

    header$path <- paths_to_load[match(header$trkid, trackids)]


    read_vital_csv_by_row <- function(row, tz) {
        res <- read_vital_csv(row$path, var_name = row$name, tz = tz)

        if(row$srate == 1) {
            samp_len <- as.numeric((res$time[11] - res$time[1]), units = 'secs') / 10
            srate <- 1/samp_len
        }
        else {
            samp_len <- 1 / row$srate
            srate <- row$srate
            }

        if (recalc_time) {
            new_time <- res$time[1] + (0:(length(res$time) - 1)) * samp_len

            print(glue::glue('{row$name} \nMax difference between recalculated time and org time:  {max(res$time - new_time)}'))
            res$time <- new_time
        }

        attr(res, 'signal.unit') <- row$unit
        attr(res, 'signal.samplerate') <- srate

        res
    }



    if (nested_list) {

        # Split dataframe into a nested list of tracks inside device
        header_list_device <- split(header, header$devname)
        header_list_device_track <- lapply(header_list_device, function(x) split(x, x$name))

        # Nested lapply, to apply select_vital_loader to second layer (tracks)
        tracks <- lapply(header_list_device_track, lapply, read_vital_csv_by_row, tz = tz)
    }

    else {
        stopifnot(!anyDuplicated(header$name))

        header_list_track <- split(header, header$name)

        tracks <- lapply(header_list_track, read_vital_csv_by_row, tz = tz)

    }

    if (tracks_only) return(tracks)

    list(header = header, tracks = tracks)
}

