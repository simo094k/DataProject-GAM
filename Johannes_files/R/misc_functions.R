.datatable.aware = TRUE

#' Join Nearest
#' Left join two data frames, by aligning the nearest xkey with the nearest ykey.
#'
#' @param x left data frame
#' @param y right data frame
#' @param xkey,ykey keys for x and y, Defaults to joining by time
#' @import data.table
#' @export
join_nearest <- function(x, y, xkey = 'time', ykey = 'time') {

    x$roll <- x[[xkey]]
    y$roll <- y[[ykey]]
    if (xkey == ykey) names(y)[names(y) == ykey] <- paste0(ykey, '.y')
    setDT(x)
    setDT(y)
    setkey(x, roll)
    setkey(y, roll)

    res <- dplyr::as_tibble(y[x, roll='nearest'])

    # Reorder to start with key, then x, then y
    res <- res[,c(xkey,
        names(x)[names(x) != xkey],
        names(y)[names(y) != 'roll'])]

    res[['roll']] <- NULL
    res
}

#' Split by index
#'
#' Similar to `split()` but takes a vector of indexes instead of a group vector.
#' From SO: <https://stackoverflow.com/a/19274414/1498656>
#'
#' @param x input vector
#' @param pos index vector
#' @param trim_ends exclude first and last group
#'
#' @export
splitAt <- function(x, pos, trim_ends = FALSE) {
    assertthat::assert_that(is.vector(x))
    assertthat::assert_that(is.vector(pos))
    out <- list()
    pos2 <- if (trim_ends) pos else c(1, pos, length(x)+1)

    for (i in seq_along(pos2[-1])) {
        out[[i]] <- x[pos2[i]:(pos2[i+1]-1)]
    }
    return(out)
}


#' Return seconds since epoch
#'
#' @param x time vector
#'
#' @export
as_seconds <- function(x) {
    as.numeric(x, units = 'secs')
}

#' Return seconds since start of time vector
#'
#' @param x time vector
#'
#' @export
seconds_since_start <- function(x) {
    assertthat::assert_that('POSIXt' %in% class(x))
    as.numeric(x - x[1], units = 'secs')
}


#' Format as HMS with 3 decimals
#'
#' @param time_stamp Timestamp to format
#'
#' @export
format_hms <- function(time_stamp) {
    format(time_stamp, "%H:%M:%OS3")
}

#' Subset nested signals
#'
#' Subset a record of nested signals to a specific interval
#'
#' @param rec (nested) list of signals
#' @param starttime start of interval (or interval)
#' @param endtime end of interval
#' @param relative_time include time since subset start as time_s
#' @param absolute_time include original absolute time.
#' If absolute time is not included, relative time overwrites the time column
#'
#' @export
subset_record <- function(rec, starttime, endtime = NULL, relative_time = TRUE, absolute_time = TRUE) {

    assertthat::assert_that(relative_time | absolute_time)

    if (length(starttime) > 1) {
        assertthat::assert_that(is.null(endtime))
        assertthat::assert_that(length(starttime) == 2)

        endtime <- starttime[2]
        starttime <- starttime[1]
    }

    subset_func <- function(df, starttime, endtime, .xname) {
        if('time' %in% names(df)) {
            sub <- dplyr::filter(df, dplyr::between(time, starttime, endtime))

            if (relative_time & absolute_time) dplyr::mutate(sub, time_s = as_seconds(time - starttime))
            else if (relative_time & !absolute_time) dplyr::mutate(sub, time = as_seconds(time - starttime))
            else sub
        } else {
            message(paste(.xname, 'has no time column'))
            df
        }

    }

    rrapply::rrapply(
                rec,
                f = subset_func,
                classes = 'data.frame',
                how = 'replace',
                starttime = starttime,
                endtime = endtime
            )
    }


#' Shift record time
#'
#' Shift all time variables
#'
#' @param rec (nested) list of signals
#' @param seconds Time interval to shift by
#'
#' @export
shift_record_time <- function(rec, seconds) {

    shift_func <- function(df, time, .xname) {
        df_classes <- purrr::map_chr(df, ~class(.x)[1])
        if("POSIXct" %in% df_classes) {

            dplyr::mutate_at(df, which(df_classes == 'POSIXct'),
                             ~. + seconds)

        } else {
            message(paste(.xname, 'has no time columns'))
            df
        }

    }

    rrapply::rrapply(
        rec,
        f = shift_func,
        classes = 'data.frame',
        how = 'replace',
        time = seconds
    )
}

#' Anonymize record time
#'
#' Shift all time variables by a random amount between 0 and ~32 years (1e9 secs)
#'
#' @param rec (nested) list of signals
#'
#' @export
anonymize_record_time <- function(rec) {
    offset_s <- runif(1, 0, 1e9)

    shift_record_time(rec, -offset_s)
}
