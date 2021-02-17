#' Generate volume waveform from flow
#' @param vent_df ventilation waveform data frame
#'
#' @export
add_volume_signal <- function(vent_df) {
    samp_interval <- 1000 / attributes(vent_df$flow)$signal.samplerate

    assertthat::assert_that(assertthat::are_equal(samp_interval, 20, tol = 0.01)) # 50 Hz

    # Calculate change in volume between two flow measurements
    # L/min -> ml/ms = 1000/60e3
    res <- dplyr::mutate(vent_df, .dVol = c(flow[1], (flow[1:length(flow)-1] + flow[2:length(flow)])/2) *
                             (1000/60e3) * 20) # Trapezoid integration rule
    res <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(res, insp), volume = cumsum(.dVol)))
    dplyr::filter(res, insp != 0)
}
