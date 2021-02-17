# Functions to work with gam models

#' Get PPV from GAM model
#'
#' @param model Gam model
#' @param term term to calculate ppv from (e.g. 's(resp_index)')
#' @param ppv_only Return only PPV value. Otherwise return row of related parameters as well.
#'
#' @export
get_gam_ppv <- function(model, term, ppv_only = FALSE) {
  if (!("gam" %in% class(model))) {
    if (length(model) == 1 && is.na(model)) {
      min_est <- NA
      max_est <- NA
      min_se <- NA
      max_se <- NA
      PPV <- NA
      PPV_se <- NA
      model_intercept <- NA
    }
    else {
      stop("model is not a GAM model")
    }
  } else {
    smooth <- gratia::evaluate_smooth(model, term)
    model_intercept <- model$coefficients[1]

    max_i <- which.max(smooth$est)
    min_i <- which.min(smooth$est)
    max_est <- smooth$est[max_i]
    max_se <- smooth$se[max_i]
    min_est <- smooth$est[min_i]
    min_se <- smooth$se[min_i]

    PPV <- (max_est - min_est) / model_intercept

    # Standard error of difference https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1125071/

    PPV_se <- sqrt(max_se^2 + min_se^2) / model_intercept
  }

  if (ppv_only) {
    return(PPV)
  }




  dplyr::tibble(min_est,
    max_est,
    min_se,
    max_se,
    PPV,
    PPV_se,
    PPV_lower = PPV - 1.96 * PPV_se,
    PPV_upper = PPV + 1.96 * PPV_se,
    model_intercept
  )
}

#' Create data frame of data, fit and residuals of pulse pressure GAM
#'
#' @param PP_gam GAM of PP
#'
#' @return data frame
#' @export
get_PP_gam_predictions <- function(PP_gam) {
  params <- purrr::map_chr(PP_gam$smooth, purrr::pluck, "label")

  dplyr::mutate(PP_gam$model,
    PP_predict = mgcv::predict.gam(PP_gam),
    PP_res = mgcv::residuals.gam(PP_gam),
    PP_insp = as.numeric(mgcv::predict.gam(PP_gam, type = "terms", terms = params[1])),
    PP_trend = as.numeric(mgcv::predict.gam(PP_gam, type = "terms", terms = params[2])),
    PP_detrend = PP_res + PP_insp,
    PP_mean = coef(PP_gam)[1]
  )
}

#' Plot pulse pressure GAM smooths and partial residuals
#'
#' @param PP_gam PP GAM
#' @param return_list If true, return list of plots instead of combined plot.
#'
#' @return Patchworked ggplot
#' @export
plot_PP_gam <- function(PP_gam, return_list = FALSE) {

  # Get smooth labels
  params <- purrr::map_chr(PP_gam$smooth, purrr::pluck, "label")
  terms <- purrr::map_chr(PP_gam$smooth, purrr::pluck, "term")



  beats_p <- get_PP_gam_predictions(PP_gam)

  # Model visualizations

  # Get representations of smooths
  insp_smooth <- gratia::evaluate_smooth(PP_gam,
    smooth = params[1],
    newdata = seq(0, 1, length.out = 100)
  )

  time_smooth <- gratia::evaluate_smooth(PP_gam,
                                         smooth = params[2])


  geom_ci_ribbon <- ggplot2::geom_ribbon(ggplot2::aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se),
    fill = ggplot2::alpha("black", 0),
    colour = "black", linetype = 2, outline.type = "both"
  )

  insp_smooth_plot <- ggplot2::ggplot(insp_smooth, ggplot2::aes_string(terms[1], "est")) +
    geom_ci_ribbon +
    ggplot2::geom_line() +
    ggplot2::geom_point(ggplot2::aes(y = PP_detrend), data = beats_p, alpha = 0.7) +
    ggplot2::labs(
      x = "Index to inspiration start (relative)",
      y = "Additive effect on mean PP [mmHg]"
    ) +
    ggplot2::ggtitle("Position in respiratory cycle") +
    ggplot2::scale_x_continuous(labels = scales::percent)

  time_smooth_plot <- ggplot2::ggplot(time_smooth, ggplot2::aes_string(terms[2], "est")) +
    geom_ci_ribbon +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Time [s]",
      y = "Additive effect on mean PP [mmHg]"
    ) +
    ggplot2::ggtitle("Trend over time")

  if (return_list) {
    list(insp_smooth_plot, time_smooth_plot)
  } else {
    if(requireNamespace('patchwork', quietly = TRUE)) {
      patchwork::wrap_plots(insp_smooth_plot, time_smooth_plot)
    } else {
      warning("returning list. Install {patchwork} to return combined plot")
      list(insp_smooth_plot, time_smooth_plot)
    }
  }

}

#' Plot model predictions and observations using Dygraphs
#'
#' @param ... One or more Gam models.
#' @param resid include residuals
#'
#' @export
dygraph_gam <- function(..., resid = FALSE) {
  gams <- rlang::list2(...)
  mod_names <- rlang::exprs(...)

  pred <- lapply(gams, predict)
  resids <- NULL

  if (resid) {
    resids <- lapply(gams, resid)
    names(resids) <- paste0(mod_names, "_resid")
  }

  names(pred) <- paste0(mod_names, "_pred")

  assertthat::assert_that(length(pred) == 1 || sd(sapply(pred, length)) == 0) # All lengths equal

  df <- dplyr::bind_cols(
    time = gams[[1]]$model$time,
    obs = gams[[1]]$model[,1],
    pred,
    resids
  )

  dygraph_signal(df, main = '')
}
