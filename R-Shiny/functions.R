gen_predictions <- function(cvp_gam) {
  # a function that generates predictions according to the given model
  # the code is a slightly modified version of Johannes'
  if (is.null(cvp_gam)) return(NULL)
  
  #find max of insp_index smooth
  #insp_rel_index_smooth <- gratia::evaluate_smooth(cvp_gam, 's(insp_rel_index)')
  #insp_rel_index_max <- insp_rel_index_smooth$insp_rel_index[which.max(insp_rel_index_smooth$est)]
  
  structured_data <- expand.grid(insp_rel_index = c(
    # Generate predictions for 10 levels of inspiration between start and max inspiration
    seq(0, 0.9, length.out = 10)), #insp_rel_index_max, length.out = 10)), 
    # And 1 prediction 1 sec after max inspiration                                      
    #insp_rel_index_max + c(1)), # uncommented for our use case
    # Generate 200 predictions uniformaly across the qrs_rel_index (to get a smooth line)
    qrs_rel_index = seq(0,1, length.out = 500),
    # Choose 'some' position for the trend smooth (or better, do not include it in the prediction)
    time_s = 580)
  
  struct_pred <- predict(cvp_gam, newdata = structured_data, se.fit = TRUE, exclude='s(time_s)') 
  
  mutate(structured_data, 
         CVP = struct_pred$fit,
         CVP_low = CVP - 1.96 * struct_pred$se.fit,
         CVP_up = CVP + 1.96 * struct_pred$se.fit)
}


create_patient_preds <- function(ptopen, ptclosed){
  # a function that creates predictions and combines data from closed and open chest into one data frame 
  pred_open <- gen_predictions(ptopen)
  pred_closed <- gen_predictions(ptclosed)
  pred_both <- cbind(indicator=c(rep("open_chest", nrow(pred_open)),
                                 rep("closed_chest", nrow(pred_closed))),
                     pred_both <- rbind(pred_open, pred_closed)) # makes an indicator variable wether it is open or closed chest
  return(pred_both)
}



stacked_plot <- function(io, ic, qo, qc, patient_number = NULL){
  # This function takes (inspiration/qrs) and (open/closed-chest) and outputs the plot for each combination.
  
  colnames(io)[3] <- "rel_index"
  colnames(ic)[3] <- "rel_index"
  colnames(qo)[3] <- "rel_index"
  colnames(qc)[3] <- "rel_index"
  
  comb <- cbind(indicator=c(rep("open_chest", nrow(io)),
                            rep("closed_chest", nrow(ic))),
                both <- rbind(io, ic, qo, qc)) # combines the data to make it easier to work with
  
  #creates the plots with ggplot
  viz <- ggplot(data=comb, aes(x=rel_index, y=est, group=indicator)) +
    geom_line(aes(color=indicator), size = 0.9) +
    facet_wrap(~ smooth) +
    labs(title = "The marginal effects of respiration and QRS, respectively, on CVP",
         y="Estimated effect on CVP",
         x="Relative respiration/QRS index",
         color = "Time period") +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
  return(viz)
}



stacked_diff_plot <- function(io, ic, qo, qc){
  # This function takes (inspiration/qrs) and (open/closed-chest) and plots the difference in open/closed-chest.
  
  colnames(io)[3] <- "rel_index"
  colnames(ic)[3] <- "rel_index"
  colnames(qo)[3] <- "rel_index"
  colnames(qc)[3] <- "rel_index"
  
  combo <- rbind(io, qo)
  combo <- combo[order(combo$smooth,combo$rel_index),] # sorts the data
  combc <- rbind(ic, qc)
  combc <- combc[order(combc$smooth,combc$rel_index),]
  
  combo$est_diff <- combo$est - combc$est # takes the difference
  
  # plots
  viz <- ggplot(data=combo, aes(x=rel_index, y=est_diff)) +
    geom_line(aes(), size = 0.9) +
    facet_wrap(~ smooth) +
    labs(title = "The difference in marginal effects (respiration and QRS) before and after thoractomy",
         y="Effect during open chest minus effect during closed chest",
         x="Relative respiration/QRS index") +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
  return(viz)
}


create_contour_patient <- function(odata, cdata){
  # this function takes open- and closed data. It is a helper function that combines the two datasets in one table 
  # with an added column indicating whether the chest is open or closed
  
  gather_in_one <- rbind(odata, cdata)
  gather_in_one$indicator <- cbind(c(rep("open_chest", nrow(odata)),
                                     rep("closed_chest", nrow(cdata))))
  return(gather_in_one)
}



move_contour <- function(data_contour, input_scale_open, input_scale_closed){
  # this function takes the combined dataset for open and closed and outputs a new dataset where the order has been rearranged
  data_contour <- mutate(data_contour,
                         est = c(est[((100 - input_scale_open) * 100 + 1):10000],
                                 est[1:((100 - input_scale_open) * 100)],
                                 est[((100 - input_scale_closed) * 100 + 10001):20000],
                                 est[10001:((100 - input_scale_closed) * 100 + 10000)]))
  return(data_contour)
}


create_contour_diff_data <- function(data_set){
  # this function takes the combined dataset and outputs the difference in open and closed chest.
  # this is done by subtracting the the open estimate from the closed estimate.
  data_set <- data_set[order(data_set$indicator, data_set$qrs_rel_index, data_set$insp_rel_index),]
  new_diff_data <- data.frame(qrs_rel_index = data_set$qrs_rel_index[1:10000],
                              insp_rel_index = data_set$insp_rel_index[1:10000],
                              est = data_set$est[10001:20000] - data_set$est[1:10000], # open minus closed est
                              indicator = rep("Difference", 10000))
  return(new_diff_data)
}


create_contour_viz <- function(data_set,
                               lab.title = "Effect of interaction term between respiration and QRS on CVP",
                               lab.y = "Relative respiration index",
                               lab.x = "Relative QRS index",
                               lab.legend = "Estimated effect\non CVP",
                               wrap = TRUE,
                               left_margin = 5){
  
  # This functions plots a contour plot for a dataset with x-axis as QRS and y as repiration index
  viz <- ggplot(data_set, aes(qrs_rel_index, insp_rel_index, z = est)) + # x,y and z dimensions
    geom_raster(aes(fill = est)) + # this creates the contour-plot
    geom_contour(colour = "black") +
    scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
    labs(title = lab.title,
         y = lab.y,
         x = lab.x,
         fill = lab.legend) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = left_margin, b = 0, l = 0)))
  if(wrap){
    viz <- viz + facet_wrap(~ indicator)}
  return(viz)
}


create_viz <- function(pred_both, patient_number = NULL){
  # a function that creates a visualization comparing closed and open predictions with 10 lines
  
  viz <- ggplot(data = pred_both, aes(x=qrs_rel_index, y=CVP, group = insp_rel_index)) + 
    # takes x and y values. Where it groups over 10 different insp_rel_index indekses 
    geom_line(aes(color=insp_rel_index), size = 0.9) +
    scale_color_gradientn(colours = turbo(10)) +
    labs(title = "Change in CVP as a function of relative QRS index at selected relative respiration index times",
         y = "CVP",
         x = "Relative QRS index",
         colour = "Relative        \nrespiration index") +
    facet_wrap(~ indicator) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
  return(viz)
}



diff_snap_plot <- function(preds_data, lab.legend = "Relative        \nrespiration index"){
  # this function takes the combined shifted data, sorts it, takes the difference in open and closed where it plots this difference for 10 different lines of indeks' 
  
  preds_data <- preds_data[order(preds_data$indicator, preds_data$insp_rel_index, preds_data$qrs_rel_index),]
  closed_and <- preds_data[1:5000,]
  open_and <- preds_data[5001:10000,]
  
  data_to_plot <- tibble(insp_rel_index = (closed_and$insp_rel_index + open_and$insp_rel_index)/2,
                         qrs_rel_index = closed_and$qrs_rel_index,
                         CVP = open_and$CVP - closed_and$CVP, # takes the difference in cvp
                         indicator = rep("Difference", 5000))
  
  # plots the difference
  viz <- ggplot(data = data_to_plot, aes(x=qrs_rel_index, y=CVP, group = insp_rel_index)) +
    geom_line(aes(color=insp_rel_index), size = 0.9) +
    scale_color_gradientn(colours = turbo(10)) +
    labs(title = "The difference in CVP as a function of relative QRS index at selected relative respiration index times",
         y = "Difference in CVP (open minus close)",
         x = "Relative QRS index",
         colour = lab.legend) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) + 
    facet_wrap(~ indicator)
  return(viz)
}

insp_shift_data <- function(data_set, shift_value_open, shift_value_closed){
  # This function also shifts data. This is done for the datasets containing 10 different lines
  data_set <- mutate(data_set, insp_rel_index = ifelse(indicator=="open_chest",
                                                       (data_set$insp_rel_index + shift_value_open) %% 1,
                                                       (data_set$insp_rel_index + shift_value_closed) %% 1))
  return(data_set)
}


give_it_center <- function(data_set, center = 50){
  # this function takes a dataset and calculates whether the data should be shifted 
  # to the left or right so the largest value is the middel index
  if (which.max(data_set$est) > center ){ # this is executed if the largest input is to the right of the center
    biggest <- which.max(data_set$est)
    value_to_subtract <- data_set$insp_rel_index[biggest] - data_set$insp_rel_index[center]
    data_set$insp_rel_index <- data_set$insp_rel_index - value_to_subtract
    for (index in 1:length(data_set$insp_rel_index)) {
      if (data_set$insp_rel_index[index] < 0 ){
        data_set$insp_rel_index[index] <- data_set$insp_rel_index[index] +1
      }
    }
  } else{ # this is executed if the largest input is to the left (or exactly on) the center
    smallest <- which.max(data_set$est)
    value_to_add <- data_set$insp_rel_index[center] - data_set$insp_rel_index[smallest]
    data_set$insp_rel_index <- data_set$insp_rel_index + value_to_add
    for (index in 1:length(data_set$insp_rel_index)) {
      if (data_set$insp_rel_index[index] > 1 ){
        data_set$insp_rel_index[index] <- data_set$insp_rel_index[index] -1
      }
    }
  }
  return(data_set)
}


insp_add_sub_val <- function(data_set, center = 50){
  # takes dataset and outputs value to add/subtract from index of a dataset depending on given center value. 
  if (which.max(data_set$est) > center ){
    biggest <- which.max(data_set$est)
    value_to_subtract <- data_set$insp_rel_index[center] - data_set$insp_rel_index[biggest]
    return(value_to_subtract)
  }else{
    smallest <- which.max(data_set$est)
    value_to_add <- data_set$insp_rel_index[center] - data_set$insp_rel_index[smallest]
    return(value_to_add)
  }
}




create_dygraphs <- function(obj, closed, ptnumber){
  # this is a helper function that styles the residual plots
  if (closed==TRUE){
    res1 <- dygraph_gam(obj, resid=TRUE)
    res1$x$group <- NULL
    res1$x$attrs$ylabel <- "CVP"
    res1$x$attrs$title <- paste0(ptnumber, " - Closed chest")
    res1$x$attrs$xlabel <- "Time"
  }
  else{
    res1 <- dygraph_gam(obj, resid=TRUE)
    res1$x$attrs$ylabel <- "CVP"
    res1$x$attrs$title <- paste0(ptnumber, " - Open chest")
    res1$x$attrs$xlabel <- "Time"
  }
  return(res1)
}