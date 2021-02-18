library(gratia)

round_to_interval <- function(x) {
    interval <- last(x) - first(x)
    
    signif_interval <- signif(interval * 0.05, 2)
    
    round(x/signif_interval) * signif_interval
}

# UI module
visualize_GAM_UI <- function(id) {
    tagList(
        
        sidebarPanel(
            #selectInput('x_scale', 'X scale'),
            
            uiOutput(NS(id,"varsInput"))
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput(NS(id, "gamPlot"))
        )
    )
    
}

# Module server
visualize_GAM_server <- function(id, gam_obj) {
    stopifnot(is.reactive(gam_obj))
    
    moduleServer(id, function(input, output, session) {
        
        output$varsInput <- renderUI({
            relevant_vars <- gam_obj()$var.summary[names(gam_obj()$var.summary) != "qrs_rel_index"]
            
            var_range <- lapply(relevant_vars, round_to_interval)
            
            var_names <- names(var_range)
            input_names <- paste0("slider_", var_names)
            
            output <- tagList()
            
            for(i in seq_along(var_range)){
                output[[i]] <-  sliderInput(session$ns(input_names[i]), var_names[i], 
                                            min = var_range[[i]][1],
                                            max = var_range[[i]][3],
                                            value = var_range[[i]][1])
            } 
            
            output
        })
        
        output$gamPlot <- renderPlot({
            req(length(reactiveValuesToList(input)) > 0)
            
            gam_abp_range <- range(gam_obj()$fitted.values)
            
            n_grid <- 101
            
            slider_data <- reactiveValuesToList(input)
            slider_data_df <- as_tibble(lapply(slider_data, rep, n_grid))
            names(slider_data_df) <- str_remove(names(slider_data_df), 'slider_')
            newdata <- tibble(
                qrs_rel_index = seq(0, 1, length.out = n_grid),
                
            ) %>% 
                bind_cols(slider_data_df)
            
            plot_data <- tibble(qrs_rel_index = newdata$qrs_rel_index, 
                                fit = predict(gam_obj(), newdata = newdata))
            
            ggplot(plot_data, aes(qrs_rel_index, fit)) +
                geom_line() +
                ylim(c(floor(gam_abp_range[1] / 10) * 10, 
                       ceiling(gam_abp_range[2] / 10) * 10))
            
        })
        
    })
}
