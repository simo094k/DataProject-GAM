#
# This is a Shiny web application. 
# It is meant as a visualisation tool for researchers to use, where they in an easy manner 
# can inspect their GAM-objects and then compare effect of the heart (QRS-complex) and 
# lungs (respiration) as well as their interaction on the CVP before and after thorax has been opened.
#
# You can run the application by clicking the 'Run App' button above.
#
#

### Dependencies ###

rm(list = ls()) # clearing workspace
library(viridis)
library(waveformtools)
library(tidyverse)
library(mgcv)
library(hms)
library(shiny)
library(shinythemes)
library(shinyTime)
library(data.table)
library(readr)
library(dygraphs)
library(dplyr)
library(gridExtra)

### Load the functions from the functions.R file ###
source("functions.R")

### Load data ###
table_values <- read_csv("../GAM_RDS/CV_MAE.csv") 
bic_values <- read_csv("../GAM_RDS/CV_BIC.csv") 
mae_values <- read_csv("../GAM_RDS/CV_MAE_BIC.csv") 

### Load patients ###
patient_list <- paste('pt', str_pad(c(2, 3, 4, 5, 7, 9, 10, 11, 12, 14, 15, 16, 17, 18, 21, 22, 23, 24, 28, 31, 32, 
                                      35, 37, 38, 41, 42, 44, 45, 46, 47, 48, 51, 52, 53, 54, 55, 57, 58), 2, pad = "0"), sep = "")

# Define UI for application 
ui <- fluidPage(
    tags$head(
        tags$script(
            HTML("
            $(document).ready(function(){
              // Mark columns we want to toggle
              $('body').find('div [class=col-sm-4]').addClass('sidebarPanel');
              $('body').find('div [class=col-sm-8]').addClass('mainPanel');
            })


            Shiny.addCustomMessageHandler ('resize',function (message) {
              $('.sidebarPanel').toggle();
              $('.mainPanel').toggleClass('col-sm-8 col-sm-12');
              $(window).trigger('resize')
            });

           "))
    ),
    
    actionButton("showpanel", "Show/hide sidebar"),
    theme = shinytheme("sandstone"),
    navbarPage("GAM on heart-lung interaction during thoractomy", # Page header
               
               # First page begins
               tabPanel("Main page",
                        sidebarLayout(
                            
                            # Parameters
                            sidebarPanel(selectInput("patient", label = "Select patient",
                                                     choices = list(patient_list = patient_list),  selected = "pt02"),
                                         sliderInput("center_closed", label = "Index placement of largest respiration effect during closed chest",
                                                     min = 1, max = 99, value = 50),
                                         sliderInput("center_open", label = "Index placement of largest respiration effect during open chest",
                                                     min = 1, max = 99, value = 50),
                                         tags$div(class="multicol", checkboxGroupInput("showhide", "Plots to show",
                                                                                       choices = c("Regular plots" = "reg", "Closed/open difference plots" = "diff"),
                                                                                       selected = c("reg"))),
                                         h3("Model details"),
                                         tableOutput('table'),
                                         tags$style(type="text/css", "#table tr:first-child, td:first-child { font-weight: bold;}"),
                                         tags$style(type="text/css",
                                                    "#table tr:nth-child(1), tr:nth-child(5), tr:nth-child(10) {border-bottom: solid 2px black;}"),
                                         tags$style(type="text/css", "#table tr:hover {background-color: #f8fff5;"),
                                         htmlOutput("modeltext"),
                                         uiOutput("ex0"),
                                         uiOutput("ex1"),
                                         uiOutput("ex2")),
                           
                            
                            # Main panel for displaying outputs
                            mainPanel(plotOutput(outputId = "all_plots", height = "1200px", width = "100%")))
               ), 
               
               # Second page begins
               tabPanel("Residuals",
                        sidebarPanel(h3("Description"), htmlOutput("text2")),
                        mainPanel(dygraphOutput("dygraph1"), dygraphOutput("dygraph2"))
               ),
               
               # Third page begins
               tabPanel("Info",  # Skriv her
                        "Created by Andreas, Casper, Mads & Simon")
               
               
               
    )
) # UI ends

# Define server logic 
server <- function(input, output, session) {
    
    # Show/hide sidepanel
    observeEvent(input$showpanel,{
        session$sendCustomMessage(type = 'resize', message = 1)
    })
    
    
    # Loading data
    patient_path_open <- reactive({paste("../GAM_RDS/", input$patient, "OPEN.rds",sep="") })
    pt_data_open <- reactive({readRDS(patient_path_open(), refhook = NULL) })
    
    patient_path_closed <- reactive({paste("../GAM_RDS/", input$patient, "CLOSED.rds",sep="") })
    pt_data_closed <- reactive({readRDS(patient_path_closed(), refhook = NULL) })
    
    open_add_sub_value <- reactive({insp_add_sub_val(gratia::evaluate_smooth(pt_data_open(),'s(insp_rel_index)'),
                                                     input$center_open)})
    closed_add_sub_value <- reactive({insp_add_sub_val(gratia::evaluate_smooth(pt_data_closed(),'s(insp_rel_index)'), 
                                                       input$center_closed)})
    
    # For snap plots
    preds <- reactive({insp_shift_data(create_patient_preds(pt_data_open(), pt_data_closed()),
                                       open_add_sub_value(), 
                                       closed_add_sub_value())})
    snapplot <- reactive({create_viz(pred_both = preds())})
    
    # For snap diff plots
    snapplot_diff <- reactive({diff_snap_plot(preds(), lab.legend = "Relative \nrespiration index                   ")})
    
    # For stacked plots
    stackedp <- reactive({stacked_plot(
        give_it_center(gratia::evaluate_smooth(pt_data_open(),'s(insp_rel_index)'),input$center_open),
        give_it_center(gratia::evaluate_smooth(pt_data_closed(),'s(insp_rel_index)'),input$center_closed),
        gratia::evaluate_smooth(pt_data_open(),'s(qrs_rel_index)'),
        gratia::evaluate_smooth(pt_data_closed(),'s(qrs_rel_index)'))})
    
    # For stacked diff plots
    stackeddiffp <- reactive({stacked_diff_plot(
        give_it_center(gratia::evaluate_smooth(pt_data_open(),'s(insp_rel_index)'),input$center_open),
        give_it_center(gratia::evaluate_smooth(pt_data_closed(),'s(insp_rel_index)'),input$center_closed),
        gratia::evaluate_smooth(pt_data_open(),'s(qrs_rel_index)'),
        gratia::evaluate_smooth(pt_data_closed(),'s(qrs_rel_index)'))})
    
    
    # For contour plots
    data_open_contour <- reactive({gratia::evaluate_smooth(pt_data_open(), 'ti(qrs_rel_index,insp_rel_index)')})
    data_closed_contour <- reactive({gratia::evaluate_smooth(pt_data_closed(), 'ti(qrs_rel_index,insp_rel_index)')})
    
    contourp <- reactive({create_contour_viz(move_contour(create_contour_patient(data_open_contour(), 
                                                                                 data_closed_contour()),
                                                          input$center_open,
                                                          input$center_closed))})
    # For contour difference plots
    contourdiffp <- reactive({create_contour_viz(create_contour_diff_data(move_contour(create_contour_patient(
        data_open_contour(),
        data_closed_contour()),
        input$center_open,
        input$center_closed)),
        lab.legend = "Interaction effect during   \nopen chest minus interaction \neffect during closed chest",
        lab.title = "The difference in interaction effect (respiration and QRS) before and after thoractomy",
        wrap = TRUE,
        left_margin = 12)})
    
    
    
    # Rendering plots
    output$all_plots <- renderPlot({
        
        if (length(input$showhide) == 0) {
            ggplot(data.frame())
        } else {
            if (length(input$showhide) == 2){
                plotlist <- list(stackedp(), contourp(), snapplot(), stackeddiffp(), contourdiffp(), snapplot_diff())
                plotmatrix <- rbind(c(1,1,1,1,1,4,4,4,4,4),
                                    c(2,2,2,2,2,5,5,5,5,5),
                                    c(3,3,3,3,3,6,6,6,6,6))
            } else {
                if ("reg" == input$showhide){
                    plotlist <- list(stackedp(), contourp(), snapplot())
                    plotmatrix <- rbind(c(1),
                                        c(2),
                                        c(3))
                }
                if ("diff" == input$showhide){
                    plotlist <- list(stackeddiffp(), contourdiffp(), snapplot_diff())
                    plotmatrix <- rbind(c(4),
                                        c(5),
                                        c(6))
                }
            }
            grid.arrange(grobs = plotlist, layout_matrix = plotmatrix, nrow = 3)
        }
    }) # renderPlot ends here
    
    # Dygraphs
    dy_graph1 <- reactive({create_dygraphs(pt_data_closed(), closed=TRUE, input$patient)})
    dy_graph2 <- reactive({create_dygraphs(pt_data_open(), closed=FALSE, input$patient)})
    
    output$dygraph1 <- renderDygraph({dy_graph1()})
    output$dygraph2 <- renderDygraph({dy_graph2()})
    
    # Dygraph render ends here
    
    output$text2 <- renderText({
        HTML('On this page you can see the observed values, the predicted values of the models, and the residuals.
         <br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>
         <br><br><br><br><br><br><br>')
    }) # Render Text for second page
    
    # This is a TABLE!
    data_time_open  <- reactive({gratia::evaluate_smooth(pt_data_open(), 's(time_s)')})
    data_time_closed <- reactive({gratia::evaluate_smooth(pt_data_closed(), 's(time_s)')})
    
    
    model_table <- reactive({rename(subset(table_values,
                                           Model %in% c(paste0(input$patient, "OPEN.rds"),
                                                        paste0(input$patient, "CLOSED.rds")),
                                           select = c("Model")),
                                    `Model file` = Model)})
    
    cv_mae_table <- reactive({round(rename(subset(table_values,
                                                  Model %in% c(paste0(input$patient, "OPEN.rds"),
                                                               paste0(input$patient, "CLOSED.rds")),
                                                  select = c("CV_MAE")),
                                           `Cross validation MAE` = CV_MAE),
                                    3)})
    
    mae_table <- reactive({round(rename(subset(mae_values,
                                               Model %in% c(paste0(input$patient, "OPEN.rds"),
                                                            paste0(input$patient, "CLOSED.rds")),
                                               select = c("Actual Model MAE")),
                                        `Model MAE` = `Actual Model MAE`),
                                 3)})
    
    mae_table2 <- reactive({round(subset(mae_values,
                                         Model %in% c(paste0(input$patient, "OPEN.rds"),
                                                      paste0(input$patient, "CLOSED.rds")),
                                         select = c("LM MAE","Poly1 MAE","Poly2 MAE")),
                                  3)})
    
    bic_table <- reactive({round(rename(subset(bic_values,
                                               Model %in% c(paste0(input$patient, "OPEN.rds"),
                                                            paste0(input$patient, "CLOSED.rds")),
                                               select = c("Actual Model BIC","LM BIC","Poly1 BIC","Poly2 BIC")),
                                        `Model BIC`=`Actual Model BIC`))})
    
    
    intercept <- reactive({round(c(pt_data_open()[["coefficients"]][1],
                                   pt_data_closed()[["coefficients"]][1]),3)})
    
    start_time <- reactive({formatC(round(c(min(data_time_closed()$time_s), min(data_time_open()$time_s)), 1), 1, format="f")})
    end_time <- reactive({formatC(round(c(max(data_time_closed()$time_s), max(data_time_open()$time_s)), 1), 1, format="f")})
    diff_time <- reactive({formatC(round(c(max(data_time_closed()$time_s) - min(data_time_closed()$time_s),
                                           max(data_time_open()$time_s) - min(data_time_open()$time_s)), 1), 1, format="f")})
    
    output$table <- renderTable({t(cbind(model_table(),
                                         `Model intercept` = intercept(),
                                         `Model data start time` = start_time(),
                                         `Model data end time` = end_time(),
                                         `Model data time length` = diff_time(),
                                         mae_table(),
                                         cv_mae_table(),
                                         mae_table2(),
                                         bic_table()))}, 
                                rownames=TRUE, 
                                colnames=FALSE)
    
    output$modeltext <- renderUI({
        str1 <- "LM: Linear Regression Model"
        str11 <- "CVP = &beta;<sub>0</sub> (qrs_rel_index &times; insp_rel_index) + &beta;<sub>1</sub> qrs_rel_index + 
                  &beta;<sub>2</sub> insp_rel_index + &beta;<sub>3</sub> time + &epsilon;."
        str111 <- ""
        str2 <- "Poly1: Polynomial Regression Model 1"
        str22 <- "CVP = &beta;<sub>0</sub> qrs_rel_index<sup>12</sup> + &beta;<sub>1</sub> insp_rel_index<sup>8</sup> + 
                  &beta;<sub>2</sub> time<sup>3</sup> + &epsilon;."
        str222 <- ""
        str3 <- "Poly 2: Polynomial Regression Model 2"
        str33 <- "CVP = &beta;<sub>0</sub> qrs_rel_index<sup>23</sup> + &beta;<sub>1</sub> insp_rel_index<sup>23</sup> + 
                  &beta;<sub>2</sub> time<sup>23</sup> + &epsilon;."
        str333 <- ""
        HTML(paste(str1, str11, str111, str2, str22, str222, str3, str33, str333, sep = '<br/>') )
    })   
    
    # Modeller
    output$ex0 <- renderUI({
      withMathJax('LM: Linear Regression Model $$ \\scriptsize{  \\widehat{CVP} =\\beta_{1} \\cdot \\text{qrs_rel_index} +  \\gamma_{1} \\cdot
                          \\text{insp_rel_index} + \\zeta_{1} \\cdot \\text{time} +
                           \\epsilon }   $$')
    })
    output$ex1 <- renderUI({
      withMathJax('Poly 1: Polynomial Regression Model 1 $$ \\scriptsize{  \\widehat{CVP} =\\sum_{i=1}^{12} \\beta_{i} \\cdot \\text{qrs_rel_index}^{i} +  \\sum_{j=1}^{8} \\gamma_{j} \\cdot
                          \\text{insp_rel_index}^{j} + \\sum_{k=1}^{3} \\zeta_{k} \\cdot \\text{time}^{k} +
                           \\epsilon }   $$')
    })
    output$ex2 <- renderUI({
      withMathJax('Poly 2: Polynomial Regression Model 2 $$ \\scriptsize{  \\widehat{CVP} =\\sum_{i=1}^{23} \\beta_{i} \\cdot \\text{qrs_rel_index}^{i} +  \\sum_{j=1}^{23} \\gamma_{j} \\cdot
                          \\text{insp_rel_index}^{j}  + \\sum_{k=1}^{23} \\zeta_{k} \\cdot \\text{time}^{k} +
                           \\epsilon }   $$')
    })
    
    
} # server ends

# Run the application 
shinyApp(ui = ui, server = server)
