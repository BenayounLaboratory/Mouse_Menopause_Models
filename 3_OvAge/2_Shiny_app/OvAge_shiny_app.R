# Load libraries

library(shiny)
library(shinythemes)
library(shinyjs)
library(readxl) 
library(writexl)
library(DT)     
library(caret)  
library(randomForest)

################################################################################
# Menopause-model project
# Shiny app for OvAge prediction
################################################################################

# Set paths for model files
OVAGE_MODEL_PATH <- "./OvAge.RData"
FSH_MODEL_PATH <- "./UVA_FSH_Multiplex_US_normalization_model.RData"

# Load pre-trained models
load(OVAGE_MODEL_PATH)
rf_model <- OvAge

load(FSH_MODEL_PATH)
fsh_correction_model <- fsh_poly 

# Create example data frame for template
example_data <- data.frame(
    Sample_ID = c("Mouse1", "Mouse2", "Mouse3"),
    AMH = c(2.5, 1.8, 3.1),
    FSH = c(1.2, 0.9, 1.5),
    INHBA = c(150, 120, 180),
    FSH_Assay = c("millipore", "millipore", "uva"),
    stringsAsFactors = FALSE
)

# UI Definition
ui <- fluidPage(
    theme = shinytheme("flatly"),
    useShinyjs(),
    titlePanel("Mouse Age Prediction from Hormone Levels"),
    
    # Description Panel
    wellPanel(
        h4("Assay Details"),
        tags$ul(
            tags$li("AMH (ng/ml): Mouse Anti-Müllerian Hormone ELISA Kit (Ansh Labs, AL-113)"),
            tags$li(HTML("FSH (ng/ml): Two possible assays:<br/>",
                        "&nbsp;&nbsp;&nbsp;&nbsp;• Millipore Pituitary Panel Multiplex kit (RPT86K)<br/>",
                        "&nbsp;&nbsp;&nbsp;&nbsp;• Ultra-Sensitive Mouse & Rat FSH (UVA Ligand Core, in-house)")),
            tags$li("INHBA (pg/ml): Mouse Inhibin A ELISA Kit (Ansh Labs, AL-161)")
        )
    ),
    
    # Tab panels for single/batch modes
    tabsetPanel(
        # Single Sample Tab
        tabPanel("Single Sample",
            sidebarLayout(
                sidebarPanel(
                    # Input fields
                    numericInput("amh", 
                                "AMH Level (ng/ml):", 
                                value = NA,
                                min = 0),
                    
                    numericInput("fsh", 
                                "FSH Level (ng/ml):", 
                                value = NA,
                                min = 0),
                    
                    numericInput("inhba", 
                                "INHBA Level (pg/ml):", 
                                value = NA,
                                min = 0),
                    
                    # FSH Assay type selection
                    radioButtons("fsh_assay", 
                                "FSH Assay Type:",
                                choices = c(
                                    "Millipore ELISA Kit" = "millipore",
                                    "In-house UVA assay (requires correction)" = "uva"
                                )),
                    
                    # Submit button
                    actionButton("predict", "Predict Age", 
                                class = "btn-primary btn-lg btn-block")
                ),
                
                mainPanel(
                    # Output panel
                    wellPanel(
                        h4("Prediction Results"),
                        textOutput("prediction_text"),
                        textOutput("error_text")
                    )
                )
            )
        ),
        
        # Batch Processing Tab
        tabPanel("Batch Processing",
            fluidRow(
                column(12,
                    wellPanel(
                        h4("Batch Upload Instructions"),
                        tags$ul(
                            tags$li("Upload an Excel file (.xlsx) or CSV file with the following columns:"),
                            tags$li(HTML("Required columns and units:<br/>",
                                       "&nbsp;&nbsp;&nbsp;&nbsp;• Sample_ID: Unique identifier for each sample<br/>",
                                       "&nbsp;&nbsp;&nbsp;&nbsp;• AMH: Concentration in ng/ml<br/>",
                                       "&nbsp;&nbsp;&nbsp;&nbsp;• FSH: Concentration in ng/ml<br/>",
                                       "&nbsp;&nbsp;&nbsp;&nbsp;• INHBA: Concentration in pg/ml<br/>",
                                       "&nbsp;&nbsp;&nbsp;&nbsp;• FSH_Assay: Either 'millipore' or 'uva'")),
                            tags$li(downloadLink("download_template", "Download template file"))
                        )
                    )
                )
            ),
            fluidRow(
                column(6,
                    fileInput("file1", "Choose Excel/CSV File",
                        accept = c(".xlsx", ".csv")
                    ),
                    actionButton("process_batch", "Process Batch",
                        class = "btn-primary"
                    )
                ),
                column(6,
                    downloadButton("download_results", "Download Results",
                        class = "btn-success"
                    )
                )
            ),
            fluidRow(
                column(12,
                    # Processing status
                    textOutput("batch_status"),
                    # Error messages
                    textOutput("batch_error"),
                    # Results table
                    DTOutput("batch_results")
                )
            )
        )
    )
)

# Server
server <- function(input, output, session) {
    
    # Template download handler
    output$download_template <- downloadHandler(
        filename = function() {
            "template.xlsx"
        },
        content = function(file) {
            write_xlsx(example_data, file)
        }
    )
    
    # Validation function
    validate_inputs <- function(amh, fsh, inhba) {
        if (is.na(amh) || is.na(fsh) || is.na(inhba)) {
            return("Please fill in all hormone levels")
        }
        if (amh < 0 || fsh < 0 || inhba < 0) {
            return("All hormone levels must be non-negative")
        }
        return(NULL)
    }
    
    # Single prediction function
    predict_age <- function(amh, fsh, inhba, fsh_assay) {
        # Apply FSH correction if UVA assay is selected
        if (fsh_assay == "uva") {
            # Create prediction data frame with correct column name
            correction_data <- data.frame(FSH_US = as.numeric(fsh))
            fsh <- predict(fsh_correction_model, newdata = correction_data)
        }
        
        # Create prediction data frame
        pred_data <- data.frame(
            AMH = as.numeric(amh),
            FSH = as.numeric(fsh),
            INHBA = as.numeric(inhba)
        )
        
        # Predict age using caret's predict method
        predicted_age <- predict(rf_model, newdata = pred_data)
        return(round(as.numeric(predicted_age), 1))
    }
    
    # Single sample prediction
    observeEvent(input$predict, {
        output$error_text <- renderText("")
        
        error_msg <- validate_inputs(input$amh, input$fsh, input$inhba)
        if (!is.null(error_msg)) {
            output$error_text <- renderText({
                paste("Error:", error_msg)
            })
            return()
        }
        
        tryCatch({
            predicted_age <- predict_age(input$amh, input$fsh, input$inhba, input$fsh_assay)
            output$prediction_text <- renderText({
                paste("Predicted mouse age:", predicted_age, "weeks")
            })
        }, error = function(e) {
            output$error_text <- renderText({
                paste("Error in prediction:", e$message)
            })
        })
    })
    
    # Batch processing
    batch_results <- reactiveVal(NULL)
    
    # Batch status
    output$batch_status <- renderText({
        if (is.null(input$file1)) return("Upload a file to begin")
        if (is.null(batch_results())) return("Click 'Process Batch' to analyze data")
        return("Processing complete!")
    })
    
    # Batch error messages
    output$batch_error <- renderText({
        if (is.null(input$file1)) return("")
        if (!is.null(batch_results()) && any(!is.na(batch_results()$Error))) {
            return("Warning: Some samples had errors. Check the 'Error' column for details.")
        }
        return("")
    })
    
    observeEvent(input$process_batch, {
        req(input$file1)
        
        # Show processing status
        output$batch_status <- renderText("Processing file...")
        
        # Read the uploaded file
        df <- tryCatch({
            if (tools::file_ext(input$file1$datapath) == "csv") {
                read.csv(input$file1$datapath, stringsAsFactors = FALSE)
            } else {
                read_excel(input$file1$datapath)
            }
        }, error = function(e) {
            output$batch_error <- renderText(paste("Error reading file:", e$message))
            return(NULL)
        })
        
        if (is.null(df)) return()
        
        # Validate columns
        required_cols <- c("Sample_ID", "AMH", "FSH", "INHBA", "FSH_Assay")
        missing_cols <- setdiff(required_cols, names(df))
        
        if (length(missing_cols) > 0) {
            output$batch_error <- renderText(
                paste("Missing required columns:", paste(missing_cols, collapse = ", "))
            )
            return()
        }
        
        # Process each row
        df$Predicted_Age <- NA
        df$Error <- NA
        
        for (i in 1:nrow(df)) {
            tryCatch({
                df$Predicted_Age[i] <- predict_age(
                    as.numeric(df$AMH[i]),
                    as.numeric(df$FSH[i]),
                    as.numeric(df$INHBA[i]),
                    as.character(df$FSH_Assay[i])
                )
            }, error = function(e) {
                df$Error[i] <- as.character(e$message)
            })
        }
        
        batch_results(df)
        
        # Update status
        output$batch_status <- renderText("Processing complete!")
    })
    
    # Display batch results
    output$batch_results <- renderDT({
        req(batch_results())
        datatable(batch_results(),
                 options = list(pageLength = 10,
                              scrollX = TRUE))
    })
    
    # Download handler
    output$download_results <- downloadHandler(
        filename = function() {
            paste0("age_predictions_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
        },
        content = function(file) {
            req(batch_results())
            write_xlsx(batch_results(), file)
        }
    )
}

# Run the app
shinyApp(ui = ui, server = server) 