# Required libraries
packages <- c(
    "shiny",
    "shinythemes",
    "shinyjs",
    "readxl",
    "writexl",
    "DT",
    "caret",
    "randomForest"
)

# Function to install missing packages
install_if_missing <- function(package) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package)
        library(package, character.only = TRUE)
    }
}

# Install and load all required packages
sapply(packages, install_if_missing)

# Create www directory if it doesn't exist
if (!dir.exists("www")) {
    dir.create("www")
}

# Check for required model files
required_files <- c(
    "OvAge.RData",
    "UVA_FSH_Multiplex_US_normalization_model.RData"
)
