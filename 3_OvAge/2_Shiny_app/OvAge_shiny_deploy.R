# Load libraries
library(rsconnect)

################################################################################
# Menopause-model project
# OvAge shiny app
################################################################################

# Set account info
rsconnect::setAccountInfo(
    name='minhooki',
    token='',
    secret=''
)

# Terminate all instances of the app
tryCatch({
    rsconnect::terminateApp("OvAge_Predictor", account = "minhooki")
}, error = function(e) {
    cat("Note: No existing application found to terminate\n")
})

# Create a temporary directory for deployment
temp_dir <- tempfile("shiny_deploy_")
dir.create(temp_dir)

# Copy only the necessary files
file.copy("app.R", temp_dir)
file.copy("OvAge.RData", temp_dir)
file.copy("UVA_FSH_Multiplex_US_normalization_model.RData", temp_dir)
dir.create(file.path(temp_dir, "www"))
file.copy("www/template.xlsx", file.path(temp_dir, "www"))

# Deploy the application from the temporary directory
tryCatch({
    rsconnect::deployApp(
        appDir = temp_dir,
        appName = "OvAge_Predictor",
        appTitle = "Mouse OvAge Predictor",
        forceUpdate = TRUE,
        launch.browser = FALSE
    )
}, error = function(e) {
    cat("Deployment error:", e$message, "\n")
})

# Clean up
unlink(temp_dir, recursive = TRUE) 