#' Extracting outputs from environment-wise mixed models
#'
#' @title Extracting outputs from environment-wise mixed models
#' 
#' @description
#' This function compiles the outputs of mixed models fitted using 
#' [asreml::asreml()], and computes several useful parameters, such as the 
#' coefficient of variation, Cullis' heritability, and likelihood ratio test (LRT) 
#' significance for each environment.
#' 
#' @param data A data frame containing the dataset. It must include at least 
#' the columns for genotype, block, environment, and the response variable.
#' @param response_var A string representing the response variable to be analyzed.
#' @param name.env A string representing the name of the "environment" factor in the dataset.
#' @param name.gen A string representing the name of the "genotype" factor in the dataset.
#'
#' @return The function returns a formatted table with:
#' \itemize{
#' \item \code{Environment}: Name of the environment.
#' \item \code{Statistic}: The computed parameter (Overall Mean, Coefficient of Variation, Cullis' Heritability, LRT p-value).
#' \item \code{Value}: The numerical value for each statistic.
#' \item \code{Significance}: Significance level based on LRT (\code{"***"}, \code{"**"}, \code{"*"}, \code{"."}, or blank).
#' }
#'
#' @details
#' This function fits an environment-wise mixed model with the structure:
#' \deqn{ Y = \mu + block + gen + \varepsilon }
#' where \code{gen} is treated as a random effect. The function also tests the 
#' significance of the genetic variance component using a likelihood ratio test (LRT).
#'
#' @examples
#' \dontrun{
#' # Load dataset
#' data <- read.table("data.txt", header = TRUE, na.strings = "NA")
#' data <- transform(data,
#'                   gen = factor(Genotype),
#'                   block = factor(Block),
#'                   env = factor(Env),
#'                   GY = as.numeric(GY))
#'
#' # Run function for GY
#' analyze_environments(data, response_var = "GY", name.env = "env", name.gen = "gen")
#' }
#' 
#' 
#' 
##' @author João Marcos Amario de Sousa (Amario at ufv.br)
#' 
analyze_environments <- function(data, response_var, name.env = "env", name.gen = "gen") {
  
  library(asreml)
  library(kableExtra)
  
  # Ensure response_var exists in data
  
  if (!(response_var %in% names(data))) {
    stop(paste("Error: The variable", response_var, "was not found in the dataset. Check the column names."))
  }
  
  # Ensure environment and genotype columns exist
  
  if (!(name.env %in% names(data)) | !(name.gen %in% names(data))) {
    stop("Error: The specified environment or genotype column was not found in the dataset.")
  }
  
  results_df <- data.frame(
    Environment = character(),
    Statistic = character(),
    Value = numeric(),
    Significance = character(),
    stringsAsFactors = FALSE)
  
  for (environment in unique(data[[name.env]])) {
    edat <- subset(data, data[[name.env]] == environment) |> droplevels()
    
# Fit the full model
    model_full <- asreml(
      fixed = as.formula(paste(response_var, "~ block")),
      random = as.formula(paste("~", name.gen)),
      workspace = "2gb",
      maxit = 100, 
      trace = FALSE, 
      na.action = na.method(x = "include", y = "include"), 
      data = edat)
    
    # Extract variance components and compute statistics
    
    vari1 <- summary(model_full)$varcomp
    vari_res <- vari1["units!R", "component"]
    overall_mean <- mean(edat[[response_var]], na.rm = TRUE)
    CV <- (sqrt(vari_res) / overall_mean) * 100
    
# Predictions and Cullis' Heritability
    
    predm1 <- predict(model_full, classify = name.gen, vcov = TRUE)
    H2 <- 1 - (predm1$avsed[2]^2) / (2 * vari1[1, 1])
    
# Likelihood Ratio Test (LRT) 
    
    model_reduced <- asreml(fixed = as.formula(paste(response_var, "~ block")),
                            workspace = "2gb",
                            maxit = 100, 
                            trace = FALSE, 
                            na.action = na.method(x = "include", y = "include"), 
                            data = edat)
    
    lrt_result <- lrt(model_full, model_reduced)
    
# Extract p-value
    
    p_val_str <- as.character(lrt_result[1, "Pr(Chisq)"])
    p_val <- as.numeric(sub(" .*", "", p_val_str))
    
    sig <- ifelse(p_val < 0.001, "***",
                  ifelse(p_val < 0.01, "**",
                         ifelse(p_val < 0.05, "*",
                                ifelse(p_val < 0.1, ".", " "))))
    
    results_env <- data.frame(
      Environment = environment,
      Statistic = c(paste("Overall Mean of", response_var),
                    "Coefficient of Variation (%)",
                    "Cullis' Heritability",
                    "LRT p-value"),
      Value = round(c(overall_mean, CV, H2, p_val), 4),
      Significance = c("", "", "", sig))
    
    results_df <- rbind(results_df, results_env)}
  
  envs_no_significance <- unique(results_df$Environment[results_df$Significance == " "])
  
  return(
    kbl(results_df, caption = paste("Model Results for", response_var, "by Environment"),
        align = "c", col.names = c("Environment", "Statistic", "Value", "Significance")) |>
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                    full_width = FALSE, position = "center") |>
      collapse_rows(columns = 1, valign = "middle") |>  
      row_spec(0, bold = TRUE, background = "#f5f5f5") |> 
      row_spec(seq(2, nrow(results_df), by = 2), background = "#f9f9f9") |> 
      row_spec(which(results_df$Environment %in% envs_no_significance), color = "red")
  )
}





