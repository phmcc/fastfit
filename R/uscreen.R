#' Univariable Screening for Multiple Predictors
#'
#' Performs comprehensive univariable (unadjusted) regression analyses by fitting 
#' separate models for each predictor against a single outcome. This function is 
#' designed for initial variable screening, hypothesis generation, and understanding 
#' crude associations before multivariable modeling. Returns publication-ready 
#' formatted results with optional p-value filtering.
#'
#' @param data A data.frame or data.table containing the analysis dataset. The 
#'   function automatically converts data.frames to data.tables for efficient 
#'   processing.
#'   
#' @param outcome Character string specifying the outcome variable name. For 
#'   survival analysis, use \code{Surv()} syntax from the survival package 
#'   (e.g., \code{"Surv(time, status)"} or \code{"Surv(os_months, os_status)"}).
#'   
#' @param predictors Character vector of predictor variable names to screen. Each 
#'   predictor is tested independently in its own univariable model. Can include 
#'   continuous, categorical (factor), or binary variables.
#'   
#' @param model_type Character string specifying the type of regression model to 
#'   fit. Options include:
#'   \itemize{
#'     \item \code{"glm"} - Generalized linear model [default]
#'     \item \code{"lm"} - Linear regression
#'     \item \code{"coxph"} - Cox proportional hazards (survival analysis)
#'     \item \code{"clogit"} - Conditional logistic regression
#'   }
#'   
#' @param family For GLM models, the error distribution and link function. Common 
#'   options include:
#'   \itemize{
#'     \item \code{"binomial"} - Logistic regression for binary outcomes [default]
#'     \item \code{"poisson"} - Poisson regression for count data
#'     \item \code{"gaussian"} - Normal linear regression via GLM
#'     \item \code{"Gamma"} - Gamma regression for positive continuous data
#'   }
#'   See \code{\link[stats]{family}} for all options. Ignored for non-GLM models.
#'   
#' @param p_threshold Numeric value between 0 and 1 specifying a p-value threshold 
#'   for filtering results. Only predictors with p-value ≤ threshold in their 
#'   univariable model are included in the output. Default is 1 (no filtering, 
#'   all predictors returned). Common values: 0.05, 0.10, 0.20 for screening.
#'   
#' @param conf_level Numeric confidence level for confidence intervals. Must be 
#'   between 0 and 1. Default is 0.95 (95\% confidence intervals).
#'   
#' @param add_reference_rows Logical. If \code{TRUE}, adds rows for reference 
#'   categories of factor variables with baseline values (OR/HR/RR = 1, 
#'   coefficient = 0). Makes tables complete and easier to interpret. 
#'   Default is \code{TRUE}.
#'   
#' @param show_n Logical. If \code{TRUE}, includes the sample size column in 
#'   the output table. Default is \code{TRUE}.
#'   
#' @param show_events Logical. If \code{TRUE}, includes the events column in the 
#'   output table (relevant for survival and logistic regression). Default is 
#'   \code{TRUE}.
#'   
#' @param digits Integer specifying the number of decimal places for effect 
#'   estimates (OR, HR, RR, coefficients). Default is 2.
#'   
#' @param digits_p Integer specifying the number of decimal places for p-values. 
#'   P-values smaller than \code{10^(-digits_p)} are displayed as "< 0.001", 
#'   "< 0.01", etc. Default is 3.
#'   
#' @param labels Named character vector or list providing custom display 
#'   labels for variables. Names should match predictor names, values are the 
#'   display labels. Predictors not in \code{labels} use their original names. 
#'   Default is \code{NULL}.
#'   
#' @param keep_models Logical. If \code{TRUE}, stores all fitted model objects 
#'   in the output as an attribute. This allows access to models for diagnostics, 
#'   predictions, or further analysis, but can consume significant memory for 
#'   large datasets or many predictors. Models are accessible via 
#'   \code{attr(result, "models")}. Default is \code{FALSE}.
#'   
#' @param exponentiate Logical. Whether to exponentiate coefficients (display 
#'   OR/HR/RR instead of log odds/log hazards). Default is \code{NULL}, which 
#'   automatically exponentiates for logistic, Poisson, and Cox models, and 
#'   displays raw coefficients for linear models and other GLM families. Set 
#'   to \code{TRUE} to force exponentiation or \code{FALSE} to force coefficients.
#'   
#' @param ... Additional arguments passed to the underlying model fitting functions 
#'   (\code{\link[stats]{glm}}, \code{\link[stats]{lm}}, 
#'   \code{\link[survival]{coxph}}, etc.). Common options include \code{weights}, 
#'   \code{subset}, \code{na.action}, and model-specific control parameters.
#'
#' @return A data.table with S3 class \code{"uscreen_result"} containing formatted 
#'   univariable screening results. The table structure includes:
#'   \describe{
#'     \item{Variable}{Character. Predictor name or custom label (from \code{labels})}
#'     \item{Group}{Character. For factor variables: category level. For continuous 
#'       variables: typically empty or descriptive statistic label}
#'     \item{n}{Integer. Sample size used in the model (if \code{show_n = TRUE})}
#'     \item{n_group}{Integer. Sample size for this specific factor level 
#'       (factor variables only)}
#'     \item{events}{Integer. Total number of events in the model for survival 
#'       or logistic regression (if \code{show_events = TRUE})}
#'     \item{events_group}{Integer. Number of events for this specific factor 
#'       level (factor variables only)}
#'     \item{Univariable OR/HR/RR/Estimate (95\% CI)}{Character. Formatted effect 
#'       estimate with confidence interval. Column name depends on model type:
#'       "Univariable OR (95\% CI)" for logistic, "Univariable HR (95\% CI)" for 
#'       Cox, "Univariable RR (95\% CI)" for Poisson, "Univariable Estimate (95\% CI)" 
#'       for linear models}
#'     \item{p-value}{Character. Formatted p-value from the Wald test}
#'   }
#'   
#'   The returned object includes the following attributes accessible via \code{attr()}:
#'   \describe{
#'     \item{raw_data}{data.table. Unformatted numeric results with separate 
#'       columns for coefficients, standard errors, confidence interval bounds, 
#'       etc. Suitable for further statistical analysis or custom formatting}
#'     \item{models}{list (if \code{keep_models = TRUE}). Named list of fitted 
#'       model objects, with predictor names as list names. Access specific models 
#'       via \code{attr(result, "models")[["predictor_name"]]}}
#'     \item{outcome}{Character. The outcome variable name used}
#'     \item{model_type}{Character. The regression model type used}
#'     \item{model_scope}{Character. Always "Univariable" for screening results}
#'     \item{screening_type}{Character. Always "univariable" to identify the 
#'       analysis type}
#'   }
#'
#' @details
#' \strong{Analysis Approach:}
#' 
#' The function implements a comprehensive univariable screening workflow:
#' \enumerate{
#'   \item For each predictor in \code{predictors}, fits a separate model: 
#'     \code{outcome ~ predictor}
#'   \item Extracts coefficients, confidence intervals, and p-values from each model
#'   \item Combines results into a single table for easy comparison
#'   \item Optionally filters results based on \code{p_threshold}
#'   \item Formats output for publication with appropriate effect measures
#' }
#' 
#' Each predictor is tested \emph{independently} - these are crude (unadjusted) 
#' associations that do not account for confounding or interaction effects.
#' 
#' \strong{When to Use Univariable Screening:}
#' \itemize{
#'   \item \strong{Initial variable selection}: Identify predictors associated 
#'     with the outcome before building multivariable models
#'   \item \strong{Hypothesis generation}: Explore potential associations in 
#'     exploratory analyses
#'   \item \strong{Understanding crude associations}: Report unadjusted effects 
#'     alongside adjusted estimates
#'   \item \strong{Variable reduction}: Use p-value thresholds (e.g., p < 0.20) 
#'     to reduce the number of candidates for multivariable modeling
#'   \item \strong{Checking multicollinearity}: Compare univariable and 
#'     multivariable effects to identify potential collinearity
#' }
#' 
#' \strong{Factor Variables and Reference Categories:}
#' 
#' When \code{add_reference_rows = TRUE} (default):
#' \itemize{
#'   \item Reference categories are explicitly shown with OR/HR/RR = 1.00
#'   \item The reference row displays "(Reference)" instead of an effect estimate
#'   \item P-values are shown only for non-reference categories
#'   \item Group-specific sample sizes and event counts are calculated
#' }
#' 
#' \strong{P-value Filtering:}
#' 
#' The \code{p_threshold} parameter enables automatic filtering:
#' \itemize{
#'   \item Only predictors with at least one significant term (p ≤ threshold) 
#'     are retained
#'   \item For factor variables, if any level is significant, all levels are kept
#'   \item Common thresholds: 0.05 (strict), 0.10 (moderate), 0.20 (liberal screening)
#'   \item When \code{keep_models = TRUE}, non-significant models are also removed 
#'     from the models list
#' }
#' 
#' \strong{Effect Measures by Model Type:}
#' \itemize{
#'   \item \strong{Logistic regression} (\code{model_type = "glm"}, 
#'     \code{family = "binomial"}): Odds ratios (OR)
#'   \item \strong{Cox regression} (\code{model_type = "coxph"}): Hazard ratios (HR)
#'   \item \strong{Poisson regression} (\code{model_type = "glm"}, 
#'     \code{family = "poisson"}): Rate/risk ratios (RR)
#'   \item \strong{Linear regression} (\code{model_type = "lm"} or GLM with 
#'     identity link): Raw coefficient estimates
#' }
#' 
#' \strong{Memory Considerations:}
#' 
#' When \code{keep_models = FALSE} (default), fitted models are discarded after 
#' extracting results to conserve memory. Set \code{keep_models = TRUE} only when 
#' you need:
#' \itemize{
#'   \item Model diagnostic plots
#'   \item Predictions from individual models
#'   \item Additional model statistics not extracted by default
#'   \item Further analysis of specific models
#' }
#'
#' @seealso 
#' \code{\link{fit}} for fitting a single multivariable model,
#' \code{\link{fastfit}} for complete univariable-to-multivariable workflow,
#' \code{\link{compfit}} for comparing multiple models,
#' \code{\link{m2dt}} for converting individual models to tables
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' # Example 1: Basic logistic regression screening
#' screen1 <- uscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "hypertension"),
#'     model_type = "glm",
#'     family = "binomial"
#' )
#' print(screen1)
#' 
#' # Example 2: With custom variable labels
#' screen2 <- uscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "treatment"),
#'     labels = clintrial_labels
#' )
#' print(screen2)
#' 
#' # Example 3: Filter by p-value threshold
#' # Only keep predictors with p < 0.20 (common for screening)
#' screen3 <- uscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "hypertension", 
#'                   "diabetes", "ecog", "stage"),
#'     p_threshold = 0.20,
#'     labels = clintrial_labels
#' )
#' print(screen3)
#' # Only significant predictors are shown
#' 
#' # Example 4: Cox proportional hazards screening
#' library(survival)
#' cox_screen <- uscreen(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "treatment", "stage", "grade"),
#'     model_type = "coxph",
#'     labels = clintrial_labels
#' )
#' print(cox_screen)
#' # Returns hazard ratios (HR) instead of odds ratios
#' 
#' # Example 5: Keep models for diagnostics
#' screen5 <- uscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "bmi", "creatinine"),
#'     keep_models = TRUE
#' )
#' 
#' # Access stored models
#' models <- attr(screen5, "models")
#' summary(models[["age"]])
#' plot(models[["age"]])  # Diagnostic plots
#' 
#' # Example 6: Linear regression screening
#' linear_screen <- uscreen(
#'     data = clintrial,
#'     outcome = "bmi",
#'     predictors = c("age", "sex", "smoking", "creatinine", "hemoglobin"),
#'     model_type = "lm",
#'     labels = clintrial_labels
#' )
#' print(linear_screen)
#' 
#' # Example 7: Poisson regression for count outcomes
#' poisson_screen <- uscreen(
#'     data = clintrial,
#'     outcome = "los_days",
#'     predictors = c("age", "sex", "treatment", "surgery", "stage"),
#'     model_type = "glm",
#'     family = "poisson",
#'     labels = clintrial_labels
#' )
#' print(poisson_screen)
#' # Returns rate ratios (RR)
#' 
#' # Example 8: Hide reference rows for factor variables
#' screen8 <- uscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("treatment", "stage", "grade"),
#'     add_reference_rows = FALSE
#' )
#' print(screen8)
#' # Reference categories not shown
#' 
#' # Example 9: Customize decimal places
#' screen9 <- uscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "bmi", "creatinine"),
#'     digits = 3,      # 3 decimal places for OR
#'     digits_p = 4     # 4 decimal places for p-values
#' )
#' print(screen9)
#' 
#' # Example 10: Hide sample size and event columns
#' screen10 <- uscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi"),
#'     show_n = FALSE,
#'     show_events = FALSE
#' )
#' print(screen10)
#' 
#' # Example 11: Access raw numeric data
#' screen11 <- uscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment")
#' )
#' raw_data <- attr(screen11, "raw_data")
#' print(raw_data)
#' # Contains unformatted coefficients, SEs, CIs, etc.
#' 
#' # Example 12: Force coefficient display instead of OR
#' screen12 <- uscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "bmi"),
#'     model_type = "glm",
#'     family = "binomial",
#'     exponentiate = FALSE  # Show log odds instead of OR
#' )
#' print(screen12)
#' 
#' # Example 13: Screening with weights
#' screen13 <- uscreen(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "bmi"),
#'     model_type = "coxph",
#'     weights = runif(nrow(clintrial), min = 0.5, max = 2)  # Random numbers for example
#' )
#' 
#' # Example 14: Strict significance filter (p < 0.05)
#' sig_only <- uscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "hypertension", 
#'                   "diabetes", "ecog", "treatment", "stage", "grade"),
#'     p_threshold = 0.05,
#'     labels = clintrial_labels
#' )
#' 
#' # Check how many predictors passed the filter
#' n_significant <- length(unique(sig_only$Variable[sig_only$Variable != ""]))
#' cat("Significant predictors:", n_significant, "\n")
#' 
#' # Example 15: Complete workflow - screen then use in multivariable
#' # Step 1: Screen with liberal threshold
#' candidates <- uscreen(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "hypertension",
#'                   "diabetes", "treatment", "stage", "grade"),
#'     p_threshold = 0.20
#' )
#' 
#' # Step 2: Extract significant predictor names from raw data
#' sig_predictors <- unique(attr(candidates, "raw_data")$variable)
#' 
#' # Step 3: Fit multivariable model with selected predictors
#' multi_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = sig_predictors,
#'     labels = clintrial_labels
#' )
#' print(multi_model)
#'
#' @export
uscreen <- function(data,
                    outcome,
                    predictors,
                    model_type = "glm",
                    family = "binomial",
                    p_threshold = 1,
                    conf_level = 0.95,
                    add_reference_rows = TRUE,
                    show_n = TRUE,
                    show_events = TRUE,
                    digits = 2,
                    digits_p = 3,
                    labels = NULL,
                    keep_models = FALSE,
                    exponentiate = NULL,
                    ...) {
    
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
    ## Store models only if requested
    if (keep_models) {
        models <- list()
    }
    raw_results <- list()
    
    ## Fit univariable model for each predictor
    for (pred in predictors) {
        ## Build formula
        formula_str <- paste(outcome, "~", pred)
        formula <- stats::as.formula(formula_str)
        
        ## Fit model based on type
        if (model_type == "glm") {
            model <- stats::glm(formula, data = data, family = family, ...)
        } else if (model_type == "lm") {
            model <- stats::lm(formula, data = data, ...)
        } else if (model_type == "coxph") {
            if (!requireNamespace("survival", quietly = TRUE)) 
                stop("Package 'survival' required for Cox models")
            model <- survival::coxph(formula, data = data, ...)
        } else if (model_type == "clogit") {
            if (!requireNamespace("survival", quietly = TRUE))
                stop("Package 'survival' required for conditional logistic regression")
            model <- survival::clogit(formula, data = data, ...)
        } else {
            stop("Unsupported model type: ", model_type)
        }

        ## Store the data directly in the model
        model$data <- data
        
        ## Store the model only if requested
        if (keep_models) {
            models[[pred]] <- model
        }
        
        ## Get raw results using m2dt
        raw_result <- m2dt(data = data,
                           model = model,
                           conf_level = conf_level,
                           keep_qc_stats = FALSE,
                           include_intercept = FALSE,
                           add_reference_rows = add_reference_rows)
        
        ## Add predictor name for tracking
        raw_result[, predictor := pred]
        raw_results[[pred]] <- raw_result
        
        ## Clear model from memory if not keeping
        if (!keep_models) {
            rm(model)
        }
    }
    
    ## Combine all raw results
    combined_raw <- rbindlist(raw_results, fill = TRUE)
    
    ## Filter by p-value if requested
    if (p_threshold < 1) {
        passing_predictors <- unique(combined_raw[p_value <= p_threshold]$predictor)
        combined_raw <- combined_raw[predictor %in% passing_predictors]
        
        ## If filtering and keeping models, remove non-passing models
        if (keep_models) {
            models <- models[names(models) %in% passing_predictors]
        }
    }
    
    ## Format the combined results
    formatted <- format_model_table(combined_raw,
                                    show_n = show_n,
                                    show_events = show_events,
                                    digits = digits,
                                    digits_p = digits_p,
                                    labels = labels,
                                    exponentiate = exponentiate
                                    )
    
    ## Attach raw data
    setattr(formatted, "raw_data", combined_raw)
    
    ## Attach models only if kept
    if (keep_models) {
        setattr(formatted, "models", models)
    }
    
    setattr(formatted, "outcome", outcome)
    setattr(formatted, "model_type", unique(combined_raw$model_type)[1])
    setattr(formatted, "model_scope", "Univariable")
    setattr(formatted, "screening_type", "univariable")
    
    ## Add class for methods
    class(formatted) <- c("uscreen_result", class(formatted))
    
    return(formatted)
}

#' Print method for uscreen results
#' @keywords internal
#' @export
print.uscreen_result <- function(x, ...) {
    cat("\nUnivariable Screening Results\n")
    cat("Outcome: ", attr(x, "outcome"), "\n", sep = "")
    cat("Model Type: ", attr(x, "model_type"), "\n", sep = "")
    
    ## Get unique predictors from the formatted table
    n_predictors <- length(unique(x$Variable[x$Variable != ""]))
    cat("Predictors Screened: ", n_predictors, "\n", sep = "")
    
    ## Count significant from raw data
    raw <- attr(x, "raw_data")
    if (!is.null(raw) && "p_value" %in% names(raw)) {
        sig_predictors <- unique(raw[p_value < 0.05]$predictor)
        n_sig <- length(sig_predictors)
        cat("Significant (p < 0.05): ", n_sig, "\n", sep = "")
    }
    
    ## Note if models are stored
    if (!is.null(attr(x, "models"))) {
        cat("Models stored: Yes (", length(attr(x, "models")), ")\n", sep = "")
    }
    
    cat("\n")
    
    NextMethod("print", x)
    invisible(x)
}
