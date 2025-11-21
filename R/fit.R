#' Fit Regression Model with Publication-Ready Output
#'
#' Provides a unified interface for fitting various types of regression models 
#' with automatic formatting of results for publication. Supports generalized 
#' linear models, linear models, survival models, and mixed effects models with 
#' consistent syntax and output formatting. Handles both univariable and 
#' multivariable models automatically.
#'
#' @param data A data.frame or data.table containing the analysis dataset. The 
#'   function automatically converts data.frames to data.tables for efficient 
#'   processing.
#'   
#' @param outcome Character string specifying the outcome variable name. For 
#'   survival analysis, use \code{Surv()} syntax from the survival package 
#'   (e.g., \code{"Surv(time, status)"} or \code{"Surv(os_months, os_status)"}).
#'   
#' @param predictors Character vector of predictor variable names to include in 
#'   the model. All predictors are included simultaneously (multivariable model). 
#'   For univariable models, provide a single predictor. Can include continuous, 
#'   categorical (factor), or binary variables.
#'   
#' @param model_type Character string specifying the type of regression model. 
#'   Options include:
#'   \itemize{
#'     \item \code{"glm"} - Generalized linear model [default]
#'     \item \code{"lm"} - Linear regression
#'     \item \code{"coxph"} - Cox proportional hazards (survival)
#'     \item \code{"clogit"} - Conditional logistic regression
#'     \item \code{"coxme"} - Mixed effects Cox model
#'     \item \code{"glmer"} - Generalized linear mixed effects model
#'   }
#'   
#' @param family For GLM models, the error distribution and link function. Common 
#'   options include:
#'   \itemize{
#'     \item \code{"binomial"} - Logistic regression [default]
#'     \item \code{"poisson"} - Poisson regression for count data
#'     \item \code{"gaussian"} - Normal linear regression via GLM
#'     \item \code{"Gamma"} - Gamma regression for positive continuous data
#'   }
#'   See \code{\link[stats]{family}} for all options. Ignored for non-GLM models.
#'   
#' @param interactions Character vector of interaction terms using colon 
#'   notation (e.g., \code{c("age:sex", "treatment:stage")}). Interaction terms 
#'   are added to the model in addition to main effects. Default is \code{NULL} 
#'   (no interactions).
#'   
#' @param strata For Cox or conditional logistic models, character string naming 
#'   the stratification variable. Creates separate baseline hazards for each 
#'   stratum level without estimating stratum effects. Default is \code{NULL}.
#'   
#' @param cluster For Cox models, character string naming the variable for 
#'   robust clustered standard errors. Accounts for within-cluster correlation 
#'   (e.g., patients within hospitals). Default is \code{NULL}.
#'   
#' @param weights Character string naming the weights variable in \code{data}. 
#'   The specified column should contain numeric values for observation weights. 
#'   Used for weighted regression, survey data, or inverse probability weighting. 
#'   Default is \code{NULL}.
#'   
#' @param conf_level Numeric confidence level for confidence intervals. Must be 
#'   between 0 and 1. Default is 0.95 (95\% confidence intervals).
#'   
#' @param add_reference_rows Logical. If \code{TRUE}, adds rows for reference 
#'   categories of factor variables with baseline values (OR/HR/RR = 1, 
#'   coefficient = 0). Default is \code{TRUE}.
#'   
#' @param show_n Logical. If \code{TRUE}, includes the sample size column in 
#'   the output. Default is \code{TRUE}.
#'   
#' @param show_events Logical. If \code{TRUE}, includes the events column in 
#'   the output (for survival and logistic regression). Default is \code{TRUE}.
#'   
#' @param digits Integer specifying the number of decimal places for effect 
#'   estimates (OR, HR, RR, coefficients). Default is 2.
#'   
#' @param digits_p Integer specifying the number of decimal places for p-values. 
#'   P-values smaller than \code{10^(-digits_p)} are displayed as "< 0.001", 
#'   etc. Default is 3.
#'   
#' @param labels Named character vector or list providing custom display 
#'   labels for variables. Names should match variable names, values are display 
#'   labels. Default is \code{NULL}.
#'   
#' @param keep_qc_stats Logical. If \code{TRUE}, includes model quality statistics 
#'   (AIC, BIC, R-squared, concordance, etc.) in the raw data attribute for 
#'   model diagnostics and comparison. Default is \code{TRUE}.
#'   
#' @param exponentiate Logical. Whether to exponentiate coefficients. Default 
#'   is \code{NULL}, which automatically exponentiates for logistic, Poisson, 
#'   and Cox models, and displays raw coefficients for linear models. Set to 
#'   \code{TRUE} to force exponentiation or \code{FALSE} to force coefficients.
#'   
#' @param ... Additional arguments passed to the underlying model fitting 
#'   function (\code{\link[stats]{glm}}, \code{\link[stats]{lm}}, 
#'   \code{\link[survival]{coxph}}, \code{\link[lme4]{glmer}}, etc.). Common 
#'   options include \code{subset}, \code{na.action}, and model-specific control 
#'   parameters.
#'
#' @return A data.table with S3 class \code{"fit_result"} containing formatted 
#'   regression results. The table structure includes:
#'   \describe{
#'     \item{Variable}{Character. Predictor name or custom label}
#'     \item{Group}{Character. For factor variables: category level. For 
#'       interactions: interaction term. For continuous: typically empty}
#'     \item{n}{Integer. Total sample size (if \code{show_n = TRUE})}
#'     \item{n_group}{Integer. Sample size for this factor level}
#'     \item{events}{Integer. Total number of events (if \code{show_events = TRUE})}
#'     \item{events_group}{Integer. Events for this factor level}
#'     \item{Univariable/Multivariable OR/HR/RR/Estimate (95\% CI)}{Character. 
#'       Formatted effect estimate with confidence interval. Label depends on 
#'       model scope (univariable vs multivariable) and type. Multivariable 
#'       models use adjusted notation: aOR, aHR, aRR}
#'     \item{p-value}{Character. Formatted p-value from Wald test}
#'   }
#'   
#'   The returned object includes the following attributes accessible via \code{attr()}:
#'   \describe{
#'     \item{model}{The fitted model object (glm, lm, coxph, etc.). Access for 
#'       diagnostics, predictions, or further analysis}
#'     \item{raw_data}{data.table. Unformatted numeric results with columns for 
#'       coefficients, standard errors, confidence bounds, quality statistics, etc.}
#'     \item{outcome}{Character. The outcome variable name}
#'     \item{predictors}{Character vector. The predictor variable names}
#'     \item{formula_str}{Character. The complete model formula as a string}
#'     \item{model_scope}{Character. "Univariable" (one predictor) or 
#'       "Multivariable" (multiple predictors)}
#'     \item{model_type}{Character. The regression model type used}
#'     \item{interactions}{Character vector (if interactions specified). 
#'       The interaction terms included}
#'     \item{strata}{Character (if stratification used). The stratification variable}
#'     \item{cluster}{Character (if clustering used). The cluster variable}
#'     \item{weights}{Character (if weighting used). The weights variable}
#'   }
#'
#' @details
#' \strong{Model Scope Detection:}
#' 
#' The function automatically detects whether the model is:
#' \itemize{
#'   \item \strong{Univariable}: Single predictor (e.g., \code{predictors = "age"})
#'     - Effect estimates labeled as "Univariable OR", "Univariable HR", etc.
#'     - Represents crude (unadjusted) association
#'   \item \strong{Multivariable}: Multiple predictors (e.g., 
#'     \code{predictors = c("age", "sex", "treatment")})
#'     - Effect estimates labeled as "Multivariable aOR", "Multivariable aHR", etc.
#'     - "a" prefix indicates "adjusted" for other variables in the model
#'     - Represents associations adjusted for confounding
#' }
#' 
#' \strong{Interaction Terms:}
#' 
#' Interactions are specified using colon notation and added to the model:
#' \itemize{
#'   \item \code{interactions = c("age:treatment")} creates interaction 
#'     between age and treatment
#'   \item Main effects for both variables are automatically included
#'   \item Multiple interactions can be specified: 
#'     \code{c("age:sex", "treatment:stage")}
#'   \item For interactions between categorical variables, separate terms are 
#'     created for each combination of levels
#' }
#' 
#' \strong{Stratification (Cox/Conditional Logistic):}
#' 
#' The \code{strata} parameter creates separate baseline hazards:
#' \itemize{
#'   \item Allows baseline hazard to vary across strata without estimating 
#'     stratum effects
#'   \item Useful when proportional hazards assumption violated across strata
#'   \item Example: \code{strata = "center"} for multicenter studies
#'   \item Stratification variable is not included as a predictor
#' }
#' 
#' \strong{Clustering (Cox Models):}
#' 
#' The \code{cluster} parameter computes robust standard errors:
#' \itemize{
#'   \item Accounts for within-cluster correlation (e.g., multiple observations 
#'     per patient)
#'   \item Uses sandwich variance estimator
#'   \item Does not change point estimates, only standard errors and p-values
#' }
#' 
#' \strong{Weighting:}
#' 
#' The \code{weights} parameter enables weighted regression:
#' \itemize{
#'   \item For survey data with sampling weights
#'   \item Inverse probability weighting for causal inference
#'   \item Frequency weights for aggregated data
#'   \item Weights should be in a column of \code{data}
#' }
#' 
#' \strong{Effect Measures by Model Type:}
#' \itemize{
#'   \item \strong{Logistic} (\code{family = "binomial"}): Odds ratios (OR/aOR)
#'   \item \strong{Cox} (\code{model_type = "coxph"}): Hazard ratios (HR/aHR)
#'   \item \strong{Poisson} (\code{family = "poisson"}): Rate ratios (RR/aRR)
#'   \item \strong{Linear}: Raw coefficient estimates
#' }
#'
#' @seealso 
#' \code{\link{uscreen}} for univariable screening of multiple predictors,
#' \code{\link{fastfit}} for complete univariable-to-multivariable workflow,
#' \code{\link{compfit}} for comparing multiple models,
#' \code{\link{m2dt}} for model-to-table conversion
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' library(survival)
#' 
#' # Example 1: Univariable logistic regression
#' uni_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = "age"
#' )
#' print(uni_model)
#' # Labeled as "Univariable OR"
#' 
#' # Example 2: Multivariable logistic regression
#' multi_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "treatment"),
#'     labels = clintrial_labels
#' )
#' print(multi_model)
#' # Labeled as "Multivariable aOR" (adjusted OR)
#' 
#' # Example 3: Cox proportional hazards model
#' cox_model <- fit(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     model_type = "coxph",
#'     labels = clintrial_labels
#' )
#' print(cox_model)
#' 
#' # Example 4: Model with interaction terms
#' interact_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "treatment", "sex"),
#'     interactions = c("age:treatment"),
#'     labels = clintrial_labels
#' )
#' print(interact_model)
#' 
#' # Example 5: Cox model with stratification
#' strat_model <- fit(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "treatment"),
#'     model_type = "coxph",
#'     strata = "site",  # Separate baseline hazards by site
#'     labels = clintrial_labels
#' )
#' print(strat_model)
#' 
#' # Example 6: Cox model with clustering
#' cluster_model <- fit(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "treatment"),
#'     model_type = "coxph",
#'     cluster = "site",  # Robust SEs accounting for site clustering
#'     labels = clintrial_labels
#' )
#' print(cluster_model)
#' 
#' # Example 7: Linear regression
#' linear_model <- fit(
#'     data = clintrial,
#'     outcome = "bmi",
#'     predictors = c("age", "sex", "smoking"),
#'     model_type = "lm",
#'     labels = clintrial_labels
#' )
#' print(linear_model)
#' 
#' # Example 8: Poisson regression for count data
#' poisson_model <- fit(
#'     data = clintrial,
#'     outcome = "los_days",
#'     predictors = c("age", "treatment", "surgery", "stage"),
#'     model_type = "glm",
#'     family = "poisson",
#'     labels = clintrial_labels
#' )
#' print(poisson_model)
#' # Returns rate ratios (RR/aRR)
#' 
#' # Example 9: Access the underlying fitted model
#' result <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi")
#' )
#' 
#' # Get the model object
#' model_obj <- attr(result, "model")
#' summary(model_obj)
#' 
#' # Model diagnostics
#' plot(model_obj)
#' 
#' # Predictions
#' preds <- predict(model_obj, type = "response")
#' 
#' # Example 10: Access raw numeric data
#' raw_data <- attr(result, "raw_data")
#' print(raw_data)
#' # Contains unformatted coefficients, SEs, CIs, AIC, BIC, etc.
#' 
#' # Example 11: Multiple interactions
#' complex_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment", "bmi"),
#'     interactions = c("age:treatment", "sex:bmi"),
#'     labels = clintrial_labels
#' )
#' print(complex_model)
#' 
#' # Example 12: Customize output columns
#' minimal <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment"),
#'     show_n = FALSE,
#'     show_events = FALSE,
#'     add_reference_rows = FALSE
#' )
#' print(minimal)
#' 
#' # Example 13: Different confidence levels
#' ci90 <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "treatment"),
#'     conf_level = 0.90  # 90% confidence intervals
#' )
#' print(ci90)
#' 
#' # Example 14: Force coefficient display instead of OR
#' coef_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "bmi"),
#'     exponentiate = FALSE  # Show log odds instead of OR
#' )
#' print(coef_model)
#' 
#' # Example 15: Check model quality statistics
#' result <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     keep_qc_stats = TRUE
#' )
#' 
#' raw <- attr(result, "raw_data")
#' cat("AIC:", raw$AIC[1], "\n")
#' cat("BIC:", raw$BIC[1], "\n")
#' cat("C-statistic:", raw$c_statistic[1], "\n")
#'
#' @export
fit <- function(data, 
                outcome,
                predictors,
                model_type = "glm",
                family = "binomial",
                interactions = NULL,
                strata = NULL,
                cluster = NULL,
                weights = NULL,
                conf_level = 0.95,
                add_reference_rows = TRUE,
                show_n = TRUE,
                show_events = TRUE,
                digits = 2,
                digits_p = 3,
                labels = NULL,
                keep_qc_stats = TRUE,
                exponentiate = NULL,
                ...) {
    
    ## Don't create internal variables - work with 'data' directly
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
    ## Store original data name for reference
    data_name <- deparse(substitute(data))
    
    ## Build formula
    if (!is.null(interactions)) {
        formula_str <- paste(outcome, "~", 
                             paste(c(predictors, interactions), collapse = " + "))
    } else {
        formula_str <- paste(outcome, "~", paste(predictors, collapse = " + "))
    }
    
    ## Add strata if provided (for Cox models)
    if (!is.null(strata) && model_type %in% c("coxph", "clogit")) {
        formula_str <- paste(formula_str, "+ strata(", strata, ")")
    }
    
    formula <- stats::as.formula(formula_str)
    
    ## Fit model based on type - use 'data' directly
    if (model_type == "glm") {
        if (!is.null(weights)) {
            model <- stats::glm(formula, data = data, family = family, 
                                weights = data[[weights]], ...)
        } else {
            model <- stats::glm(formula, data = data, family = family, ...)
        }
        
    } else if (model_type == "lm") {
        if (!is.null(weights)) {
            model <- stats::lm(formula, data = data, weights = data[[weights]], ...)
        } else {
            model <- stats::lm(formula, data = data, ...)
        }
        
    } else if (model_type == "coxph") {
        if (!requireNamespace("survival", quietly = TRUE)) 
            stop("Package 'survival' required for Cox models")
        
        if (!is.null(cluster)) {
            model <- survival::coxph(formula, data = data, 
                                     cluster = data[[cluster]], ...)
        } else {
            model <- survival::coxph(formula, data = data, ...)
        }
        
    } else if (model_type == "clogit") {
        if (!requireNamespace("survival", quietly = TRUE)) 
            stop("Package 'survival' required for conditional logistic regression")
        model <- survival::clogit(formula, data = data, ...)
        
    } else if (model_type == "coxme") {
        if (!requireNamespace("coxme", quietly = TRUE))
            stop("Package 'coxme' required for mixed effects Cox models")
        model <- coxme::coxme(formula, data = data, ...)
        
    } else if (model_type == "glmer") {
        if (!requireNamespace("lme4", quietly = TRUE))
            stop("Package 'lme4' required for mixed effects models")
        model <- lme4::glmer(formula, data = data, family = family, ...)
        
    } else {
        stop("Unsupported model type: ", model_type)
    }
    
    ## Store the data directly in the model
    if (isS4(model)) {
        attr(model, "data") <- data
    } else {
        model$data <- data
        attr(model, "data") <- data
    }
    
    ## Attach metadata as attributes
    data.table::setattr(model, "formula_str", formula_str)
    data.table::setattr(model, "predictors", predictors)
    data.table::setattr(model, "model_type", model_type)
    data.table::setattr(model, "data_name", data_name)
    
    if (!is.null(interactions)) {
        data.table::setattr(model, "interactions", interactions)
    }
    if (!is.null(strata)) {
        data.table::setattr(model, "strata", strata)
    }
    if (!is.null(cluster)) {
        data.table::setattr(model, "cluster", cluster)
    }
    if (!is.null(weights)) {
        data.table::setattr(model, "weights", weights)
    }
    
    ## Convert to readable format using m2dt
    raw_data <- m2dt(data,
                     model, 
                     conf_level = conf_level,
                     keep_qc_stats = keep_qc_stats,
                     include_intercept = FALSE,
                     add_reference_rows = add_reference_rows)

    ## Format results for publication
    formatted <- format_model_table(raw_data,
                                    show_n = show_n,
                                    show_events = show_events,
                                    digits = digits,
                                    digits_p = digits_p,
                                    labels = labels,
                                    exponentiate = exponentiate)

    ## Attach metadata
    data.table::setattr(formatted, "model", model)
    data.table::setattr(formatted, "raw_data", raw_data)
    data.table::setattr(formatted, "outcome", outcome)
    data.table::setattr(formatted, "predictors", predictors)
    data.table::setattr(formatted, "formula_str", formula_str)
    data.table::setattr(formatted, "model_scope", raw_data$model_scope[1])
    data.table::setattr(formatted, "model_type", raw_data$model_type[1])
    
    ## Add metadata from model fitting
    if (!is.null(interactions)) {
        setattr(formatted, "interactions", interactions)
    }
    if (!is.null(strata)) {
        setattr(formatted, "strata", strata)
    }
    if (!is.null(cluster)) {
        setattr(formatted, "cluster", cluster)
    }
    if (!is.null(weights)) {
        setattr(formatted, "weights", weights)
    }

    class(formatted) <- c("fit_result", class(formatted))
    
    return(formatted)
}

#' Print method for fit results
#' @keywords internal
#' @export
print.fit_result <- function(x, ...) {
    cat("\n", attr(x, "model_scope"), " ", attr(x, "model_type"), " Model\n", sep = "")
    cat("Formula: ", attr(x, "formula_str"), "\n", sep = "")
    
    ## Get sample size from raw results
    raw <- attr(x, "raw_data")
    if (!is.null(raw)) {
        cat("N = ", raw$n[1], sep = "")
        if (!is.na(raw$events[1])) cat(", Events = ", raw$events[1], sep = "")
        cat("\n")
    }
    cat("\n")
    
    NextMethod("print", x)
    invisible(x)
}
