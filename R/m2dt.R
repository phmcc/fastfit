#' Convert Model to Data Table
#'
#' Extracts coefficients, confidence intervals, and comprehensive model statistics 
#' from fitted regression models and converts them to a standardized data.table 
#' format suitable for further analysis or publication. This is a core utility 
#' function in the fastfit package for working with regression results.
#'
#' @param model A fitted model object. Supported classes include:
#'   \itemize{
#'     \item \code{glm} - Generalized linear models (logistic, Poisson, etc.)
#'     \item \code{lm} - Linear models
#'     \item \code{coxph} - Cox proportional hazards models
#'     \item \code{clogit} - Conditional logistic regression
#'     \item \code{coxme} - Mixed effects Cox models
#'     \item \code{glmer} - Generalized linear mixed effects models
#'   }
#'   Also accepts models wrapped with \code{mmodel()}.
#'   
#' @param conf_level Numeric confidence level for confidence intervals. Must be
#'   between 0 and 1. Default is 0.95 (95\% CI).
#'   
#' @param keep_qc_stats Logical. If \code{TRUE}, includes model quality statistics 
#'   such as AIC, BIC, R-squared, concordance, and model fit tests. These appear 
#'   as additional columns in the output. Default is \code{TRUE}.
#'   
#' @param include_intercept Logical. If \code{TRUE}, includes the model intercept 
#'   in output. If \code{FALSE}, removes the intercept row from results. Useful 
#'   for creating cleaner presentation tables. Default is \code{TRUE}.
#'   
#' @param terms_to_exclude Character vector of term names to exclude from output.
#'   Useful for removing specific unwanted parameters (e.g., nuisance variables,
#'   spline terms). Default is \code{NULL}. Note: If \code{include_intercept = FALSE}, 
#'   "(Intercept)" is automatically added to this list.
#'   
#' @param add_reference_rows Logical. If \code{TRUE}, adds rows for reference 
#'   categories of factor variables with appropriate labels and baseline values 
#'   (OR/HR = 1, Estimate = 0). This makes tables more complete and easier to 
#'   interpret. Default is \code{TRUE}.
#'   
#' @param reference_label Character string used to label reference category rows
#'   in the output. Appears in the \code{reference} column. Default is \code{"reference"}.
#'
#' @return A \code{data.table} containing extracted model information with the 
#'   following standard columns:
#'   \describe{
#'     \item{model_scope}{Character. Either "Univariable" (unadjusted model with 
#'       single predictor) or "Multivariable" (adjusted model with multiple predictors)}
#'     \item{model_type}{Character. Type of regression (e.g., "Logistic", "Linear", 
#'       "Cox PH", "Poisson")}
#'     \item{variable}{Character. Variable name (for factor variables, the base 
#'       variable name without the level)}
#'     \item{group}{Character. Group/level name for factor variables; empty string 
#'       for continuous variables}
#'     \item{n}{Integer. Total sample size used in the model}
#'     \item{n_group}{Integer. Sample size for this specific variable level 
#'       (factor variables only)}
#'     \item{events}{Integer. Total number of events in the model (for survival 
#'       and logistic models)}
#'     \item{events_group}{Integer. Number of events for this specific variable 
#'       level (for survival and logistic models with factor variables)}
#'     \item{coefficient}{Numeric. Raw regression coefficient (log odds, log hazard, 
#'       etc.)}
#'     \item{se}{Numeric. Standard error of the coefficient}
#'     \item{OR/HR/RR/Estimate}{Numeric. Effect estimate - column name depends on 
#'       model type:
#'       \itemize{
#'         \item \code{OR} for logistic regression (odds ratio)
#'         \item \code{HR} for Cox models (hazard ratio)
#'         \item \code{RR} for Poisson regression (rate/risk ratio)
#'         \item \code{Estimate} for linear models or other GLMs
#'       }}
#'     \item{CI_lower}{Numeric. Lower bound of confidence interval for effect estimate}
#'     \item{CI_upper}{Numeric. Upper bound of confidence interval for effect estimate}
#'     \item{statistic}{Numeric. Test statistic (z-value for GLM/Cox, t-value for LM)}
#'     \item{p_value}{Numeric. P-value for coefficient test}
#'     \item{sig}{Character. Significance markers: "***" (p<0.001), "**" (p<0.01), 
#'       "*" (p<0.05), "." (p<0.10), "" (p≥0.10)}
#'     \item{sig_binary}{Logical. Binary indicator: \code{TRUE} if p<0.05, 
#'       \code{FALSE} otherwise}
#'     \item{reference}{Character. Contains \code{reference_label} for reference 
#'       category rows when \code{add_reference_rows = TRUE}, empty string otherwise}
#'   }
#'
#' @export
m2dt <- function(data,
                 model,
                 conf_level = 0.95,
                 keep_qc_stats = TRUE,
                 include_intercept = TRUE,
                 terms_to_exclude = NULL,
                 add_reference_rows = TRUE,
                 reference_label = "reference") {

    ## Validate inputs
    if (missing(data)) {
        stop("data parameter is required. Usage: m2dt(data, model)")
    }
    
    if (!inherits(data, c("data.frame", "data.table"))) {
        stop("data must be a data.frame or data.table")
    }
    
    ## Store data for use throughout function
    model_data <- data.table::as.data.table(data)
    
    ## Store as attribute for helper functions
    attr(model, "data") <- model_data

    ## Set model class
    model_class <- class(model)[1]
    
    ## Null operator for handling missing values
    `%||%` <- function(a, b) if (is.null(a)) b else a
    
    ## Add intercept to exclusion list if requested
    if (!include_intercept) {
        terms_to_exclude <- unique(c(terms_to_exclude, "(Intercept)"))
    }
    
    ## Get model type (Univariable vs Multivariable)
    model_scope <- detect_model_type(model)
    
    ## Get readable model type name
    model_type_name <- get_model_type_name(model)
    
    ## Extract results based on model class
    if (model_class %in% c("glm", "lm")) {
        
        coef_summary <- summary(model)$coefficients
        conf_int <- stats::confint.default(model, level = conf_level)
        
        dt <- data.table::data.table(
                              model_scope = model_scope,
                              model_type = model_type_name,
                              term = rownames(coef_summary),
                              n = stats::nobs(model),
                              events = NA_integer_,
                              coefficient = coef_summary[, "Estimate"],
                              se = coef_summary[, "Std. Error"]
                          )
        
        ## Calculate confidence intervals
        z_score <- stats::qnorm((1 + conf_level) / 2)
        dt[, `:=`(
            coef = coefficient,
            coef_lower = coefficient - z_score * se,
            coef_upper = coefficient + z_score * se
        )]
        
        ## Special handling for logistic regression
        if (model_class == "glm" && !isS4(model)) {
            if (model$family$family == "binomial") {
                if (!is.null(model$y)) {
                    dt[, events := sum(model$y, na.rm = TRUE)]
                } else if (!is.null(model$model)) {
                    outcome_col <- model$model[[1]]
                    if (is.factor(outcome_col)) {
                        dt[, events := sum(as.numeric(outcome_col) == 2, na.rm = TRUE)]
                    } else {
                        dt[, events := sum(outcome_col, na.rm = TRUE)]
                    }
                }
            }
        }
        
        ## Determine if should exponentiate
        should_exp <- FALSE
        is_logistic <- FALSE
        is_poisson <- FALSE
        
        if (model_class == "glm" && !isS4(model)) {
            family_name <- model$family$family
            link_name <- model$family$link
            
            is_logistic <- family_name == "binomial"
            is_poisson <- family_name == "poisson" || 
                (family_name == "quasipoisson")
            
            should_exp <- is_logistic || is_poisson || link_name == "log"
        }
        
        ## Add exponentiated coefficients
        dt[, `:=`(
            exp_coef = if (should_exp) exp(coefficient) else coefficient,
            exp_lower = if (should_exp) exp(coef_lower) else coef_lower,
            exp_upper = if (should_exp) exp(coef_upper) else coef_upper
        )]
        
        if (is_logistic) {
            dt[, `:=`(
                OR = exp_coef,
                CI_lower = exp_lower,
                CI_upper = exp_upper
            )]
        } else if (is_poisson) {
            dt[, `:=`(
                RR = exp_coef,
                CI_lower = exp_lower,
                CI_upper = exp_upper
            )]
        } else if (should_exp) {
            ## Other log-link models
            dt[, `:=`(
                Estimate = exp_coef,
                CI_lower = exp_lower,
                CI_upper = exp_upper
            )]
        } else {
            ## Linear models - use raw coefficients
            dt[, `:=`(
                Estimate = coef,
                CI_lower = coef_lower,
                CI_upper = coef_upper
            )]
        }
        
        ## Add test statistics
        stat_col <- if ("z value" %in% colnames(coef_summary)) "z value" else "t value"
        dt[, `:=`(
            statistic = coef_summary[, stat_col],
            p_value = coef_summary[, ncol(coef_summary)]
        )]
        
        ## Add QC stats if requested
        if (keep_qc_stats) {
            if (model_class == "glm") {
                dt[, `:=`(
                    AIC = stats::AIC(model),
                    BIC = stats::BIC(model),
                    deviance = stats::deviance(model),
                    null_deviance = if (!isS4(model)) model$null.deviance else NA,
                    df_residual = stats::df.residual(model)
                )]
                
                ## Add R-squared for non-binomial GLMs
                if (!isS4(model) && model$family$family != "binomial") {
                    dt[, R2 := 1 - (deviance / null_deviance)]
                }
                
                ## For binomial, add discrimination/calibration metrics
                if (is_logistic && keep_qc_stats) {
                    ## C-statistic (if pROC available)
                    if (requireNamespace("pROC", quietly = TRUE)) {
                        if (!isS4(model)) {
                            roc_obj <- pROC::roc(model$y, stats::fitted(model), quiet = TRUE)
                            dt[, c_statistic := as.numeric(pROC::auc(roc_obj))]
                        }
                    }
                    
                    ## Hosmer-Lemeshow test (if ResourceSelection available)
                    if (requireNamespace("ResourceSelection", quietly = TRUE)) {
                        if (!isS4(model)) {
                            hl <- ResourceSelection::hoslem.test(model$y, stats::fitted(model), g = 10)
                            dt[, `:=`(
                                hoslem_chi2 = hl$statistic,
                                hoslem_p = hl$p.value
                            )]
                        }
                    }
                }
            } else if (model_class == "lm") {
                summ <- summary(model)
                dt[, `:=`(
                    R2 = summ$r.squared,
                    adj_R2 = summ$adj.r.squared,
                    AIC = stats::AIC(model),
                    BIC = stats::BIC(model),
                    sigma = summ$sigma,
                    df_residual = stats::df.residual(model)
                )]
            }
        }
        
    } else if (model_class %in% c("coxph", "clogit")) {
        
        if (!requireNamespace("survival", quietly = TRUE))
            stop("Package 'survival' required")
        
        summ <- summary(model)
        coef_summary <- summ$coefficients
        conf_int <- stats::confint(model, level = conf_level)
        
        dt <- data.table::data.table(
                              model_scope = model_scope %||% "Multivariable",
                              model_type = model_type_name,
                              term = rownames(coef_summary),
                              n = if (!is.null(model$n)) model$n[1] else summ$n,
                              events = if (!is.null(model$nevent)) model$nevent 
                                       else if (!is.null(model$n)) model$n[2] 
                                       else summ$nevent,
                              coefficient = coef_summary[, "coef"],
                              se = coef_summary[, "se(coef)"],
                              ## Store both versions
                              coef = coef_summary[, "coef"],
                              coef_lower = conf_int[, 1],
                              coef_upper = conf_int[, 2],
                              exp_coef = coef_summary[, "exp(coef)"],
                              exp_lower = exp(conf_int[, 1]),
                              exp_upper = exp(conf_int[, 2]),
                              ## Primary display columns
                              HR = coef_summary[, "exp(coef)"],
                              CI_lower = exp(conf_int[, 1]),
                              CI_upper = exp(conf_int[, 2]),
                              statistic = coef_summary[, "z"],
                              p_value = coef_summary[, "Pr(>|z|)"]
                          )

        ## Add QC stats
        if (keep_qc_stats) {
            dt[, `:=`(
                concordance = summ$concordance[1],
                concordance_se = summ$concordance[2],
                rsq = summ$rsq[1],
                rsq_max = summ$rsq[2],
                likelihood_ratio_test = summ$logtest[1],
                likelihood_ratio_df = summ$logtest[2],
                likelihood_ratio_p = summ$logtest[3],
                wald_test = summ$waldtest[1],
                wald_df = summ$waldtest[2],
                wald_p = summ$waldtest[3],
                score_test = summ$sctest[1],
                score_df = summ$sctest[2],
                score_p = summ$sctest[3]
            )]
        }

    } else if (model_class %in% c("coxme", "lme", "lmer", "glmer", "glmerMod")) {
        
        ## Mixed effects models
        summ <- summary(model)
        
        if (model_class == "coxme") {
            coef_vec <- coxme::fixef(model)
            vcov_mat <- as.matrix(stats::vcov(model))
            se_vec <- sqrt(diag(vcov_mat))
            z_vec <- coef_vec / se_vec
            p_vec <- 2 * (1 - stats::pnorm(abs(z_vec)))
            
            dt <- data.table::data.table(
                                  model_scope = model_scope %||% "Multivariable",
                                  model_type = model_type_name,
                                  term = names(coef_vec),
                                  n = model$n[2],
                                  events = model$n[1],
                                  n_group = NA_real_,
                                  events_group = NA_real_, 
                                  coefficient = coef_vec,
                                  se = se_vec,
                                  coef = coef_vec,
                                  coef_lower = coef_vec - stats::qnorm((1 + conf_level) / 2) * se_vec,
                                  coef_upper = coef_vec + stats::qnorm((1 + conf_level) / 2) * se_vec,
                                  exp_coef = exp(coef_vec),
                                  exp_lower = exp(coef_vec - stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  exp_upper = exp(coef_vec + stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  HR = exp(coef_vec),
                                  CI_lower = exp(coef_vec - stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  CI_upper = exp(coef_vec + stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  statistic = z_vec,
                                  p_value = p_vec
                              )
            
        } else {
            
            ## lmer/glmer from lme4
            if (!requireNamespace("lme4", quietly = TRUE))
                stop("Package 'lme4' required")
            
            coef_summary <- stats::coef(summ)
            
            ## Determine if should exponentiate
            should_exp <- (model_class %in% c("glmer", "glmerMod")) && 
                (summ$family == "binomial" || summ$link == "log")
            
            dt <- data.table::data.table(
                                  model_scope = model_scope %||% "Multivariable",
                                  model_type = model_type_name,
                                  term = rownames(coef_summary),
                                  n = stats::nobs(model),
                                  events = NA_integer_,
                                  coefficient = coef_summary[, "Estimate"],
                                  se = coef_summary[, "Std. Error"]
                              )
            
            ## For glmer binomial, calculate total events from the response
            if (model_class %in% c("glmer", "glmerMod") && inherits(model, "merMod")) {
                if (summ$family == "binomial") {
                    response_var <- model@resp$y  
                    if (!is.null(response_var)) {
                        dt[, events := sum(response_var, na.rm = TRUE)]
                    }
                }
            }
            
            z_score <- stats::qnorm((1 + conf_level) / 2)
            dt[, `:=`(
                coef = coefficient,
                coef_lower = coefficient - z_score * se,
                coef_upper = coefficient + z_score * se,
                exp_coef = if (should_exp) exp(coefficient) else coefficient,
                exp_lower = if (should_exp) exp(coefficient - z_score * se) else coefficient - z_score * se,
                exp_upper = if (should_exp) exp(coefficient + z_score * se) else coefficient + z_score * se
            )]
            
            ## Add appropriate effect column
            if (model_class %in% c("glmer", "glmerMod") && summ$family == "binomial") {
                dt[, `:=`(
                    OR = exp_coef,
                    CI_lower = exp_lower,
                    CI_upper = exp_upper
                )]
            } else if (should_exp) {
                dt[, `:=`(
                    RR = exp_coef,
                    CI_lower = exp_lower,
                    CI_upper = exp_upper
                )]
            } else {
                dt[, `:=`(
                    Estimate = coef,
                    CI_lower = coef_lower,
                    CI_upper = coef_upper
                )]
            }
            
            ## Add test statistics
            if (ncol(coef_summary) >= 3) {
                ## Use z-values if available
                stat_col <- if ("z value" %in% colnames(coef_summary)) {
                                "z value"
                            } else if ("t value" %in% colnames(coef_summary)) {
                                "t value"
                            } else {
                                NULL
                            }
                
                if (!is.null(stat_col)) {
                    dt[, statistic := coef_summary[, stat_col]]
                } else {
                    ## Calculate z-statistics manually
                    dt[, statistic := coefficient / se]
                }
                
                ## Check if p-values are provided
                p_col <- grep("^Pr\\(", colnames(coef_summary))
                if (length(p_col) > 0) {
                    dt[, p_value := coef_summary[, p_col[1]]]
                } else {
                    ## Calculate p-values from z/t statistics
                    dt[, p_value := 2 * (1 - stats::pnorm(abs(statistic)))]
                }
            } else {
                ## No test statistics - calculate manually
                dt[, `:=`(
                    statistic = coefficient / se,
                    p_value = 2 * (1 - stats::pnorm(abs(coefficient / se)))
                )]
            }
        }
        
    } else {
        stop("Unsupported model class: ", model_class)
    }
    
    ## Parse terms into variable and group AFTER creating dt
    xlevels <- get_model_xlevels(model)
    
    ## Parse terms - pass model for coxme special handling
    parsed <- parse_term(dt$term, xlevels, model)
    dt[, `:=`(variable = parsed$variable, group = parsed$group)]

    ## Prepare data for group counting
    if (!is.null(xlevels) || model_class == "coxme") {
        data_source <- model_data
        data_dt <- data.table::as.data.table(data_source)
        
        ## For coxme, reconstruct xlevels from the data if needed
        if (is.null(xlevels) && model_class == "coxme") {
            ## Use formulaList$fixed for coxme
            formula_to_use <- model$formulaList$fixed
            
            ## Get variables
            all_formula_vars <- all.vars(formula_to_use)
            term_vars <- all_formula_vars[-1]  # Exclude response
            
            xlevels <- list()
            
            for (var_name in term_vars) {
                if (var_name %in% names(data_dt) && is.factor(data_dt[[var_name]])) {
                    xlevels[[var_name]] <- levels(data_dt[[var_name]])
                }
            }
        }
        
        ## Determine outcome and event variables
        outcome_var <- NULL
        event_var <- NULL
        
        if (model_class == "glm" && !isS4(model)) {
            outcome_var <- all.vars(model$formula)[1]
        } else if (model_class %in% c("glmer", "glmerMod") && inherits(model, "merMod")) {
            ## For mixed models, get the response variable name from the formula
            formula_obj <- model@call$formula
            if (!is.null(formula_obj)) {
                outcome_var <- all.vars(formula_obj)[1]
            }
        } else if (model_class %in% c("coxph", "clogit", "coxme")) {
            ## For all survival models, use the unified helper function
            event_var <- get_event_variable(model, model_class)
        }
    
        ## Calculate n_group and events_group for all factor variables at once
        if (length(xlevels) > 0 && !is.null(data_dt)) {
            
            ## Get factor variables that exist in the data
            factor_vars <- names(xlevels)[names(xlevels) %in% names(data_dt)]
            
            if (length(factor_vars) > 0) {
                
                if (!is.null(event_var) && event_var %in% names(data_dt)) {
                    ## For survival models: calculate all counts at once
                    
                    ## Stack all factor variables into long format
                    all_counts <- rbindlist(lapply(factor_vars, function(var) {
                        data_dt[!is.na(get(var)), .(
                                                      variable = var,
                                                      group = as.character(get(var)),
                                                      n_group = .N,
                                                      events_group = sum(get(event_var), na.rm = TRUE)
                                                  ), by = get(var)][, get := NULL]
                    }))
                    
                    ## Single join to update all counts at once
                    dt[all_counts, `:=`(
                                       n_group = i.n_group,
                                       events_group = i.events_group
                                   ), on = .(variable, group)]
                    
                } else if (!is.null(outcome_var) && outcome_var %in% names(data_dt)) {
                    ## For GLM/GLMER models: calculate all counts at once
                    
                    ## Prepare outcome calculation
                    outcome_col <- data_dt[[outcome_var]]
                    if (is.factor(outcome_col)) {
                        data_dt[, .events_calc := as.numeric(outcome_col) == 2]
                    } else {
                        data_dt[, .events_calc := outcome_col]
                    }
                    
                    ## Stack all factor variables into long format
                    all_counts <- rbindlist(lapply(factor_vars, function(var) {
                        data_dt[!is.na(get(var)), .(
                                                      variable = var,
                                                      group = as.character(get(var)),
                                                      n_group = .N,
                                                      events_group = sum(.events_calc, na.rm = TRUE)
                                                  ), by = get(var)][, get := NULL]
                    }))
                    
                    ## Clean up temporary column
                    data_dt[, .events_calc := NULL]
                    
                    ## Single join to update all counts
                    dt[all_counts, `:=`(
                                       n_group = i.n_group,
                                       events_group = i.events_group
                                   ), on = .(variable, group)]
                    
                } else {
                    ## No outcome/event variable: just count n
                    all_counts <- rbindlist(lapply(factor_vars, function(var) {
                        data_dt[!is.na(get(var)), .(
                                                      variable = var,
                                                      group = as.character(get(var)),
                                                      n_group = .N
                                                  ), by = get(var)][, get := NULL]
                    }))
                    
                    ## Single join to update counts
                    dt[all_counts, n_group := i.n_group, on = .(variable, group)]
                }
            }
        }
    }
    
    ## Add reference rows for factor variables while maintaining original order
    xlevels_ref <- xlevels

    ## Create all reference rows at once (vectorized approach)
    if (!is.null(xlevels_ref) && add_reference_rows) {

        ## Add reference column to existing data
        dt[, reference := ""]
        
        ## Pre-calculate all reference level counts at once for efficiency
        ref_counts <- list()

        ## Pre-calculate all reference level counts from all_counts
        ref_counts <- list()

        ## Extract reference level counts from all_counts if it exists
        if (exists("all_counts") && !is.null(all_counts) && nrow(all_counts) > 0) {
            for (var in names(xlevels_ref)) {
                ref_level <- xlevels_ref[[var]][1]
                ref_row_data <- all_counts[variable == var & group == ref_level]
                if (nrow(ref_row_data) > 0) {
                    ref_counts[[var]] <- ref_row_data[, .(n_group = n_group, events_group = events_group)]
                }
            }
        } else {
            ## Fallback: calculate directly from data if all_counts doesn't exist
            if (!is.null(data_dt)) {
                for (var in names(xlevels_ref)) {
                    if (var %in% names(data_dt)) {
                        ref_level <- xlevels_ref[[var]][1]
                        
                        if (!is.null(event_var) && event_var %in% names(data_dt)) {
                            ref_counts[[var]] <- data_dt[get(var) == ref_level & !is.na(get(var)), .(
                                                                                                       n_group = .N,
                                                                                                       events_group = sum(get(event_var), na.rm = TRUE)
                                                                                                   )]
                        } else {
                            ref_counts[[var]] <- data_dt[get(var) == ref_level & !is.na(get(var)), .(
                                                                                                       n_group = .N
                                                                                                   )]
                        }
                    }
                }
            }
        }

        ## Build all reference rows in a list
        ref_rows_list <- lapply(names(xlevels_ref), function(var) {
            ref_level <- xlevels_ref[[var]][1]
            
            ## Get counts for this reference level
            if (var %in% names(ref_counts) && nrow(ref_counts[[var]]) > 0) {
                n_val <- ref_counts[[var]]$n_group[1]
                events_val <- if ("events_group" %in% names(ref_counts[[var]])) {
                                  ref_counts[[var]]$events_group[1]
                              } else {
                                  NA_real_
                              }
            } else {
                n_val <- NA_real_
                events_val <- NA_real_
            }
            
            ## Find first matching row to use as template
            template_idx <- which(dt$variable == var)[1]
            
            if (!is.na(template_idx)) {
                ## Copy the template row
                ref_row <- dt[template_idx, ]
                
                ## Update all necessary columns at once
                ref_row[, `:=`(
                    term = paste0(var, ref_level),
                    group = ref_level,
                    n = ifelse(is.na(n_val), n[1], n_val),
                    events = ifelse(is.na(events_val), events[1], events_val),
                    n_group = n_val,
                    events_group = events_val,
                    coefficient = 0,
                    se = NA_real_,
                    coef = 0,
                    coef_lower = NA_real_,
                    coef_upper = NA_real_,
                    exp_coef = 1,
                    exp_lower = NA_real_,
                    exp_upper = NA_real_,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    reference = reference_label
                )]
                
                ## Update model-specific columns
                if ("HR" %in% names(ref_row)) {
                    ref_row[, `:=`(HR = 1, CI_lower = NA_real_, CI_upper = NA_real_)]
                }
                if ("OR" %in% names(ref_row)) {
                    ref_row[, `:=`(OR = 1, CI_lower = NA_real_, CI_upper = NA_real_)]
                }
                
                return(ref_row)
            }
            return(NULL)
        })
        
        ## Remove NULL entries and bind all reference rows at once
        ref_rows_list <- ref_rows_list[!sapply(ref_rows_list, is.null)]

        if (length(ref_rows_list) > 0) {
            ref_rows_dt <- rbindlist(ref_rows_list, use.names = TRUE, fill = TRUE)
            
            ## Combine main dt and reference rows
            dt <- rbind(dt, ref_rows_dt, use.names = TRUE, fill = TRUE)
            
            ## Sort to put reference rows in correct positions
            ## Within each variable, reference row comes first
            dt[, .is_ref := reference == reference_label]
            setorder(dt, variable, -.is_ref, group)
            dt[, .is_ref := NULL]
        }
    }

    ## Update n and events with group-specific counts where available  
    dt[!is.na(n_group), n := n_group]
    dt[!is.na(events_group), events := events_group]
    
    ## Filter excluded terms
    if (!is.null(terms_to_exclude)) {
        dt <- dt[!term %in% terms_to_exclude]
    }
    
    ## Add significance indicators
    dt[, `:=`(
        sig = data.table::fcase(
                              is.na(p_value), "",
                              p_value < 0.001, "***",
                              p_value < 0.01, "**",
                              p_value < 0.05, "*",
                              p_value < 0.1, ".",
                              default = ""
                          ),
        sig_binary = !is.na(p_value) & p_value < 0.05
    )]
    
    ## Add attributes
    data.table::setattr(dt, "model_class", model_class)
    data.table::setattr(dt, "formula_str", deparse(stats::formula(model)))
    
    if (model_class == "glm") {
        data.table::setattr(dt, "model_family", model$family$family)
        data.table::setattr(dt, "model_link", model$family$link)
    }
    
    dt[]
    return(dt)
}
