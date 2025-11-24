#' Get factor levels from model (works with S3 and S4)
#' @keywords internal
get_model_xlevels <- function(model) {
    
    ## For regular S3 models (glm, lm, coxph, clogit)
    if (!isS4(model)) {
        ## Special handling for coxme - it doesn't store xlevels
        if (inherits(model, "coxme")) {
            return(NULL)  # Will be handled separately
        }
        return(model$xlevels)
    }
    
    ## For S4 lme4 models (glmer, lmer)
    if (inherits(model, "merMod")) {
        ## Extract from the model frame
        frame_data <- model@frame
        
        xlevels <- list()
        for (var_name in names(frame_data)) {
            if (is.factor(frame_data[[var_name]])) {
                xlevels[[var_name]] <- levels(frame_data[[var_name]])
            }
        }
        
        return(if (length(xlevels) > 0) xlevels else NULL)
    }
    
    ## Default: return NULL
    return(NULL)
}

#' Get data from model object (works with S3 and S4)
#' @keywords internal
get_model_data <- function(model) {
    
    ## Try attribute first
    data_attr <- attr(model, "data")
    if (!is.null(data_attr)) {
        return(data_attr)
    }
    
    ## For S3 objects, try model$data
    if (!isS4(model)) {
        if (!is.null(model$data)) {
            return(model$data)
        }
        if (!is.null(model$model)) {
            return(model$model)
        }
    }
    
    ## For S4 lme4 objects, use model frame
    if (inherits(model, "merMod")) {
        return(model@frame)
    }
    
    ## Default
    return(NULL)
}

#' Detect if model is univariable or multivariable
#' @keywords internal
detect_model_type <- function(model) {
    ## For glmer/glmerMod models, handle random effects specially
    if (inherits(model, c("glmerMod", "merMod"))) {
                                        # Get the fixed effects formula (excluding random effects)
        if (!is.null(model@call$formula)) {
            formula_str <- as.character(model@call$formula)[3]  # RHS of formula
                                        # Remove random effects terms (anything with |)
            formula_str <- gsub("\\([^|]*\\|[^)]*\\)", "", formula_str)
                                        # Clean up extra spaces and plus signs
            formula_str <- gsub("\\s*\\+\\s*\\+", " +", formula_str)
            formula_str <- gsub("^\\s*\\+\\s*", "", formula_str)
            formula_str <- gsub("\\s*\\+\\s*$", "", formula_str)
            
                                        # Split by + to get terms
            terms <- trimws(strsplit(formula_str, "\\+")[[1]])
            terms <- terms[terms != ""]
            
                                        # Count unique base variables (handle interactions)
            unique_vars <- character()
            for (term in terms) {
                if (grepl(":", term)) {
                                        # Interaction term - split and add components
                    components <- trimws(strsplit(term, ":")[[1]])
                    unique_vars <- c(unique_vars, components)
                } else {
                    unique_vars <- c(unique_vars, term)
                }
            }
            unique_vars <- unique(unique_vars)
            
            return(data.table::fifelse(length(unique_vars) == 1, "Univariable", "Multivariable"))
        }
    }
    
    ## Original code for other model types
    ## Get coefficient names once
    coef_names <- names(stats::coef(model))
    
    ## Remove intercept from count
    if ("(Intercept)" %in% coef_names) {
        term_names <- coef_names[coef_names != "(Intercept)"]
    } else {
        term_names <- coef_names
    }
    
    ## Identify interaction terms (contain ":")
    is_interaction <- grepl(":", term_names, fixed = TRUE)
    interaction_terms <- term_names[is_interaction]
    main_terms <- term_names[!is_interaction]
    
    ## For models without factor variables, count unique base variables from main effects
    xlevels <- get_model_xlevels(model)
    
    ## Special handling for coxme - reconstruct factor info
    if (inherits(model, "coxme") && is.null(xlevels)) {
        data_source <- get_model_data(model)
        if (!is.null(data_source)) {
            xlevels <- list()
            formula_vars <- all.vars(stats::formula(model))
            if (length(formula_vars) >= 2) {
                predictor_vars <- formula_vars[-(1:2)]  # Remove Surv components
                predictor_vars <- predictor_vars[!grepl("\\|", predictor_vars)]  # Remove random effects
                
                for (var in predictor_vars) {
                    if (var %in% names(data_source) && is.factor(data_source[[var]])) {
                        xlevels[[var]] <- levels(data_source[[var]])
                    }
                }
            }
            if (length(xlevels) == 0) xlevels <- NULL
        }
    }
    
    if (is.null(xlevels) || length(xlevels) == 0) {
        ## Count main effect terms
        n_terms <- length(main_terms)
        
        ## If we have interactions, extract unique base variables from them too
        if (length(interaction_terms) > 0) {
            ## Extract variables from interactions (split on ":")
            interaction_vars <- unique(unlist(strsplit(interaction_terms, ":", fixed = TRUE)))
            ## Count any interaction variables not already in main effects
            new_vars <- setdiff(interaction_vars, main_terms)
            n_terms <- n_terms + length(new_vars)
        }
        
        return(data.table::fifelse(n_terms == 1, "Univariable", "Multivariable"))
    }
    
    ## For models with factor variables, count unique base variables
    factor_pattern <- paste0("^(", paste(names(xlevels), collapse = "|"), ")")
    is_factor_term <- grepl(factor_pattern, main_terms)
    factor_vars_present <- unique(sapply(main_terms[is_factor_term], function(term) {
        for (var in names(xlevels)) {
            if (startsWith(term, var)) return(var)
        }
        return(NA_character_)
    }))
    factor_vars_present <- factor_vars_present[!is.na(factor_vars_present)]
    n_vars <- length(factor_vars_present)
    
    ## Count unique continuous variables IN MAIN EFFECTS
    continuous_terms <- main_terms[!is_factor_term]
    n_continuous <- length(unique(continuous_terms))
    
    ## Total unique variables (from main effects only)
    ## Interactions don't add new variables to the count
    n_terms <- n_vars + n_continuous
    
    return(data.table::fifelse(n_terms == 1, "Univariable", "Multivariable"))
}

#' Get readable model type name
#' @keywords internal
get_model_type_name <- function(model) {
    model_class <- class(model)[1]
    
    ## Remove wrapper classes
    if (model_class == "mmodel") {
        model_class <- class(model)[2]
    }
    
    ## For GLM, be more specific based on family
    if (model_class == "glm") {
        family <- model$family$family
        
        ## OPTIMIZED: Use switch instead of multiple if-else
        return(switch(family,
                      binomial = "Logistic",
                      poisson = "Poisson",
                      gaussian = "Linear (GLM)",
                      Gamma = "Gamma",
                      quasibinomial = "Quasi-Binomial",
                      quasipoisson = "Quasi-Poisson",
                      paste0(stringr::str_to_title(family), " GLM")
                      ))
    }
    
    ## Map to readable names for other model types
    ## OPTIMIZED: Use switch instead of named vector lookup
    return(switch(model_class,
                  lm = "Linear",
                  coxph = "Cox PH",
                  clogit = "Conditional Logistic",
                  coxme = "Mixed Effects Cox",
                  glmer = "glmerMod",  # Keep full name for clarity
                  glmerMod = "glmerMod",
                  lmer = "Mixed Effects Linear",
                  model_class  # default
                  ))
}

#' Parse term into variable and group
#' @keywords internal
parse_term <- function(terms, xlevels = NULL, model = NULL) {
    n_terms <- length(terms)
    
    ## Initialize result vectors
    variable <- character(n_terms)
    group <- character(n_terms)
    
    ## Interactions should not be parsed as factor variables
    is_interaction <- grepl(":", terms, fixed = TRUE)

    ## Special handling for coxme models - reconstruct xlevels if needed
    if (!is.null(model) && inherits(model, "coxme") && is.null(xlevels)) {
        data_source <- attr(model, "data")
        if (is.null(data_source)) {
            data_source <- get_model_data(model)
        }
        
        if (!is.null(data_source)) {
            xlevels <- list()
            
            ## Extract factor structure from coefficient names
            coef_names <- names(coxme::fixef(model))
            
            ## Check each column in the data
            for (col_name in names(data_source)) {
                ## Check if any coefficient starts with this column name
                if (any(grepl(paste0("^", col_name), coef_names))) {
                    ## This is a factor variable - get its levels
                    if (is.factor(data_source[[col_name]])) {
                        xlevels[[col_name]] <- levels(data_source[[col_name]])
                    } else if (is.character(data_source[[col_name]])) {
                        xlevels[[col_name]] <- sort(unique(data_source[[col_name]]))
                    }
                }
            }
            
            if (length(xlevels) == 0) xlevels <- NULL
        }
    }
    
    if (!is.null(xlevels) && length(xlevels) > 0) {
        ## OPTIMIZED: Vectorized approach
        xlevel_names <- names(xlevels)
        
        ## For each factor variable, find all matching terms at once
        ## BUT SKIP INTERACTION TERMS
        for (var in xlevel_names) {
            ## Find which terms start with this variable name
            pattern <- paste0("^", var)
            matches <- grepl(pattern, terms) & !is_interaction  # Skip interactions
            
            if (any(matches)) {
                ## Extract levels for all matching terms at once
                variable[matches] <- var
                group[matches] <- sub(pattern, "", terms[matches])
            }
        }
        
        ## Any remaining unmatched terms (including interactions) are treated as-is
        unmatched <- variable == ""
        if (any(unmatched)) {
            variable[unmatched] <- terms[unmatched]
            ## group already initialized to "" for these
        }
    } else {
        ## No factor variables - all continuous (including interactions)
        variable <- terms
        ## group already initialized to ""
    }
    
    return(data.table::data.table(variable = variable, group = group))
}

#' Extract event variable from survival model
#' @keywords internal
get_event_variable <- function(model, model_class) {
    event_var <- NULL
    
    if (model_class %in% c("coxph", "clogit", "coxme")) {
        
        ## Get the formula string differently for each model type
        outcome_str <- NULL
        
        if (model_class == "coxme") {
            ## For coxme, use formulaList$fixed
            if (!is.null(model$formulaList$fixed)) {
                outcome_str <- tryCatch({
                    as.character(model$formulaList$fixed)[2]
                }, error = function(e) NULL)
            }
        } else {
            ## For coxph and clogit, use standard formula
            outcome_str <- tryCatch({
                as.character(stats::formula(model))[2]
            }, error = function(e) NULL)
        }
        
        ## Extract event variable from Surv()
        if (!is.null(outcome_str) && !is.na(outcome_str)) {
            if (grepl("Surv\\(", outcome_str)) {
                ## Remove "Surv(" from beginning and ")" from end
                surv_expr <- gsub("^Surv\\(", "", outcome_str)
                surv_expr <- gsub("\\)$", "", surv_expr)
                
                ## Split by comma to get time and event
                surv_parts <- trimws(strsplit(surv_expr, ",")[[1]])
                if (length(surv_parts) >= 2) {
                    event_var <- surv_parts[2]
                }
            }
        }
    }
    
    return(event_var)
}
