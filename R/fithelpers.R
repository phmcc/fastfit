#' Format model results for publication-ready display
#' @keywords internal
format_model_table <- function(data, 
                               effect_col = NULL,
                               digits = 2, 
                               digits_p = 3, 
                               labels = NULL,
                               show_n = TRUE,
                               show_events = TRUE,
                               reference_label = "reference",
                               exponentiate = NULL) {
    
    result <- data.table::copy(data)

    ## Disallow "Events" column if linear model 
    if ("model_type" %in% names(result)) {
        model_type <- unique(result$model_type)[1]
        if (grepl("Linear", model_type, ignore.case = TRUE)) {
            show_events <- FALSE
        }
    }
    
    ## Standardize column names
    if ("variable" %in% names(result)) {
        setnames(result, "variable", "Variable")
    }
    if ("group" %in% names(result)) {
        setnames(result, "group", "Group")
    }

    ## Handle the exponentiate parameter to choose which columns to use
    if (!is.null(exponentiate)) {
        if (exponentiate && "exp_coef" %in% names(result)) {
            ## Check model type
            if ("OR" %in% names(result)) {
                effect_col <- "OR"
            } else if ("HR" %in% names(result)) {
                effect_col <- "HR"
            } else if ("RR" %in% names(result)) {
                effect_col <- "RR"
            } else {
                ## Generic model - create OR/RR based on model type
                model_type <- unique(result$model_type)[1]
                if (grepl("Logistic", model_type)) {
                    result[, `:=`(
                        OR = exp_coef,
                        CI_lower = exp_lower,
                        CI_upper = exp_upper
                    )]
                    effect_col <- "OR"
                } else if (grepl("Poisson", model_type)) {
                    result[, `:=`(
                        RR = exp_coef,
                        CI_lower = exp_lower,
                        CI_upper = exp_upper
                    )]
                    effect_col <- "RR"
                } else {
                    result[, `:=`(
                        Estimate = exp_coef,
                        CI_lower = exp_lower,
                        CI_upper = exp_upper
                    )]
                    effect_col <- "Estimate"
                }
            }
        } else if (!exponentiate && "coef" %in% names(result)) {
            result[, `:=`(
                Estimate = coef,
                CI_lower = coef_lower,
                CI_upper = coef_upper
            )]
            effect_col <- "Estimate"
        }
    }
    
    ## Auto-detect effect column if not specified
    if (is.null(effect_col)) {
        effect_col <- intersect(c("OR", "HR", "RR", "Estimate"), names(result))[1]
        if (length(effect_col) == 0) {
            stop("No effect measure column found (OR, HR, RR, or Estimate)")
        }
    }
    
    ## Apply variable labels if provided - with interaction term support
    if (!is.null(labels) && "Variable" %in% names(result) && length(labels) > 0) {
        
        ## Create lookup table for main effects
        label_dt <- data.table::data.table(
                                    var_orig = names(labels),
                                    var_new = unname(unlist(labels))
                                )
        
        ## Update main effect variable names using merge
        result[label_dt, Variable := i.var_new, on = .(Variable = var_orig)]
        
        ## Handle interaction terms (contain ":")
        interaction_rows <- grep(":", result$Variable, fixed = TRUE)
        
        if (length(interaction_rows) > 0) {
            for (idx in interaction_rows) {
                original_var <- result$Variable[idx]
                
                ## Check if there's a direct custom label for this exact interaction
                if (original_var %in% names(labels)) {
                    result$Variable[idx] <- labels[[original_var]]
                    next
                }
                
                ## Split on ":" to get components
                components <- strsplit(original_var, ":", fixed = TRUE)[[1]]
                
                labeled_parts <- character(length(components))
                
                for (j in seq_along(components)) {
                    comp <- components[j]
                    found_label <- FALSE
                    
                    ## Try to match against base variable names
                    ## Sort by length (longest first) to match more specific names first
                    sorted_vars <- names(labels)[order(-nchar(names(labels)))]
                    
                    for (base_var in sorted_vars) {
                        if (startsWith(comp, base_var)) {
                            ## Get the level/category suffix
                            suffix <- substring(comp, nchar(base_var) + 1)
                            
                            if (nchar(suffix) == 0) {
                                ## Just the variable name (continuous)
                                labeled_parts[j] <- labels[[base_var]]
                            } else {
                                ## Variable + category level
                                labeled_parts[j] <- paste0(suffix, " ", labels[[base_var]])
                            }
                            found_label <- TRUE
                            break
                        }
                    }
                    
                    ## If no match, keep original
                    if (!found_label) {
                        labeled_parts[j] <- comp
                    }
                }
                
                ## Combine with " x "
                result$Variable[idx] <- paste(labeled_parts, collapse = " * ")
            }
        }
    }
    
    ## Clean up Group display (handle empty groups for continuous vars)
    if ("Group" %in% names(result)) {
        result[Group == "", Group := "-"]
    }

    ## IMPORTANT FIX: Format sample size - prioritize n_group where available
    ## This is the key fix for mixed models
    if ("n" %in% names(result)) {
        ## Create a display_n column that uses n_group where available, n otherwise
        if ("n_group" %in% names(result)) {
            ## For rows with n_group data (factor levels), use n_group
            ## For rows without n_group (continuous vars, reference rows), use total n
            result[, display_n := data.table::fifelse(!is.na(n_group), n_group, n)]
        } else {
            ## No n_group column, just use n
            result[, display_n := n]
        }
        
        ## Format the display_n column with commas
        result[, n := data.table::fcase(
                                      is.na(display_n), NA_character_,
                                      display_n >= 1000, format(display_n, big.mark = ","),
                                      default = as.character(display_n)
                                  )]
        
        ## Clean up temporary column
        result[, display_n := NULL]
    }

    ## Similar fix for events
    if ("events" %in% names(result)) {
        ## Create a display_events column that uses events_group where available
        if ("events_group" %in% names(result)) {
            ## For rows with events_group data (factor levels), use events_group
            ## For rows without events_group (continuous vars, reference rows), use total events
            result[, display_events := data.table::fifelse(!is.na(events_group), events_group, events)]
        } else {
            ## No events_group column, just use events
            result[, display_events := events]
        }
        
        ## Format the display_events column with commas
        result[, events := data.table::fcase(
                                           is.na(display_events), NA_character_,
                                           display_events >= 1000, format(display_events, big.mark = ","),
                                           default = as.character(display_events)
                                       )]
        
        ## Clean up temporary column
        result[, display_events := NULL]
    }
    
    ## Eliminate repeated variable names - OPTIMIZED VECTORIZED VERSION
    if ("Variable" %in% names(result) && nrow(result) > 1) {
        ## Create a logical vector indicating where variable changes
        ## TRUE for first occurrence, FALSE for repeats
        n_rows <- nrow(result)
        is_new_var <- c(TRUE, result$Variable[-1] != result$Variable[-n_rows])
        
        ## Set non-first occurrences to empty string
        result[!is_new_var, Variable := ""]
    }
    
    ## Create effect column label based on model scope
    model_scope <- if ("model_scope" %in% names(result)) {
                       unique(result$model_scope)[1]
                   } else {
                       "Effect"
                   }
    
    ## Create appropriate label
    if (model_scope == "Univariable") {
        effect_label <- paste0("Univariable ", effect_col, " (95% CI)")
    } else if (model_scope == "Multivariable") {
        adjusted_col <- if (effect_col == "OR") "aOR" 
                        else if (effect_col == "HR") "aHR"
                        else if (effect_col == "RR") "aRR"
                        else effect_col
        effect_label <- paste0("Multivariable ", adjusted_col, " (95% CI)")
    } else {
        effect_label <- paste0(effect_col, " (95% CI)")
    }
    
    ## Format effect sizes with CI
    if ("CI_lower" %in% names(result) && "CI_upper" %in% names(result)) {
        ## Initialize is_reference as a vector with same length as result
        is_reference <- rep(FALSE, nrow(result))
        if ("reference" %in% names(result)) {
            is_reference <- !is.na(result$reference) & result$reference == reference_label
        }
        
        if (effect_col %in% c("OR", "HR", "RR")) {
            result[, (effect_label) := data.table::fcase(
                                                       is_reference, reference_label,
                                                       !is.na(get(effect_col)), sprintf("%.*f (%.*f-%.*f)",
                                                                                        digits, get(effect_col),
                                                                                        digits, CI_lower,
                                                                                        digits, CI_upper),
                                                       default = ""
                                                   )]
        } else {
            result[, (effect_label) := data.table::fcase(
                                                       is_reference, reference_label,
                                                       !is.na(get(effect_col)), sprintf("%.*f (%.*f, %.*f)",
                                                                                        digits, get(effect_col),
                                                                                        digits, CI_lower,
                                                                                        digits, CI_upper),
                                                       default = ""
                                                   )]
        }
    }
    
    ## Format p-values
    if ("p_value" %in% names(result)) {
        result[, `p-value` := format_pvalues_fit(p_value, digits_p)]
        
        if ("reference" %in% names(result)) {
            result[!is.na(reference) & reference == reference_label, `p-value` := "-"]
        }
    }

    ## Select columns for final output
    display_cols <- character()
    
    if ("Variable" %in% names(result)) display_cols <- c(display_cols, "Variable")
    if ("Group" %in% names(result)) display_cols <- c(display_cols, "Group")

    if ((show_n) && ("n" %in% names(result))) {
        display_cols <- c(display_cols, "n")
    }
    
    if ((show_events) && ("events" %in% names(result))) {
        display_cols <- c(display_cols, "events")
    }
    
    if (effect_label %in% names(result)) display_cols <- c(display_cols, effect_label)
    if ("p-value" %in% names(result)) display_cols <- c(display_cols, "p-value")
    
    formatted <- result[, ..display_cols]

    if ("events" %in% names(formatted)) {
        setnames(formatted, "events", "Events")
    }

    return(formatted)
}

#' Format p-values for display - OPTIMIZED
#' @keywords internal
format_pvalues_fit <- function(p, digits = 3) {
    data.table::fcase(
                    is.na(p), "-",
                    p < 0.001, "< 0.001",
                    default = sprintf(paste0("%.", digits, "f"), p)
                )
}
