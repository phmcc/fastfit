#' Create Forest Plot for Cox Proportional Hazards Models
#'
#' Generates a publication-ready forest plot that combines a formatted data table 
#' with a graphical representation of hazard ratios from a Cox proportional hazards 
#' survival model. The plot integrates variable names, group levels, sample sizes, 
#' event counts, hazard ratios with confidence intervals, p-values, and model 
#' diagnostics in a single comprehensive visualization designed for manuscripts 
#' and presentations.
#'
#' @param model A fitted Cox proportional hazards model object of class \code{coxph} 
#'   from the \code{survival} package. Should be fitted using 
#'   \code{\link[survival]{coxph}}.
#'   
#' @param data A data.frame or data.table containing the original data used to 
#'   fit the model. If \code{NULL} (default), the function attempts to extract 
#'   data from the model object (\code{model$data}, evaluation of \code{model$call$data}, 
#'   or \code{model$model}). Providing data explicitly is strongly recommended.
#'   
#' @param title Character string specifying the plot title displayed at the top. 
#'   Default is \code{"Cox Proportional Hazards Model"}. Use descriptive titles 
#'   for publication.
#'   
#' @param effect_label Character string for the effect measure label on the 
#'   forest plot axis. Default is \code{"Hazard Ratio"}. Alternatives might include 
#'   "HR" for space-constrained plots.
#'   
#' @param digits Integer specifying the number of decimal places for hazard ratios 
#'   and confidence intervals. Default is 2.
#'   
#' @param font_size Numeric multiplier controlling the base font size for all 
#'   text elements. Default is 1.0.
#'   
#' @param annot_size Numeric value controlling the relative font size for data 
#'   annotations. Default is 3.88.
#'   
#' @param header_size Numeric value controlling the relative font size for column 
#'   headers. Default is 5.82.
#'   
#' @param title_size Numeric value controlling the relative font size for the 
#'   main plot title. Default is 23.28.
#'   
#' @param tbl_width Numeric value between 0 and 1 specifying the proportion of 
#'   total plot width allocated to the data table. Default is 0.6 (60\% table, 
#'   40\% forest plot).
#'   
#' @param plot_width Numeric value specifying the intended output width in 
#'   specified \code{units}. Used for optimizing layout. Default is \code{NULL} 
#'   (automatic). Recommended: 10-16 inches.
#'   
#' @param plot_height Numeric value specifying the intended output height in 
#'   specified \code{units}. Default is \code{NULL} (automatic based on number 
#'   of rows).
#'   
#' @param show_n Logical. If \code{TRUE}, includes a column showing group-specific 
#'   sample sizes. Default is \code{TRUE}.
#'   
#' @param show_events Logical. If \code{TRUE}, includes a column showing the 
#'   number of events (deaths, failures) for each group. Critical for survival 
#'   analysis interpretation. Default is \code{TRUE}.
#'   
#' @param indent_groups Logical. If \code{TRUE}, indents factor levels under 
#'   their parent variable, creating hierarchical structure. The "Group" column 
#'   is hidden when \code{TRUE}. Default is \code{FALSE}.
#'   
#' @param condense_table Logical. If \code{TRUE}, condenses binary categorical 
#'   variables into single rows. Automatically sets \code{indent_groups = TRUE}. 
#'   Default is \code{FALSE}.
#'   
#' @param center_padding Numeric value specifying horizontal spacing between 
#'   table and forest plot. Default is 4.
#'   
#' @param zebra_stripes Logical. If \code{TRUE}, applies alternating gray 
#'   background shading to different variables for improved readability. 
#'   Default is \code{TRUE}.
#'   
#' @param ref_label Character string to display for reference categories. 
#'   Default is \code{"reference"}.
#'   
#' @param labels Named character vector providing custom display labels for 
#'   variables. Example: \code{c(age = "Age (years)", stage = "Disease Stage")}. 
#'   Default is \code{NULL}.
#'   
#' @param color Character string specifying the color for hazard ratio point 
#'   estimates in the forest plot. Default is \code{"#8A61D8"} (purple). Use 
#'   hex codes or R color names.
#'   
#' @param units Character string specifying units for plot dimensions: 
#'   \code{"in"} (inches), \code{"cm"}, or \code{"mm"}. Default is \code{"in"}.
#'
#' @return A \code{ggplot} object containing the complete forest plot. The plot 
#'   can be displayed, saved, or further customized.
#'   
#'   The returned object includes an attribute \code{"recommended_dims"} 
#'   accessible via \code{attr(plot, "recommended_dims")}, containing:
#'   \describe{
#'     \item{width}{Numeric. Recommended plot width in specified units}
#'     \item{height}{Numeric. Recommended plot height in specified units}
#'   }
#'   
#'   Recommendations are printed to console if dimensions are not specified.
#'
#' @details
#' \strong{Survival-Specific Features:}
#' 
#' The Cox forest plot includes several survival analysis-specific components:
#' \itemize{
#'   \item \strong{Event counts}: Number of events (deaths, failures) shown for 
#'     each predictor category, critical for assessing statistical power
#'   \item \strong{Hazard ratios}: Always exponentiated coefficients (never raw), 
#'     interpreted as the multiplicative change in hazard
#'   \item \strong{Log scale}: Forest plot uses log scale for HR (reference line at 1)
#'   \item \strong{Model diagnostics}: Includes concordance (C-index), global 
#'     log-rank test p-value, and AIC
#' }
#' 
#' \strong{Plot Components:}
#' 
#' \enumerate{
#'   \item \strong{Title}: Centered at top
#'   \item \strong{Data Table} (left): Contains:
#'     \itemize{
#'       \item Variable and Group columns
#'       \item n: Sample sizes by group
#'       \item Events: Event counts by group (critical for survival)
#'       \item aHR (95\% CI); \emph{p}-value: Adjusted hazard ratios with CIs and p-values
#'     }
#'   \item \strong{Forest Plot} (right):
#'     \itemize{
#'       \item Point estimates (squares sized by sample size)
#'       \item 95\% confidence intervals
#'       \item Reference line at HR = 1
#'       \item Log scale for hazard ratios
#'     }
#'   \item \strong{Model Statistics} (footer):
#'     \itemize{
#'       \item Events analyzed (with percentage of total)
#'       \item Global log-rank test p-value
#'       \item Concordance (C-index) with standard error
#'       \item AIC
#'     }
#' }
#' 
#' \strong{Interpreting Hazard Ratios:}
#' 
#' \itemize{
#'   \item \strong{HR = 1}: No effect on hazard (reference)
#'   \item \strong{HR > 1}: Increased hazard (worse survival)
#'   \item \strong{HR < 1}: Decreased hazard (better survival)
#'   \item Example: HR = 2.0 means twice the hazard of the event at any time
#' }
#' 
#' \strong{Event Counts:}
#' 
#' The "Events" column is particularly important in survival analysis:
#' \itemize{
#'   \item Indicates the number of actual events (not censored observations) 
#'     in each group
#'   \item Essential for assessing statistical power
#'   \item Categories with very few events may have unreliable HR estimates
#'   \item The footer shows total events analyzed and percentage of all events 
#'     in the original data
#' }
#' 
#' \strong{Concordance (C-index):}
#' 
#' The concordance statistic displayed in the footer indicates discrimination:
#' \itemize{
#'   \item Range: 0.5 to 1.0
#'   \item 0.5 = random prediction (coin flip)
#'   \item 0.7-0.8 = acceptable discrimination
#'   \item >0.8 = excellent discrimination
#'   \item Standard error provided for confidence interval calculation
#' }
#' 
#' \strong{Global Log-Rank Test:}
#' 
#' The global p-value tests the null hypothesis that all coefficients are zero:
#' \itemize{
#'   \item Significant p-value (< 0.05) indicates the model as a whole predicts survival
#'   \item Non-significant global test doesn't preclude significant individual 
#'     predictors
#'   \item Based on the score (log-rank) test
#' }
#' 
#' \strong{Stratification and Clustering:}
#' 
#' If the model includes stratification (\code{strata()}) or clustering 
#' (\code{cluster()}):
#' \itemize{
#'   \item Stratified variables are not shown in the forest plot (they don't 
#'     have HRs)
#'   \item Clustering affects standard errors but not point estimates
#'   \item Both are handled automatically by the function
#' }
#' 
#' \strong{Proportional Hazards Assumption:}
#' 
#' The forest plot assumes proportional hazards (constant HR over time). Users 
#' should verify this assumption using:
#' \itemize{
#'   \item \code{cox.zph(model)} for testing
#'   \item Stratification for variables violating the assumption
#'   \item Time-dependent coefficients if needed
#' }
#'
#' @seealso 
#' \code{\link{glmforest}} for logistic/GLM forest plots,
#' \code{\link{lmforest}} for linear model forest plots,
#' \code{\link[survival]{coxph}} for fitting Cox models,
#' \code{\link[survival]{cox.zph}} for testing proportional hazards assumption,
#' \code{\link{fit}} for fastfit regression modeling
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' library(survival)
#' 
#' # Example 1: Basic Cox model forest plot
#' model1 <- coxph(Surv(os_months, os_status) ~ age + sex + treatment,
#'                 data = clintrial)
#' 
#' plot1 <- coxforest(model1, data = clintrial)
#' print(plot1)
#' 
#' # Example 2: With custom labels and title
#' plot2 <- coxforest(
#'     model = model1,
#'     data = clintrial,
#'     title = "Prognostic Factors for Overall Survival",
#'     labels = clintrial_labels
#' )
#' print(plot2)
#' 
#' # Example 3: Comprehensive multivariable model
#' model3 <- coxph(
#'     Surv(os_months, os_status) ~ age + sex + bmi + smoking + 
#'         treatment + stage + grade,
#'     data = clintrial
#' )
#' 
#' plot3 <- coxforest(
#'     model = model3,
#'     data = clintrial,
#'     labels = clintrial_labels,
#'     indent_groups = TRUE
#' )
#' print(plot3)
#' 
#' # Example 4: Condensed layout for many binary predictors
#' model4 <- coxph(
#'     Surv(os_months, os_status) ~ age + sex + smoking + 
#'         hypertension + diabetes + surgery,
#'     data = clintrial
#' )
#' 
#' plot4 <- coxforest(
#'     model = model4,
#'     data = clintrial,
#'     condense_table = TRUE,
#'     labels = clintrial_labels
#' )
#' print(plot4)
#' 
#' # Example 5: Custom color scheme
#' plot5 <- coxforest(
#'     model = model1,
#'     data = clintrial,
#'     color = "#E74C3C",  # Red
#'     zebra_stripes = FALSE,
#'     labels = clintrial_labels
#' )
#' print(plot5)
#' 
#' # Example 6: Hide sample sizes, show only events
#' plot6 <- coxforest(
#'     model = model1,
#'     data = clintrial,
#'     show_n = FALSE,
#'     show_events = TRUE,
#'     labels = clintrial_labels
#' )
#' print(plot6)
#' 
#' # Example 7: Adjust table width
#' plot7 <- coxforest(
#'     model = model3,
#'     data = clintrial,
#'     tbl_width = 0.65,  # More space for long variable names
#'     labels = clintrial_labels
#' )
#' print(plot7)
#' 
#' # Example 8: Specify exact dimensions
#' plot8 <- coxforest(
#'     model = model1,
#'     data = clintrial,
#'     plot_width = 14,
#'     plot_height = 8,
#'     labels = clintrial_labels
#' )
#' 
#' # Example 9: Use recommended dimensions for saving
#' plot9 <- coxforest(model1, data = clintrial)
#' dims <- attr(plot9, "recommended_dims")
#' 
#' # Save to PDF
#' # ggsave("survival_forest.pdf", plot9,
#' #        width = dims$width, height = dims$height)
#' 
#' # Example 10: Different units (centimeters)
#' plot10 <- coxforest(
#'     model = model1,
#'     data = clintrial,
#'     plot_width = 35,
#'     plot_height = 25,
#'     units = "cm",
#'     labels = clintrial_labels
#' )
#' 
#' # Example 11: Stratified Cox model
#' # Stratification variable (site) won't appear in plot
#' model11 <- coxph(
#'     Surv(os_months, os_status) ~ age + sex + treatment + strata(site),
#'     data = clintrial
#' )
#' 
#' plot11 <- coxforest(
#'     model = model11,
#'     data = clintrial,
#'     title = "Stratified by Study Site",
#'     labels = clintrial_labels
#' )
#' print(plot11)
#' 
#' # Example 12: With clustering for robust SE
#' model12 <- coxph(
#'     Surv(os_months, os_status) ~ age + sex + treatment + cluster(site),
#'     data = clintrial
#' )
#' 
#' plot12 <- coxforest(
#'     model = model12,
#'     data = clintrial,
#'     title = "Clustered by Study Site",
#'     labels = clintrial_labels
#' )
#' print(plot12)
#' # Standard errors account for site clustering
#' 
#' # Example 13: Custom reference label
#' plot13 <- coxforest(
#'     model = model1,
#'     data = clintrial,
#'     ref_label = "1.00 (ref)",
#'     labels = clintrial_labels
#' )
#' print(plot13)
#' 
#' # Example 14: Presentation-sized fonts
#' plot14 <- coxforest(
#'     model = model1,
#'     data = clintrial,
#'     font_size = 1.4,
#'     title_size = 28,
#'     labels = clintrial_labels
#' )
#' print(plot14)
#' 
#' # Example 15: Publication-ready final plot
#' final_model <- coxph(
#'     Surv(os_months, os_status) ~ age + sex + bmi + smoking + 
#'         hypertension + diabetes + ecog + treatment + stage + grade,
#'     data = clintrial
#' )
#' 
#' final_plot <- coxforest(
#'     model = final_model,
#'     data = clintrial,
#'     title = "Multivariable Cox Regression: Prognostic Factors for Overall Survival",
#'     labels = clintrial_labels,
#'     indent_groups = TRUE,
#'     zebra_stripes = TRUE,
#'     show_n = TRUE,
#'     show_events = TRUE,
#'     color = "#8A61D8"
#' )
#' 
#' # Check model fit
#' summary(final_model)
#' 
#' # Verify proportional hazards assumption
#' # cox.zph(final_model)
#' 
#' # Save for publication
#' dims <- attr(final_plot, "recommended_dims")
#' # ggsave("figure2_survival.pdf", final_plot,
#' #        width = dims$width, height = dims$height)
#' # ggsave("figure2_survival.png", final_plot,
#' #        width = dims$width, height = dims$height, dpi = 300)
#'
#' @export
coxforest <- function(model, data = NULL,
                      title = "Cox Proportional Hazards Model",
                      effect_label = "Hazard Ratio",
                      digits = 2,
                      font_size = 1.0,
                      annot_size = 3.88,
                      header_size = 5.82,
                      title_size = 23.28,
                      plot_width = NULL,
                      plot_height = NULL,
                      tbl_width = 0.6,
                      show_n = TRUE,
                      show_events = TRUE,
                      indent_groups = FALSE,
                      condense_table = FALSE,
                      center_padding = 4,
                      zebra_stripes = TRUE,
                      ref_label = "reference",
                      labels = NULL,
                      color = "#8A61D8",
                      units = "in") {

    ## Check for required packages
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required but not installed.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required but not installed.")
    }
    if (!requireNamespace("grid", quietly = TRUE)) {
        stop("Package 'grid' is required but not installed.")
    }
    
    stopifnot(inherits(model, c("coxph", "coxme")))

    ## Internally work in inches
    if (!is.null(plot_width) && units != "in") {
        plot_width <- convert_units(plot_width, from = units, to = "in")
    }

    ## Get data
    if(is.null(data)){
        warning("The `data` argument is not provided. Data will be extracted from model fit.")
        
        ## First try to get data from model$data (stored by mmodel)
        if (!is.null(model$data)) {
            data <- model$data
        } else {
            ## Try to evaluate the call
            tryCatch({
                data <- eval(model$call$data)
            }, error = function(e) {
                data <- NULL
            })
            
            ## Try model$model as fallback
            if (is.null(data))
                data <- model$model
        }
        
        if (is.null(data))
            stop("The `data` argument should be provided either to coxforest or stored in the model.")
    }
    
    ## Convert to data.table if not already
    if(!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
    terms <- attr(model$terms, "dataClasses")[-1]

    ## Filter out strata() and cluster() terms - these are not predictors to plot
    terms <- terms[!grepl("strata|cluster", names(terms))]

    ## Extract coefficients based on model type
    model_class <- class(model)[1]

    ## For coxme models, data is mandatory
    if (model_class == "coxme" && is.null(data)) {
        stop("The 'data' argument is required for coxme models.\n",
             "Usage: coxforest(model, data = your_data, ...)")
    }

    ## Initialize coef variables
    coef_summary <- NULL
    conf_int <- NULL
    
    if (model_class == "coxme") {
        ## Special handling for coxme models
        
        ## Extract fixed effects coefficients
        fixed_effects <- coxme::fixef(model)
        
        ## Extract standard errors from variance matrix
        n_fixed <- length(fixed_effects)
        se_values <- sqrt(diag(as.matrix(model$variance))[1:n_fixed])
        
        ## Calculate z-statistics and p-values
        z_values <- fixed_effects / se_values
        p_values <- 2 * pnorm(abs(z_values), lower.tail = FALSE)
        
        ## Calculate 95% confidence intervals
        z_mult <- qnorm(0.975)
        conf_low <- fixed_effects - z_mult * se_values
        conf_high <- fixed_effects + z_mult * se_values
        
        ## Create coefficient summary in expected format
        coef_summary <- cbind(
            coef = fixed_effects,
            "se(coef)" = se_values,
            z = z_values,
            "Pr(>|z|)" = p_values
        )
        rownames(coef_summary) <- names(fixed_effects)
        
        ## Create conf_int in expected format
        conf_int <- cbind(conf_low, conf_high)
        colnames(conf_int) <- c("2.5 %", "97.5 %")
        rownames(conf_int) <- names(fixed_effects)
        
        ## Also need to set nevent and create xlevels as before...
        model$nevent <- model$n[1]

        ## For coxme, systematically extract from formula
        ## Parse formula to get predictor variables
        formula_chr <- as.character(model$formulaList$fixed)
        rhs_formula <- formula_chr[3]

        ## Split by + to get individual terms
        predictor_terms <- trimws(strsplit(rhs_formula, "\\+")[[1]])

        ## Remove random effects terms like (1|site)
        predictor_vars <- predictor_terms[!grepl("\\(.*\\|.*\\)", predictor_terms)]

        ## Create terms
        terms <- sapply(predictor_vars, function(v) {
            if (v %in% names(data)) {
                if (is.factor(data[[v]])) "factor"
                else if (is.numeric(data[[v]])) "numeric"
                else class(data[[v]])[1]
            } else {
                "numeric"  # Default
            }
        })

        ## Create xlevels (keep your existing code)
        if (is.null(model$xlevels)) {
            model$xlevels <- list()
            for (var in predictor_vars) {
                if (var %in% names(data) && is.factor(data[[var]])) {
                    model$xlevels[[var]] <- levels(data[[var]])
                }
            }
        }
        
    } else {
        ## Standard coxph extraction
        coef_summary <- summary(model)$coefficients
        conf_int <- stats::confint(model)
    }    

    
    coef <- data.table::data.table(
                            term = rownames(coef_summary),
                            estimate = coef_summary[, "coef"],
                            std_error = coef_summary[, "se(coef)"],
                            statistic = coef_summary[, "z"],
                            p_value = coef_summary[, "Pr(>|z|)"],
                            conf_low = conf_int[, 1],
                            conf_high = conf_int[, 2]
                        )
    
    ## Get model statistics
    if (model_class == "coxme") {
        ## Calculate AIC for coxme
        ## AIC = -2 * loglik + 2 * df
        coxme_aic <- if (!is.null(model$loglik) && length(model$loglik) >= 3) {
                         -2 * model$loglik[3] + 2 * model$df[2]  # Using penalized loglik
                  } else {
                         NA
                     }
        
        ## Calculate likelihood ratio test p-value
        p_value_log <- if (!is.null(model$loglik) && length(model$loglik) >= 2) {
                           lr_stat <- 2 * (model$loglik[2] - model$loglik[1])
                           lr_df <- model$df[1]
                           pchisq(lr_stat, df = lr_df, lower.tail = FALSE)
                       } else {
                           NA
                       }
        
        gmodel <- list(
            nevent = model$n[1],
            p_value_log = p_value_log,
            AIC = coxme_aic,
            concordance = NA,  # Not easily available for coxme
            concordance.se = NA
        )
    } else {
        ## Standard coxph statistics
        conc_values <- summary(model)$concordance
        gmodel <- list(
            nevent = model$nevent,
            p_value_log = as.numeric(summary(model)$sctest["pvalue"]),
            AIC = stats::AIC(model),
            concordance = as.numeric(conc_values["C"]),
            concordance.se = as.numeric(conc_values["se(C)"])
        )
    }

    ## Calculate total events in original data and percentage analyzed
    formula_terms <- all.vars(model$formula)
    if (length(formula_terms) >= 2) {
        event_var <- formula_terms[2]
        event_data <- data[[event_var]]
        
        ## Handle factor event indicators (rare but possible)
        if (is.factor(event_data)) {
            event_binary <- as.numeric(event_data) == 2
        } else {
            event_binary <- event_data
        }
        
        total_events <- sum(event_binary, na.rm = TRUE)
        gmodel$pct_events_analyzed <- (gmodel$nevent / total_events) * 100
    } else {
        total_events <- NA
        gmodel$pct_events_analyzed <- NA
    }
    
    ## Format events and AIC with commas
    gmodel$nevent_formatted <- format(gmodel$nevent, big.mark = ",", scientific = FALSE)
    gmodel$nevent_with_pct <- if(!is.na(gmodel$pct_events_analyzed)) {
                                  paste0(gmodel$nevent_formatted, " (", 
                                         sprintf("%.1f%%", gmodel$pct_events_analyzed), ")")
                              } else {
                                  gmodel$nevent_formatted
                              }
    gmodel$AIC_formatted <- format(round(gmodel$AIC, 2), big.mark = ",", scientific = FALSE, nsmall = 2)

    ## Calculate 95% CI for concordance if both concordance and SE are available
    if(!is.null(gmodel$concordance) && !is.na(gmodel$concordance) &&
       !is.null(gmodel$concordance.se) && !is.na(gmodel$concordance.se)) {
        gmodel$concordance.lower <- gmodel$concordance - 1.96 * gmodel$concordance.se
        gmodel$concordance.upper <- gmodel$concordance + 1.96 * gmodel$concordance.se
    } else {
        gmodel$concordance.lower <- NA
        gmodel$concordance.upper <- NA
    }

    ## Format the global p-value for display
    global_p_formatted <- if(as.numeric(gmodel$p_value_log) < 0.001) {
                              "< 0.001"
                          } else {
                              format.pval(gmodel$p_value_log, eps = ".001")
                          }
    
    ## Build concordance string based on whether CI is available
    concordance_string <- if(!is.null(gmodel$concordance) && !is.na(gmodel$concordance)) {
                              if(!is.null(gmodel$concordance.lower) && !is.na(gmodel$concordance.lower) && 
                                 !is.null(gmodel$concordance.upper) && !is.na(gmodel$concordance.upper)) {
                                  paste0("Concordance Index: ", round(gmodel$concordance, 2), 
                                         " (95% CI ", round(gmodel$concordance.lower, 2), "-", 
                                         round(gmodel$concordance.upper, 2), ")")
                              } else {
                                  paste0("Concordance Index: ", round(gmodel$concordance, 2))
                              }
                          } else {
                              "Concordance Index: Not available"
                          }

    ## Extract statistics for every variable - preserving order
    all_terms <- lapply(seq_along(terms), function(i){
        var <- names(terms)[i]

        if (terms[i] %in% c("factor", "character")) {
            ## Get the factor levels from the model's xlevels (proper order)
            if(var %in% names(model$xlevels)) {
                factor_levels <- model$xlevels[[var]]
                
                ## Create data table with proper levels
                level_counts <- data[!is.na(get(var)), .N, by = var]
                data.table::setnames(level_counts, c("level", "Freq"))
                
                ## Ensure all levels are present
                all_levels_dt <- data.table::data.table(
                                                 level = factor_levels,
                                                 Freq = 0
                                             )
                
                ## Update with actual counts
                all_levels_dt[level_counts, Freq := i.Freq, on = "level"]
                all_levels_dt[, var := var]
                all_levels_dt[, pos := .I]
                all_levels_dt[, var_order := i]
                all_levels_dt[, .(var, level, Freq, pos, var_order)]
            } else {
                ## Fallback for variables not in xlevels
                adf <- data[!is.na(get(var)), .N, by = var]
                data.table::setnames(adf, old = c(var, "N"), new = c("level", "Freq"))
                adf[, var := var]
                adf[, pos := .I]
                adf[, var_order := i]
                adf[, .(var, level, Freq, pos, var_order)]
            }
        }
        else if (terms[i] == "numeric") {
            data.table::data.table(var = var, level = "-", Freq = nrow(data), pos = 1, var_order = i)
        }
        else {
            vars = grep(paste0("^", var, "*."), coef$term, value=TRUE)
            data.table::data.table(var = vars, level = "-", Freq = nrow(data),
                                   pos = seq_along(vars), var_order = i)
        }
    })
    
    all_terms_df <- data.table::rbindlist(all_terms)

    data.table::setnames(all_terms_df, c("var", "level", "N", "pos", "var_order"))

    ## Add events
    if (model_class == "coxme") {
        formula_terms <- all.vars(model$formulaList$fixed)
    } else {
        formula_terms <- all.vars(model$formula)
    }
    
    if (length(formula_terms) >= 2) {
        ## Assuming standard Surv(time, status) format
        event_var <- formula_terms[2]
        
        ## Calculate events for each term
        all_terms_df[, Events := {
            event_data <- data[[event_var]]
            
            ## Handle factor event indicators (rare but possible)
            if (is.factor(event_data)) {
                ## Convert to binary (assuming second level is the event)
                event_binary <- as.numeric(event_data) == 2
            } else {
                ## Already numeric
                event_binary <- event_data
            }
            
            if (level == "-") {
                ## Continuous variable - total events
                sum(event_binary, na.rm = TRUE)
            } else {
                ## Factor level - events in that group
                sum(event_binary[data[[var]] == level], na.rm = TRUE)
            }
        }, by = .(var, level)]
    } else {
        ## No events column if can't extract event variable
        all_terms_df[, Events := NA_integer_]
    }

    ## Apply condensing and indenting if requested
    if (condense_table || indent_groups) {
        if (condense_table) {
            indent_groups <- TRUE
        }
        
        all_terms_df[, inds := data.table::fifelse(level == "-", var, paste0(var, level))]
        orig_inds_map <- data.table::copy(all_terms_df[, .(var, level, inds, N, Events)])
        
        processed_rows <- list()
        row_counter <- 1
        unique_vars <- unique(all_terms_df[, var])
        
        for (v in unique_vars) {
            var_rows <- all_terms_df[var == v]
            
            if (nrow(var_rows) == 1) {
                processed_rows[[row_counter]] <- var_rows
                row_counter <- row_counter + 1
            } else {
                is_binary <- nrow(var_rows) == 2
                
                if (condense_table && is_binary) {
                    non_ref_row <- var_rows[2]
                    condensed_row <- data.table::copy(non_ref_row)
                    condensed_row[, var := paste0(v, " (", level, ")")]
                    condensed_row[, level := "-"]
                    processed_rows[[row_counter]] <- condensed_row
                    row_counter <- row_counter + 1
                } else {
                    if (indent_groups) {
                        header_row <- data.table::data.table(
                                                      var = v,
                                                      level = "-",
                                                      N = NA_integer_,
                                                      Events = NA_integer_,
                                                      pos = var_rows$pos[1],
                                                      var_order = var_rows$var_order[1],
                                                      inds = NA_character_
                                                  )
                        processed_rows[[row_counter]] <- header_row
                        row_counter <- row_counter + 1
                        
                        for (i in 1:nrow(var_rows)) {
                            group_row <- data.table::copy(var_rows[i])
                            group_row[, var := paste0("    ", level)]
                            group_row[, level := "-"]
                            processed_rows[[row_counter]] <- group_row
                            row_counter <- row_counter + 1
                        }
                    } else {
                        for (i in 1:nrow(var_rows)) {
                            processed_rows[[row_counter]] <- var_rows[i]
                            row_counter <- row_counter + 1
                        }
                    }
                }
            }
        }
        
        all_terms_df <- data.table::rbindlist(processed_rows, fill = TRUE)
        
        for (i in 1:nrow(all_terms_df)) {
            current_var <- all_terms_df$var[i]
            
            if (is.na(all_terms_df$inds[i]) && !grepl("^    ", current_var) && !grepl("\\(", current_var)) {
                next
            }
            
            if (grepl("\\(", current_var)) {
                orig_var <- gsub(" \\(.*\\)", "", current_var)
                orig_level <- gsub(".*\\((.*)\\)", "\\1", current_var)
                
                matching <- orig_inds_map[var == orig_var & level == orig_level]
                if (nrow(matching) > 0) {
                    all_terms_df[i, `:=`(inds = matching$inds[1],
                                         N = matching$N[1],
                                         Events = matching$Events[1])]
                }
            }
            else if (grepl("^    ", current_var)) {
                clean_level <- gsub("^    ", "", current_var)
                
                parent_var <- NA_character_
                for (j in (i-1):1) {
                    if (!grepl("^    ", all_terms_df$var[j]) && all_terms_df$var[j] != "") {
                        parent_var <- gsub(" \\(.*\\)", "", all_terms_df$var[j])
                        break
                    }
                }
                
                if (!is.na(parent_var)) {
                    matching <- orig_inds_map[var == parent_var & level == clean_level]
                    if (nrow(matching) > 0) {
                        all_terms_df[i, `:=`(inds = matching$inds[1],
                                             N = matching$N[1],
                                             Events = matching$Events[1])]
                    }
                }
            }
        }
    } else {
        all_terms_df[, inds := data.table::fifelse(level == "-", var, paste0(var, level))]
    }
    
    ## Process coefficients
    coef[, term := gsub(term, pattern = "`", replacement = "-")]
    coef[, inds := term]
    
    ## Merge data
    to_show <- merge(all_terms_df, coef, by.x = "inds", by.y = "inds", all.x = TRUE, sort = FALSE)
    
    ## Sort by variable order first, then position within variable
    data.table::setorder(to_show, var_order, pos)
    
    ## Add variable-based shading indicator (zebra stripes)
    if (zebra_stripes) {
        to_show[, shade_group := var_order %% 2]
        shade_colors <- c("#FFFFFF", "#EEEEEE")
    } else {
        to_show[, shade_group := 0]
        shade_colors <- c("#FFFFFF", "#FFFFFF")
    }
    
    ## Select columns
    to_show <- to_show[, .(var, level, N, Events, p_value, estimate, conf_low, conf_high, pos, var_order, shade_group)]
    
    ## Format the exponential values
    to_show_exp_clean <- data.table::copy(to_show)
    
    ## Create formatted columns for display
    to_show_exp_clean[, hr := data.table::fifelse(is.na(estimate), 
                                                  NA_real_,
                                     exp(estimate))]

                                        # For header rows (with NA N values), show empty strings instead of "reference"
    to_show_exp_clean[, hr_formatted := data.table::fifelse(is.na(N) & is.na(estimate),
                                                            "",  # Empty for header rows
                                                            data.table::fifelse(is.na(estimate), 
                                               ref_label,
                                               format(round(exp(estimate), digits), nsmall = digits)))]

    to_show_exp_clean[, conf_low_formatted := data.table::fifelse(is.na(conf_low), 
                                                     NA_character_,
                                                     format(round(exp(conf_low), digits), nsmall = digits))]
    to_show_exp_clean[, conf_high_formatted := data.table::fifelse(is.na(conf_high), 
                                                      NA_character_,
                                                      format(round(exp(conf_high), digits), nsmall = digits))]

    ## Format p-values
    to_show_exp_clean[, p_formatted := data.table::fifelse(is.na(p_value), 
                                              NA_character_,
                                       data.table::fifelse(p_value < 0.001, 
                                              "< 0.001",
                                              format(round(p_value, 3), nsmall = 3)))]

    ## Create the combined HR string with expression for italic p
    to_show_exp_clean[, hr_string_expr := data.table::fifelse(
                            is.na(N) & is.na(estimate),
                            "''",  # Empty string for header rows
                                          data.table::fifelse(
                                              is.na(estimate),
                                              paste0("'", ref_label, "'"),
                                          data.table::fifelse(p_value < 0.001,
                                                 paste0("'", hr_formatted, " (", conf_low_formatted, "-", 
                                                        conf_high_formatted, "); '*~italic(p)~'< 0.001'"),
                                                 paste0("'", hr_formatted, " (", conf_low_formatted, "-", 
                                                        conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'"))
                                          )
                        )]

    ## Format N and events with thousands separator
    to_show_exp_clean[, n_formatted := data.table::fifelse(is.na(N), "", format(N, big.mark = ",", scientific = FALSE))]
    to_show_exp_clean[, events_formatted := data.table::fifelse(is.na(Events), "", format(Events, big.mark = ",", scientific = FALSE))]
    
    ## Clean up variable names for display
    to_show_exp_clean[, var_display := as.character(var)]

    if (indent_groups || condense_table) {
        to_show_exp_clean[, var_display := var]
        
        for (v in unique(to_show_exp_clean$var)) {
            if (v != "" && !grepl("^    ", v)) {
                clean_v <- gsub(" \\(.*\\)", "", v)

                label <- if (!is.null(labels) && clean_v %in% names(labels)) {
                             labels[clean_v]
                         } else if (!is.null(attr(data[[clean_v]], "label"))) {
                             attr(data[[clean_v]], "label")
                         } else {
                             NULL
                         }

                if (!is.null(label)) {
                    if (grepl("\\(", v)) {
                        category <- gsub(".*\\((.*)\\)", "\\1", v)
                        if (category %in% c("1", "Yes", "yes")) {
                            to_show_exp_clean[var == v, var_display := label]
                        } else {
                            to_show_exp_clean[var == v, var_display := paste0(label, " (", category, ")")]
                        }
                    } else {
                        to_show_exp_clean[var == v, var_display := label]
                    }
                }
            }
        }
        
        to_show_exp_clean[, level := ""]
        
    } else {
        
        for(v in unique(to_show_exp_clean$var)) {
            if(v %in% to_show_exp_clean$var) {
                if(!is.null(labels) && v %in% names(labels)) {
                    to_show_exp_clean[var == v, var_display := labels[v]]
                }
                else if(!is.null(attr(data[[v]], "label"))) {
                    to_show_exp_clean[var == v, var_display := attr(data[[v]], "label")]
                }
            }
        }
    }

    if (!indent_groups) {
        to_show_exp_clean[duplicated(var), var_display := ""]
    }
    
    ## Handle missing estimates for plotting
    to_show_exp_clean[is.na(estimate), estimate := 0]
    
    ## Reorder (flip) - but maintain the variable grouping
    to_show_exp_clean <- to_show_exp_clean[order(nrow(to_show_exp_clean):1)]
    to_show_exp_clean[, x_pos := .I]
    
    ## Calculate plot ranges with better handling of extreme cases
    rangeb <- range(to_show$conf_low, to_show$conf_high, na.rm = TRUE)
    
    min_ci <- min(to_show$conf_low, na.rm = TRUE)
    max_ci <- max(to_show$conf_high, na.rm = TRUE)
    
    is_one_sided <- (min_ci > 0) || (max_ci < 0)
    
    ## Intelligent tick selection to prevent overlap
    range_magnitude <- diff(rangeb)
    
    if (exp(min_ci) < 0.01 && exp(max_ci) > 2) {
        ## Very wide range
        breaks <- c(0.01, 0.1, 0.5, 1, 2, 5)
    } else if (range_magnitude > 3) {
        ## Wide range - thin out the ticks
        all_breaks <- grDevices::axisTicks(rangeb/2, log = TRUE, nint = 7)
        if (length(all_breaks) > 7) {
            ## Too many - keep key values only
            important <- c(1)
            if (min(all_breaks) < 0.5) important <- c(min(all_breaks), important)
            if (max(all_breaks) > 2) important <- c(important, max(all_breaks))
            
            other_breaks <- setdiff(all_breaks, important)
            if (length(other_breaks) > 3) {
                keep_idx <- round(seq(1, length(other_breaks), length.out = 3))
                other_breaks <- other_breaks[keep_idx]
            }
            breaks <- sort(unique(c(important, other_breaks)))
        } else {
            breaks <- all_breaks
        }
    } else {
        ## Normal range - use standard calculation
        breaks <- grDevices::axisTicks(rangeb/2, log = TRUE, nint = 7)
    }
    
    if (!1 %in% breaks) {
        breaks <- sort(unique(c(breaks, 1)))
    }
    
    if (min_ci > 0) {
        rangeb[1] <- log(0.9)
    } else if(max_ci < 0) {
        rangeb[2] <- log(1.1)
    }
    
    breaks <- breaks[breaks >= exp(rangeb[1]) & breaks <= exp(rangeb[2])]
    reference_value <- 1

    ## Calculate layout using helper function
    layout <- calculate_forest_layout(
        to_show_exp_clean = to_show_exp_clean,
        show_n = show_n,
        show_events = show_events,
        indent_groups = indent_groups,
        condense_table = condense_table,
        effect_label = effect_label,
        ref_label = ref_label,
        font_size = font_size,
        tbl_width = data.table::fifelse(is.null(tbl_width), 0.6, tbl_width),
        rangeb = rangeb,
        center_padding = center_padding
    )

    ## Set up the extended range for plotting
    rangeplot <- c(layout$rangeplot_start, rangeb[2] + diff(rangeb) * 0.05)

    ## Extract positions
    y_variable <- layout$positions$var
    if (!(indent_groups || condense_table)) {
        y_level <- layout$positions$level
    }
    if (show_n) {
        y_n <- layout$positions$n
    }
    if (show_events) {
        y_events <- layout$positions$events
    }
    y_hr <- layout$positions$effect

    ## Use the effect abbreviation from layout
    effect_abbrev <- layout$effect_abbrev

    ## Calculate recommended dimensions
    rec_height <- max(5, min(20, 3 + nrow(to_show_exp_clean) * 0.25))

    if(!is.null(plot_width)) {
        rec_width <- plot_width
        if(!is.null(plot_height)) {
            rec_height <- plot_height
        }
    } else {
        ## Use the calculated total width
        rec_width <- layout$total_width + 1.0  # Add margins
        rec_width <- max(10, min(20, rec_width))  # Apply reasonable bounds
    }

    ## Font sizes
    annot_font <- font_size * annot_size
    header_font <- font_size * header_size
    
    ## Custom ticks data
    ticks_df <- data.frame(
        x = -0.5,
        xend = -0.7,
        y = breaks,
        yend = breaks
    )

    ## Create the plot
    p <- ggplot2::ggplot(to_show_exp_clean, ggplot2::aes(x_pos, exp(estimate))) +
        
        ## Shading rectangles - extend to cover table area too
        ggplot2::geom_rect(ggplot2::aes(xmin = x_pos - .5, xmax = x_pos + .5,
                                        ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                                        fill = ordered(shade_group))) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
        ggplot2::scale_size_continuous(range = c(1, 6), guide = "none") +
        ggplot2::scale_fill_manual(values = shade_colors, guide = "none") +
        
        ## Forest plot elements
        ggplot2::geom_point(ggplot2::aes(size = N), pch = 22, color = "#000000", fill = color, na.rm = TRUE) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = exp(conf_low), ymax = exp(conf_high)), width = 0.15) +
        
        ## Y-axis for forest plot
        ggplot2::annotate(geom = "segment",
                          x = -0.5, xend = -0.5,
                          y = exp(rangeb[1]),
                          yend = exp(rangeb[2]),
                          color = "#000000", linewidth = 1) +
        
        ## Reference line at HR = 1
        ggplot2::annotate(geom = "segment", 
                          x = -0.5, xend = max(to_show_exp_clean$x_pos) + 0.5, 
                          y = 1, yend = 1, 
                          linetype = "longdash") +
        
        ## Ticks
        ggplot2::geom_segment(data = ticks_df,
                              ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                              inherit.aes = FALSE,
                              color = "#000000",
                              linewidth = 1) +
        ggplot2::geom_text(data = ticks_df,
                           ggplot2::aes(x = xend - 0.05, y = y, label = sprintf("%g", y)),
                           inherit.aes = FALSE,
                           hjust = 0.5,
                           vjust = 1.3,
                           size = annot_font * 1.5) +
        
        ## Set coordinate system with extended limits
        ggplot2::coord_flip(ylim = exp(rangeplot)) +
        ggplot2::ggtitle(title) +
        ggplot2::scale_y_log10(name = "Hazard Ratio",
                               labels = sprintf("%g", breaks),
                               expand = c(0.02, 0.02),
                               breaks = breaks) +
        ggplot2::theme_light() +
        ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0),
                       panel.grid.minor.y = ggplot2::element_blank(),
                       panel.grid.minor.x = ggplot2::element_blank(),
                       panel.grid.major.y = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       legend.position = "none",
                       panel.border = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size = font_size * title_size, face = "bold", hjust = 0.5)) +
        ggplot2::xlab("") +
        
        ## Variable column
        ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = exp(y_variable),
                          label = "Variable", fontface = "bold", hjust = 0,
                          size = header_font) +

    {if (indent_groups || condense_table) {
                                        # When indented/condensed, use conditional formatting
         ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = exp(y_variable),
                           label = to_show_exp_clean$var_display, 
                           fontface = data.table::fifelse(grepl("^    ", to_show_exp_clean$var_display), "plain", "bold"), 
                           hjust = 0,
                           size = annot_font)
     } else {
                                        # Original formatting
         ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = exp(y_variable),
                           label = to_show_exp_clean$var_display, fontface = "bold", hjust = 0,
                           size = annot_font)
     }} +
    
    ## Group/level column
    
    {if (!(indent_groups || condense_table)) {
         list(
             ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = exp(y_level),
                               label = "Group", fontface = "bold", hjust = 0,
                               size = header_font),
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = exp(y_level),
                               label = to_show_exp_clean$level, hjust = 0,
                               size = annot_font)
         )
     }} +
    
    ## N column (conditional)
    {if (show_n) {
         list(
             ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = exp(y_n),
                               label = "n", fontface = "bold.italic", hjust = 0.5,
                               size = header_font),
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = exp(y_n),
                               label = to_show_exp_clean$n_formatted, hjust = 0.5,
                               size = annot_font)
         )
     }} +
    
    ## Events column (conditional)
    {if (show_events) {
         list(
             ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = exp(y_events),
                               label = "Events", fontface = "bold", hjust = 0.5,
                               size = header_font),
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = exp(y_events),
                               label = to_show_exp_clean$events_formatted, hjust = 0.5,
                               size = annot_font)
         )
     }} +
    
    ## Effect column
    ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.4, y = exp(y_hr),
                      label = "bold('aHR (95% CI); '*bolditalic(p)*'-value')",
                      hjust = 0, size = header_font, parse = TRUE) +
    
    ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = exp(y_hr),
                      label = to_show_exp_clean$hr_string_expr, hjust = 0,
                      size = annot_font, parse = TRUE) +
    
    ## X-axis label
    ggplot2::annotate(geom = "text", x = -1.5, y = 1,
                      label = "Hazard Ratio", fontface = "bold",
                      hjust = 0.5, vjust = 2, size = annot_font * 1.5) +
    
    ## Model statistics footer
    ggplot2::annotate(geom = "text", x = 0.5, y = exp(y_variable),
                      label = paste0("Events analyzed: ", gmodel$nevent_with_pct,
                                     "\nGlobal log-rank p: ", global_p_formatted,
                                     "\n", concordance_string,
                                     "\nAIC: ", gmodel$AIC_formatted),
                      size = annot_font, hjust = 0, vjust = 1.2, fontface = "italic")
    
    ## Convert units back for output if needed
    if (units != "in") {
        rec_width <- convert_units(rec_width, from = "in", to = units)
        rec_height <- convert_units(rec_height, from = "in", to = units)
    }

    ## Provide dimension recommendations
    if(is.null(plot_width) || is.null(plot_height)) {
        message(sprintf("Recommended plot dimensions: width = %.1f %s, height = %.1f %s",
                        rec_width, units, rec_height, units))
    }
    
    ## Add recommended dimensions as an attribute
    attr(p, "recommended_dims") <- list(width = rec_width, height = rec_height)

    ## Return the plot
    return(p)
}
