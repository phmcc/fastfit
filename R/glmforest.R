#' Create Forest Plot for Generalized Linear Models
#'
#' Generates a publication-ready forest plot that combines a formatted data table 
#' with a graphical representation of effect estimates (odds ratios, risk ratios, 
#' or coefficients) from a generalized linear model. The plot integrates variable 
#' names, group levels, sample sizes, effect estimates with confidence intervals, 
#' p-values, and model diagnostics in a single comprehensive visualization designed 
#' for manuscripts and presentations.
#'
#' @param model A fitted generalized linear model object of class \code{glm}. 
#'   The model should be fitted using \code{\link[stats]{glm}} with appropriate 
#'   family specification (binomial, poisson, gaussian, etc.).
#'   
#' @param data A data.frame or data.table containing the original data used to 
#'   fit the model. If \code{NULL} (default), the function attempts to extract 
#'   data from the model object itself (via \code{model$data} or \code{model$model}). 
#'   Providing the data explicitly is recommended for reliability.
#'   
#' @param title Character string specifying the plot title displayed at the top. 
#'   Default is \code{"Generalized Linear Model"}. Use descriptive titles like 
#'   "Risk Factors for Disease Outcome" for publication.
#'   
#' @param effect_label Character string for the effect measure label on the 
#'   forest plot axis. If \code{NULL} (default), automatically determined based 
#'   on model family and link function: "Odds Ratio" for logistic regression 
#'   (\code{family = binomial, link = logit}), "Risk Ratio" for log-link models, 
#'   "Exp(Coefficient)" for other exponential families, or "Coefficient" for 
#'   identity link.
#'   
#' @param digits Integer specifying the number of decimal places for effect 
#'   estimates and confidence intervals in the data table. Default is 2.
#'   
#' @param font_size Numeric multiplier controlling the base font size for all 
#'   text elements. Values > 1 increase all fonts proportionally, values < 1 
#'   decrease them. Default is 1.0. Useful for adjusting readability across 
#'   different output sizes.
#'   
#' @param annot_size Numeric value controlling the relative font size for 
#'   data annotations (variable names, values in table cells). Default is 3.88. 
#'   Adjust relative to \code{font_size}.
#'   
#' @param header_size Numeric value controlling the relative font size for 
#'   column headers ("Variable", "Group", "n", etc.). Default is 5.82. Headers 
#'   are typically larger than annotations for hierarchy.
#'   
#' @param title_size Numeric value controlling the relative font size for the 
#'   main plot title. Default is 23.28. The title is typically the largest text 
#'   element.
#'   
#' @param tbl_width Numeric value between 0 and 1 specifying the proportion of 
#'   total plot width allocated to the data table (left side). The forest plot 
#'   occupies \code{1 - tbl_width}. Default is 0.6 (60\% table, 40\% forest). 
#'   Increase for longer variable names, decrease to emphasize the forest plot.
#'   
#' @param plot_width Numeric value specifying the intended output width in 
#'   specified \code{units}. Used for optimizing layout and text sizing. 
#'   Default is \code{NULL} (automatic). Recommended: 10-16 inches for standard 
#'   publications.
#'   
#' @param plot_height Numeric value specifying the intended output height in 
#'   specified \code{units}. Default is \code{NULL} (automatic based on number 
#'   of rows). The function provides recommendations if not specified.
#'   
#' @param show_n Logical. If \code{TRUE}, includes a column showing group-specific 
#'   sample sizes for categorical variables and total sample size for continuous 
#'   variables. Default is \code{TRUE}.
#'   
#' @param show_events Logical. If \code{TRUE}, includes a column showing the 
#'   number of events for each group. Relevant for logistic regression (number 
#'   of cases) and other binary outcomes. Default is \code{TRUE}.
#'   
#' @param indent_groups Logical. If \code{TRUE}, indents factor levels under 
#'   their parent variable name, creating a hierarchical visual structure. When 
#'   \code{TRUE}, the "Group" column is hidden. Default is \code{FALSE}.
#'   
#' @param condense_table Logical. If \code{TRUE}, condenses binary categorical 
#'   variables into single rows by showing only the non-reference category. 
#'   Automatically sets \code{indent_groups = TRUE}. Useful for tables with 
#'   many binary variables. Default is \code{FALSE}.
#'   
#' @param center_padding Numeric value specifying the horizontal spacing (in 
#'   character units) between the data table and forest plot. Increase for more 
#'   separation, decrease to fit more content. Default is 4.
#'   
#' @param zebra_stripes Logical. If \code{TRUE}, applies alternating gray 
#'   background shading to different variables (not rows) to improve visual 
#'   grouping and readability. Default is \code{TRUE}.
#'   
#' @param ref_label Character string to display for reference categories of 
#'   factor variables. Typically shown in place of effect estimates. 
#'   Default is \code{"reference"}. Common alternatives: "ref", "1.00 (ref)".
#'   
#' @param labels Named character vector or list providing custom display 
#'   labels for variables. Names should match variable names in the model, 
#'   values are the labels to display. Example: 
#'   \code{c(age = "Age (years)", bmi = "Body Mass Index")}. Default is \code{NULL} 
#'   (use original variable names).
#'   
#' @param color Character string specifying the color for effect estimate point 
#'   markers in the forest plot. Use hex codes or R color names. Default is 
#'   \code{"#3A7E8A"} (teal). Choose colors that contrast well with black 
#'   error bars.
#'   
#' @param exponentiate Logical. If \code{TRUE}, exponentiates coefficients to 
#'   display odds ratios, risk ratios, etc. If \code{FALSE}, shows raw 
#'   coefficients. Default is \code{NULL}, which automatically exponentiates 
#'   for logit, log, and cloglog links, and shows raw coefficients for identity 
#'   link.
#'   
#' @param units Character string specifying the units for plot dimensions. 
#'   Options: \code{"in"} (inches), \code{"cm"} (centimeters), \code{"mm"} 
#'   (millimeters). Default is \code{"in"}. Affects interpretation of 
#'   \code{plot_width} and \code{plot_height}.
#'
#' @return A \code{ggplot} object containing the complete forest plot. The plot 
#'   can be:
#'   \itemize{
#'     \item Displayed directly: \code{print(plot)}
#'     \item Saved to file: \code{ggsave("forest.pdf", plot, width = 12, height = 8)}
#'     \item Further customized with ggplot2 functions
#'   }
#'   
#'   The returned object includes an attribute \code{"recommended_dims"} 
#'   accessible via \code{attr(plot, "recommended_dims")}, which is a list 
#'   containing:
#'   \describe{
#'     \item{width}{Numeric. Recommended plot width in specified units}
#'     \item{height}{Numeric. Recommended plot height in specified units}
#'   }
#'   
#'   These recommendations are automatically calculated based on the number of 
#'   variables, text sizes, and layout parameters, and are printed to console 
#'   if \code{plot_width} or \code{plot_height} are not specified.
#'
#' @details
#' \strong{Plot Components:}
#' 
#' The forest plot consists of several integrated components:
#' \enumerate{
#'   \item \strong{Title}: Centered at top, describes the analysis
#'   \item \strong{Data Table} (left side): Contains columns for:
#'     \itemize{
#'       \item Variable: Predictor names (or custom labels)
#'       \item Group: Factor levels (optional, hidden when indenting)
#'       \item n: Sample sizes by group (optional)
#'       \item Events: Event counts by group (optional)
#'       \item Effect (95\% CI); \emph{p}-value: Formatted estimates with p-values
#'     }
#'   \item \strong{Forest Plot} (right side): Graphical display with:
#'     \itemize{
#'       \item Point estimates (squares sized by sample size)
#'       \item 95\% confidence intervals (error bars)
#'       \item Reference line (at OR/RR = 1 or coefficient = 0)
#'       \item Log scale for odds/risk ratios
#'       \item Labeled axis
#'     }
#'   \item \strong{Model Statistics} (footer): Summary of:
#'     \itemize{
#'       \item Observations analyzed (with percentage of total data)
#'       \item Model family (Binomial, Poisson, etc.)
#'       \item Deviance statistics
#'       \item Pseudo-R^2 (McFadden)
#'       \item AIC
#'     }
#' }
#' 
#' \strong{Automatic Effect Measure Selection:}
#' 
#' When \code{effect_label = NULL} and \code{exponentiate = NULL}, the function 
#' intelligently selects the appropriate effect measure:
#' \itemize{
#'   \item \strong{Logistic regression} (\code{family = binomial(link = "logit")}): 
#'     Odds Ratios (OR)
#'   \item \strong{Log-link models} (\code{link = "log"}): Risk Ratios (RR) 
#'     or Rate Ratios
#'   \item \strong{Other exponential families}: Exp(Coefficient)
#'   \item \strong{Identity link}: Raw coefficients
#' }
#' 
#' \strong{Reference Categories:}
#' 
#' For factor variables, the first level (determined by factor ordering or 
#' alphabetically for character variables) serves as the reference category:
#' \itemize{
#'   \item Displayed with the \code{ref_label} instead of an estimate
#'   \item No confidence interval or p-value shown
#'   \item Visually aligned with other categories
#'   \item When \code{condense_table = TRUE}, reference-only variables may be 
#'     omitted entirely
#' }
#' 
#' \strong{Layout Optimization:}
#' 
#' The function automatically optimizes layout based on content:
#' \itemize{
#'   \item Calculates appropriate axis ranges to accommodate all confidence intervals
#'   \item Selects meaningful tick marks on log or linear scales
#'   \item Sizes point markers proportional to sample size (larger = more data)
#'   \item Adjusts table width based on variable name lengths when \code{tbl_width = NULL}
#'   \item Recommends overall dimensions based on number of rows
#' }
#' 
#' \strong{Visual Grouping Options:}
#' 
#' Three display modes are available:
#' \enumerate{
#'   \item \strong{Standard} (\code{indent_groups = FALSE}, 
#'     \code{condense_table = FALSE}): 
#'     Separate "Variable" and "Group" columns, all categories shown
#'   \item \strong{Indented} (\code{indent_groups = TRUE}, 
#'     \code{condense_table = FALSE}): 
#'     Hierarchical display with groups indented under variables
#'   \item \strong{Condensed} (\code{condense_table = TRUE}): 
#'     Binary variables shown in single rows, automatically indented
#' }
#' 
#' \strong{Zebra Striping:}
#' 
#' When \code{zebra_stripes = TRUE}, alternating variables (not individual rows) 
#' receive light gray backgrounds. This helps visually group all levels of a 
#' factor variable together, making the plot easier to read especially with 
#' many multi-level factors.
#' 
#' \strong{Model Statistics Display:}
#' 
#' The footer shows key diagnostic information:
#' \itemize{
#'   \item \strong{Observations analyzed}: Total N and percentage of original 
#'     data (accounting for missing values)
#'   \item \strong{Null/Residual Deviance}: Model fit improvement
#'   \item \strong{Pseudo-R^2}: McFadden R^2 = 1 - (log L_1 / log L_2)
#'   \item \strong{AIC}: For model comparison (lower is better)
#' }
#' 
#' For logistic regression, concordance (C-statistic/AUC) may also be displayed 
#' if available.
#' 
#' \strong{Saving Plots:}
#' 
#' Use \code{ggplot2::ggsave()} with recommended dimensions:
#' ```r
#' p <- glmforest(model, data)
#' dims <- attr(p, "recommended_dims")
#' ggsave("forest.pdf", p, width = dims$width, height = dims$height)
#' ```
#' 
#' Or specify custom dimensions:
#' ```r
#' ggsave("forest.png", p, width = 12, height = 8, dpi = 300)
#' ```
#'
#' @seealso 
#' \code{\link{coxforest}} for Cox proportional hazards forest plots,
#' \code{\link{lmforest}} for linear model forest plots,
#' \code{\link[stats]{glm}} for fitting GLMs,
#' \code{\link{fit}} for fastfit regression modeling,
#' \code{\link[ggplot2]{ggsave}} for saving plots
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' # Example 1: Basic logistic regression forest plot
#' model1 <- glm(os_status ~ age + sex + bmi + treatment,
#'               data = clintrial,
#'               family = binomial)
#' 
#' plot1 <- glmforest(model1, data = clintrial)
#' print(plot1)
#' 
#' # Example 2: With custom variable labels
#' plot2 <- glmforest(
#'     model = model1,
#'     data = clintrial,
#'     title = "Risk Factors for Mortality",
#'     labels = clintrial_labels
#' )
#' print(plot2)
#' 
#' # Example 3: Customize appearance
#' plot3 <- glmforest(
#'     model = model1,
#'     data = clintrial,
#'     title = "Adjusted Odds Ratios",
#'     color = "#D62728",  # Red points
#'     font_size = 1.2,    # Larger text
#'     zebra_stripes = FALSE,
#'     labels = clintrial_labels
#' )
#' print(plot3)
#' 
#' # Example 4: Indented layout for hierarchical view
#' plot4 <- glmforest(
#'     model = model1,
#'     data = clintrial,
#'     indent_groups = TRUE,
#'     labels = clintrial_labels
#' )
#' print(plot4)
#' # Group column hidden, levels indented under variables
#' 
#' # Example 5: Condensed layout for many binary variables
#' model5 <- glm(os_status ~ age + sex + smoking + hypertension + 
#'                   diabetes + surgery,
#'               data = clintrial,
#'               family = binomial)
#' 
#' plot5 <- glmforest(
#'     model = model5,
#'     data = clintrial,
#'     condense_table = TRUE,
#'     labels = clintrial_labels
#' )
#' print(plot5)
#' # Binary variables shown in single rows
#' 
#' # Example 6: Hide sample size and events columns
#' plot6 <- glmforest(
#'     model = model1,
#'     data = clintrial,
#'     show_n = FALSE,
#'     show_events = FALSE,
#'     labels = clintrial_labels
#' )
#' print(plot6)
#' 
#' # Example 7: Adjust table/forest proportions
#' plot7 <- glmforest(
#'     model = model1,
#'     data = clintrial,
#'     tbl_width = 0.7,  # More space for table
#'     labels = clintrial_labels
#' )
#' print(plot7)
#' 
#' # Example 8: Get and use recommended dimensions
#' plot8 <- glmforest(model1, data = clintrial)
#' 
#' dims <- attr(plot8, "recommended_dims")
#' cat("Recommended: ", dims$width, "x", dims$height, "inches\n")
#' 
#' # Save with recommended dimensions
#' # ggsave("forest.pdf", plot8, width = dims$width, height = dims$height)
#' 
#' # Example 9: Specify exact output dimensions
#' plot9 <- glmforest(
#'     model = model1,
#'     data = clintrial,
#'     plot_width = 14,
#'     plot_height = 10,
#'     labels = clintrial_labels
#' )
#' # No dimension recommendations printed
#' 
#' # Example 10: Use different units (centimeters)
#' plot10 <- glmforest(
#'     model = model1,
#'     data = clintrial,
#'     plot_width = 30,  # 30 cm
#'     plot_height = 20,  # 20 cm
#'     units = "cm",
#'     labels = clintrial_labels
#' )
#' 
#' # Example 11: Poisson regression for count data
#' model11 <- glm(los_days ~ age + treatment + surgery + stage,
#'                data = clintrial,
#'                family = poisson)
#' 
#' plot11 <- glmforest(
#'     model = model11,
#'     data = clintrial,
#'     title = "Rate Ratios for Length of Stay",
#'     effect_label = "Rate Ratio",
#'     labels = clintrial_labels
#' )
#' print(plot11)
#' # Shows rate ratios instead of odds ratios
#' 
#' # Example 12: Force display of raw coefficients
#' plot12 <- glmforest(
#'     model = model1,
#'     data = clintrial,
#'     exponentiate = FALSE,  # Show log odds
#'     effect_label = "Log Odds",
#'     labels = clintrial_labels
#' )
#' print(plot12)
#' # Reference line at 0 instead of 1
#' 
#' # Example 13: Custom reference label
#' plot13 <- glmforest(
#'     model = model1,
#'     data = clintrial,
#'     ref_label = "1.00 (ref)",
#'     labels = clintrial_labels
#' )
#' print(plot13)
#' 
#' # Example 14: Adjust font sizes for presentations
#' plot14 <- glmforest(
#'     model = model1,
#'     data = clintrial,
#'     font_size = 1.5,      # 50% larger
#'     title_size = 30,      # Larger title
#'     labels = clintrial_labels
#' )
#' print(plot14)
#' 
#' # Example 15: Complete publication-ready plot
#' final_model <- glm(
#'     os_status ~ age + sex + bmi + smoking + hypertension + 
#'         diabetes + treatment + stage,
#'     data = clintrial,
#'     family = binomial
#' )
#' 
#' final_plot <- glmforest(
#'     model = final_model,
#'     data = clintrial,
#'     title = "Multivariable Logistic Regression: Risk Factors for Mortality",
#'     labels = clintrial_labels,
#'     indent_groups = TRUE,
#'     zebra_stripes = TRUE,
#'     color = "#3A7E8A",
#'     font_size = 1.0,
#'     digits = 2
#' )
#' 
#' # Save for publication
#' dims <- attr(final_plot, "recommended_dims")
#' # ggsave("figure1.pdf", final_plot, 
#' #        width = dims$width, height = dims$height)
#' # ggsave("figure1.png", final_plot,
#' #        width = dims$width, height = dims$height, dpi = 300)
#'
#' @export
glmforest <- function(model, data = NULL,
                      title = "Generalized Linear Model",
                      effect_label = NULL,
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
                      color = "#3A7E8A",
                      exponentiate = NULL,
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
    
    ## Check model class - support both GLM and GLMER
    is_lme4 <- inherits(model, c("glmerMod", "merMod"))
    is_glm <- inherits(model, "glm")
    
    if (!is_glm && !is_lme4) {
        stop("Model must be either a glm or glmerMod/merMod object")
    }
    
    ## Internally work in inches
    if (!is.null(plot_width) && units != "in") {
        plot_width <- convert_units(plot_width, from = units, to = "in")
    }

    ## Determine if we should exponentiate based on link function
    if(is.null(exponentiate)) {
        if (is_lme4) {
            ## For lme4 models, extract family info differently
            family_obj <- model@resp$family
            link_func <- family_obj$link
        } else {
            ## For regular GLMs
            link_func <- model$family$link
        }
        exponentiate <- link_func %in% c("logit", "log", "cloglog")
    }
    
    ## Set effect label based on model family and link
    if(is.null(effect_label)) {
        if (is_lme4) {
            family_obj <- model@resp$family
            family_name <- family_obj$family
            link_func <- family_obj$link
        } else {
            family_name <- model$family$family
            link_func <- model$family$link
        }
        
        if(family_name == "binomial" && link_func == "logit") {
            effect_label <- "Odds Ratio"
        } else if(link_func == "log") {
            effect_label <- "Risk Ratio"
        } else if(exponentiate) {
            effect_label <- "Exp(Coefficient)"
        } else {
            effect_label <- "Coefficient"
        }
    }
    
    ## Get data - we need both the original data and the model data
    ## The model data excludes NAs and is what we should use for counts
    if (is_lme4) {
        ## For lme4 models, the model frame contains the complete cases
        model_data <- model@frame
    } else {
        ## For regular GLMs, use model$model which contains complete cases
        model_data <- model$model
    }
    
    ## Get the original data for reference (if provided)
    if(is.null(data)){
        warning("The `data` argument is not provided. Data will be extracted from model fit.")
        data <- model_data  # Use the model data as the base
    }
    
    ## Convert to data.table if not already
    if(!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    if(!data.table::is.data.table(model_data)) {
        model_data <- data.table::as.data.table(model_data)
    }
    
    ## Extract terms based on model type
    if (is_lme4) {
        ## For lme4 models, we need to work with the fixed effects only
        ## Get the fixed effects formula
        fixed_formula <- lme4::nobars(stats::formula(model))
        
        ## Get variable names from the fixed formula (excluding outcome)
        var_names <- all.vars(fixed_formula)[-1]  # Remove outcome variable
        
        ## Get the model frame for determining variable types
        model_frame <- model@frame
        
        ## Create terms object similar to GLM structure
        terms <- sapply(var_names, function(v) {
            if (v %in% names(model_frame)) {
                if (is.factor(model_frame[[v]])) "factor"
                else if (is.character(model_frame[[v]])) "character"
                else "numeric"
            } else if (v %in% names(data)) {
                if (is.factor(data[[v]])) "factor"
                else if (is.character(data[[v]])) "character"
                else "numeric"
            } else "numeric"
        }, USE.NAMES = TRUE)
        
    } else {
        ## For regular GLMs
        terms <- attr(model$terms, "dataClasses")[-1]
    }
    
    ## Extract coefficients and confidence intervals
    if (is_lme4) {
        ## For lme4 models
        coef_summary <- summary(model)$coefficients
        
        ## Calculate confidence intervals manually for lme4
        ## Using Wald confidence intervals
        coef_vals <- lme4::fixef(model)
        se_vals <- sqrt(diag(as.matrix(vcov(model))))
        z_crit <- qnorm(0.975)
        
        conf_int <- cbind(
            lower = coef_vals - z_crit * se_vals,
            upper = coef_vals + z_crit * se_vals
        )
        rownames(conf_int) <- names(coef_vals)
        
        ## Extract z values - for lme4, we calculate them
        z_values <- coef_summary[, "Estimate"] / coef_summary[, "Std. Error"]
        p_values <- 2 * pnorm(abs(z_values), lower.tail = FALSE)
        
    } else {
        ## For regular GLMs
        coef_summary <- summary(model)$coefficients
        conf_int <- stats::confint(model)
        z_values <- coef_summary[, "z value"]
        p_values <- coef_summary[, "Pr(>|z|)"]
    }
    
    ## Remove intercept if present
    if("(Intercept)" %in% rownames(coef_summary)) {
        intercept_idx <- which(rownames(coef_summary) == "(Intercept)")
        coef_summary <- coef_summary[-intercept_idx, , drop = FALSE]
        conf_int <- conf_int[-intercept_idx, , drop = FALSE]
        z_values <- z_values[-intercept_idx]
        p_values <- p_values[-intercept_idx]
    }
    
    ## Create coefficient data table
    coef <- data.table::data.table(
                            term = rownames(coef_summary),
                            estimate = coef_summary[, "Estimate"],
                            std_error = coef_summary[, "Std. Error"],
                            statistic = z_values,
                            p_value = p_values,
                            conf_low = conf_int[, 1],
                            conf_high = conf_int[, 2]
                        )
    
    ## Get model statistics
    if (is_lme4) {
        ## For lme4 models
        gmodel <- list(
            nobs = nrow(model@frame),
            null_deviance = NA,  # Not directly available in lme4
            residual_deviance = NA,  # Not directly available in lme4
            AIC = stats::AIC(model),
            family = paste0(model@resp$family$family, " (", model@resp$family$link, " link)")
        )
        
        ## Calculate deviances if possible
        tryCatch({
                                        # Some lme4 models have deviance method
            gmodel$residual_deviance <- stats::deviance(model)
        }, error = function(e) {
                                        # Keep as NA if not available
        })
        
    } else {
        ## For regular GLMs
        gmodel <- list(
            nobs = nobs(model),
            null_deviance = model$null.deviance,
            residual_deviance = model$deviance,
            AIC = stats::AIC(model),
            family = paste0(model$family$family, " (", model$family$link, " link)")
        )
    }
    
    ## Calculate total observations and percentage analyzed
    total_obs <- nrow(data)
    gmodel$pct_analyzed <- (gmodel$nobs / total_obs) * 100
    
    ## Format observations and AIC with commas
    gmodel$nobs_formatted <- format(gmodel$nobs, big.mark = ",", scientific = FALSE)
    gmodel$nobs_with_pct <- paste0(gmodel$nobs_formatted, " (", 
                                   sprintf("%.1f%%", gmodel$pct_analyzed), ")")
    gmodel$AIC_formatted <- format(round(gmodel$AIC, 2), big.mark = ",", scientific = FALSE, nsmall = 2)
    
    ## Calculate pseudo R-squared (McFadden's) - only if deviances available
    if (!is.na(gmodel$null_deviance) && !is.na(gmodel$residual_deviance)) {
        gmodel$pseudo_r2 <- 1 - (gmodel$residual_deviance / gmodel$null_deviance)
    } else {
        gmodel$pseudo_r2 <- NA
    }
    
    ## Extract xlevels (factor levels) based on model type
    if (is_lme4) {
        ## For lme4 models, extract from the model frame
        frame_data <- model@frame
        xlevels <- list()
        
        for (var_name in names(terms)) {
            if (terms[var_name] %in% c("factor", "character")) {
                if (var_name %in% names(frame_data) && is.factor(frame_data[[var_name]])) {
                    xlevels[[var_name]] <- levels(frame_data[[var_name]])
                } else if (var_name %in% names(model_data) && is.factor(model_data[[var_name]])) {
                    xlevels[[var_name]] <- levels(model_data[[var_name]])
                }
            }
        }
        
        if (length(xlevels) == 0) xlevels <- NULL
        
    } else {
        ## For regular GLMs
        xlevels <- model$xlevels
    }
    
    ## Extract statistics for every variable - preserving order
    all_terms <- lapply(seq_along(terms), function(i){
        var <- names(terms)[i]
        
        if (terms[i] %in% c("factor", "character")) {
            ## Get the factor levels from xlevels (proper order)
            if(!is.null(xlevels) && var %in% names(xlevels)) {
                factor_levels <- xlevels[[var]]
                
                ## Create data table with proper levels FROM MODEL DATA
                level_counts <- model_data[!is.na(get(var)), .N, by = var]
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
                return(all_levels_dt[, .(var, level, Freq, pos, var_order)])
            } else {
                ## Fallback for variables not in xlevels - use MODEL DATA
                adf <- model_data[!is.na(get(var)), .N, by = var]
                data.table::setnames(adf, old = c(var, "N"), new = c("level", "Freq"))
                adf[, var := var]
                adf[, pos := .I]
                adf[, var_order := i]
                return(adf[, .(var, level, Freq, pos, var_order)])
            }
        }
        else if (terms[i] == "numeric") {
            ## For numeric variables, return a single row - use MODEL DATA
            return(data.table::data.table(
                                   var = var, 
                                   level = "-", 
                                   Freq = sum(!is.na(model_data[[var]])), 
                                   pos = 1, 
                                   var_order = i
                               ))
        }
        else {
            ## Other cases
            vars = grep(paste0("^", var, "*."), coef$term, value=TRUE)
            return(data.table::data.table(
                                   var = vars, 
                                   level = "", 
                                   Freq = nrow(model_data),  # Use model_data row count
                                   pos = seq_along(vars), 
                                   var_order = i
                               ))
        }
    })
    
    ## Remove any NULL entries from the list
    all_terms <- all_terms[!sapply(all_terms, is.null)]
    
    ## Bind all terms together
    all_terms_df <- data.table::rbindlist(all_terms)
    data.table::setnames(all_terms_df, c("var", "level", "N", "pos", "var_order"))

    ## Add events for binomial models
    if (is_lme4) {
        family_name <- model@resp$family$family
    } else {
        family_name <- model$family$family
    }
    
    if (family_name == "binomial") {
        outcome_data <- NULL
        
        if (is_lme4) {
            ## For lme4 models, extract outcome from the response
            outcome_data <- model@resp$y
        } else {
            ## For regular GLMs, get outcome from model data
            outcome_var <- all.vars(stats::formula(model))[1]
            
            if (!outcome_var %in% names(model_data)) {
                warning(paste("Outcome variable", outcome_var, "not found in model data. Events column will not be created."))
            } else {
                outcome_data <- model_data[[outcome_var]]
                
                if (is.factor(outcome_data)) {
                    outcome_data <- as.numeric(outcome_data) == 2
                }
            }
        }
        
        if (!is.null(outcome_data)) {
            all_terms_df[, Events := {
                if (level == "-") {
                                        # For continuous variables, total events
                    sum(outcome_data, na.rm = TRUE)
                } else {
                                        # For factor levels, events within that level
                    if (var %in% names(model_data)) {
                        sum(outcome_data[model_data[[var]] == level], na.rm = TRUE)
                    } else {
                        NA_integer_
                    }
                }
            }, by = seq_len(nrow(all_terms_df))]
        } else {
            all_terms_df[, Events := NA_integer_]
        }
    } else {
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
        
        ## Restore inds, N, and Events values for condensed/indented rows
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
    coef[, term := gsub(term, pattern = "`", replacement = "")]
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

    print(to_show)

    ## Format the values based on exponentiate setting
    to_show_exp_clean <- data.table::copy(to_show)

    ## Create formatted columns for display
    if(exponentiate) {
        to_show_exp_clean[, effect := data.table::fifelse(is.na(estimate), 
                                                          NA_real_,
                                                          exp(estimate))]
        to_show_exp_clean[, effect_formatted := data.table::fifelse(is.na(N) & is.na(estimate),
                                                                    "",
                                                                    data.table::fifelse(is.na(estimate), 
                                                                                        ref_label,
                                                                                        format(round(exp(estimate), digits), nsmall = digits)))]
        to_show_exp_clean[, conf_low_formatted := data.table::fifelse(is.na(conf_low), 
                                                                      NA_character_,
                                                                      format(round(exp(conf_low), digits), nsmall = digits))]
        to_show_exp_clean[, conf_high_formatted := data.table::fifelse(is.na(conf_high), 
                                                                       NA_character_,
                                                                       format(round(exp(conf_high), digits), nsmall = digits))]
    } else {
        to_show_exp_clean[, effect := estimate]
        to_show_exp_clean[, effect_formatted := data.table::fifelse(is.na(N) & is.na(estimate),
                                                                    "",
                                                                    data.table::fifelse(is.na(estimate), 
                                                                                        ref_label,
                                                                                        format(round(estimate, digits), nsmall = digits)))]
        to_show_exp_clean[, conf_low_formatted := data.table::fifelse(is.na(conf_low), 
                                                                      NA_character_,
                                                                      format(round(conf_low, digits), nsmall = digits))]
        to_show_exp_clean[, conf_high_formatted := data.table::fifelse(is.na(conf_high), 
                                                                       NA_character_,
                                                                       format(round(conf_high, digits), nsmall = digits))]
    }

    ## Format p-values
    to_show_exp_clean[, p_formatted := data.table::fifelse(is.na(p_value), 
                                                           NA_character_,
                                                           data.table::fifelse(p_value < 0.001, 
                                                                               "< 0.001",
                                                                               format(round(p_value, 3), nsmall = 3)))]

    ## Create the combined effect string with expression for italic p
    effect_abbrev <- if(effect_label == "Odds Ratio") "aOR" else if(effect_label == "Risk Ratio") "aRR" else "Coef"

    to_show_exp_clean[, effect_string_expr := data.table::fifelse(
                                                              is.na(N) & is.na(estimate),
                                                              "''",
                                                              fcase(
                                                                  is.na(estimate), paste0("'", ref_label, "'"),
                                                                  
                                                                  p_value < 0.001 & !exponentiate & (conf_low < 0 | conf_high < 0),
                                                                  paste0("'", effect_formatted, " (", conf_low_formatted, ", ", 
                                                                         conf_high_formatted, "); '*~italic(p)~'< 0.001'"),
                                                                  
                                                                  p_value < 0.001,
                                                                  paste0("'", effect_formatted, " (", conf_low_formatted, "-", 
                                                                         conf_high_formatted, "); '*~italic(p)~'< 0.001'"),
                                                                  
                                                                  !exponentiate & (conf_low < 0 | conf_high < 0),
                                                                  paste0("'", effect_formatted, " (", conf_low_formatted, ", ", 
                                                                         conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'"),
                                                                  
                                                                  default = paste0("'", effect_formatted, " (", conf_low_formatted, "-", 
                                                                                   conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'")
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
    if(exponentiate) {
        to_show_exp_clean[is.na(estimate), estimate := 0]
    } else {
        to_show_exp_clean[is.na(estimate), estimate := 0]
    }

    ## Reorder (flip) - but maintain the variable grouping
    to_show_exp_clean <- to_show_exp_clean[order(nrow(to_show_exp_clean):1)]
    to_show_exp_clean[, x_pos := .I]

    ## Calculate plot ranges with better handling of extreme cases
    if(exponentiate) {
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
        
    } else {
        
        rangeb <- range(to_show$conf_low, to_show$conf_high, na.rm = TRUE)
        breaks <- pretty(rangeb, n = 7)
        breaks <- breaks[breaks >= rangeb[1] & breaks <= rangeb[2]]
        
        is_one_sided <- (min(rangeb) > 0) || (max(rangeb) < 0)
        
        if(!0 %in% breaks) {
            breaks <- sort(unique(c(breaks, 0)))
        }
        
        if(min(rangeb) > 0) {
            rangeb[1] <- -0.1 * abs(max(rangeb))
        } else if(max(rangeb) < 0) {
            rangeb[2] <- 0.1 * abs(min(rangeb))
        }
        
        reference_value <- 0
    }

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
    y_or <- layout$positions$effect

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

    ## Format deviance values
    null_dev_formatted <- format(round(gmodel$null_deviance, 2), nsmall = 2)
    resid_dev_formatted <- format(round(gmodel$residual_deviance, 2), nsmall = 2)
    pseudo_r2_formatted <- format(round(gmodel$pseudo_r2, 3), nsmall = 3)

    ## Create the plot
    if(exponentiate) {
        p <- ggplot2::ggplot(to_show_exp_clean, ggplot2::aes(x_pos, exp(estimate))) +
            
            ggplot2::geom_rect(ggplot2::aes(xmin = x_pos - .5, xmax = x_pos + .5,
                                            ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                                            fill = ordered(shade_group))) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
            ggplot2::scale_size_continuous(range = c(1, 6), guide = "none") +
            ggplot2::scale_fill_manual(values = shade_colors, guide = "none") +
            
            ggplot2::geom_point(ggplot2::aes(size = N), pch = 22, color = "#000000", fill = color, na.rm = TRUE) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = exp(conf_low), ymax = exp(conf_high)), width = 0.15) +
            
            ggplot2::annotate(geom = "segment",
                              x = -0.5, xend = -0.5,
                              y = exp(rangeb[1]),
                              yend = exp(rangeb[2]),
                              color = "#000000", linewidth = 1) +
            
            ggplot2::annotate(geom = "segment", 
                              x = -0.5, xend = max(to_show_exp_clean$x_pos) + 0.5, 
                              y = reference_value, yend = reference_value, 
                              linetype = "longdash") +
            
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
            
            ggplot2::coord_flip(ylim = exp(rangeplot)) +
            ggplot2::ggtitle(title) +
            ggplot2::scale_y_log10(name = effect_label,
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
                 ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = exp(y_variable),
                                   label = to_show_exp_clean$var_display, 
                                   fontface = data.table::fifelse(grepl("^    ", to_show_exp_clean$var_display), "plain", "bold"), 
                                   hjust = 0,
                                   size = annot_font)
             } else {
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
            ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.4, y = exp(y_or),
                              label = paste0("bold('", effect_abbrev, " (95% CI); '*bolditalic(p)*'-value')"),
                              hjust = 0, size = header_font, parse = TRUE) +
            
            ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = exp(y_or),
                              label = to_show_exp_clean$effect_string_expr, hjust = 0,
                              size = annot_font, parse = TRUE) +
            
            ## X-axis label
            ggplot2::annotate(geom = "text", x = -1.5, y = reference_value,
                              label = effect_label, fontface = "bold",
                              hjust = 0.5, vjust = 2, size = annot_font * 1.5) +
            
            ## Model statistics footer
            ggplot2::annotate(geom = "text", x = 0.5, y = exp(y_variable),
                              label = paste0("Observations analyzed: ", gmodel$nobs_with_pct,
                                             "\nModel: ", gmodel$family,
                                             "\nNull (Residual) Deviance: ", null_dev_formatted, " (", resid_dev_formatted, ")",
                                             "\nPseudo R^2: ", pseudo_r2_formatted,
                                             "\nAIC: ", gmodel$AIC_formatted),
                              size = annot_font * 0.8, hjust = 0, vjust = 1.2, fontface = "italic")
        
    } else {
        ## Non-exponentiated plot (linear scale)
        p <- ggplot2::ggplot(to_show_exp_clean, ggplot2::aes(x_pos, estimate)) +
            
            ggplot2::geom_rect(ggplot2::aes(xmin = x_pos - .5, xmax = x_pos + .5,
                                            ymin = rangeplot[1], ymax = rangeplot[2],
                                            fill = ordered(shade_group + 1))) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
            ggplot2::scale_size_continuous(range = c(1, 6), guide = "none") +
            ggplot2::scale_fill_manual(values = c("#FFFFFF", "#EEEEEE"), guide = "none") +
            
            ggplot2::geom_point(ggplot2::aes(size = N), pch = 22, color = "#000000", fill = color, na.rm = TRUE) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = conf_low, ymax = conf_high), width = 0.15) +
            
            ggplot2::annotate(geom = "segment",
                              x = -0.5, xend = -0.5,
                              y = rangeb[1],
                              yend = rangeb[2],
                              color = "#000000", linewidth = 1) +
            
            ggplot2::annotate(geom = "segment", 
                              x = -0.5, xend = max(to_show_exp_clean$x_pos) + 0.5, 
                              y = reference_value, yend = reference_value, 
                              linetype = "longdash") +
            
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
            
            ggplot2::coord_flip(ylim = rangeplot) +
            ggplot2::ggtitle(title) +
            ggplot2::scale_y_continuous(name = effect_label,
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
            ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = y_variable,
                              label = "Variable", fontface = "bold", hjust = 0,
                              size = header_font) +
            
            {if (indent_groups || condense_table) {
                 ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_variable,
                                   label = to_show_exp_clean$var_display, 
                                   fontface = data.table::fifelse(grepl("^    ", to_show_exp_clean$var_display), "plain", "bold"), 
                                   hjust = 0,
                                   size = annot_font)
             } else {
                 ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_variable,
                                   label = to_show_exp_clean$var_display, fontface = "bold", hjust = 0,
                                   size = annot_font)
             }} +
            
            ## Group/level column
            {if (!(indent_groups || condense_table)) {
                 list(
                     ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = y_level,
                                       label = "Group", fontface = "bold", hjust = 0,
                                       size = header_font),
                     ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_level,
                                       label = to_show_exp_clean$level, hjust = 0,
                                       size = annot_font)
                 )
             }} +
            
            ## N column (conditional)
            {if (show_n) {
                 list(
                     ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = y_n,
                                       label = "n", fontface = "bold.italic", hjust = 0.5,
                                       size = header_font),
                     ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_n,
                                       label = to_show_exp_clean$n_formatted, hjust = 0.5,
                                       size = annot_font)
                 )
             }} +
            
            ## Events column (conditional)
            {if (show_events) {
                 list(
                     ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = y_events,
                                       label = "Events", fontface = "bold", hjust = 0.5,
                                       size = header_font),
                     ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_events,
                                       label = to_show_exp_clean$events_formatted, hjust = 0.5,
                                       size = annot_font)
                 )
             }} +
            
            ## Effect column
            ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.4, y = y_or,
                              label = paste0("bold('", effect_abbrev, " (95% CI); '*bolditalic(p)*'-value')"),
                              hjust = 0, size = header_font, parse = TRUE) +
            
            ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_or,
                              label = to_show_exp_clean$effect_string_expr, hjust = 0,
                              size = annot_font, parse = TRUE) +
            
            ## X-axis label
            ggplot2::annotate(geom = "text", x = -1.5, y = reference_value,
                              label = effect_label, fontface = "bold",
                              hjust = 0.5, vjust = 2, size = annot_font * 1.5) +
            
            ## Model statistics at bottom
            ggplot2::annotate(geom = "text", x = 0.5, y = y_variable,
                              label = paste0("Observations analyzed: ", gmodel$nobs_with_pct,
                                             "\nModel: ", gmodel$family,
                                             "\nNull (Residual) Deviance: ", null_dev_formatted, " (", resid_dev_formatted, ")",
                                             "\nPseudo R^2: ", pseudo_r2_formatted,
                                             "\nAIC: ", gmodel$AIC_formatted),
                              size = annot_font * 0.8, hjust = 0, vjust = 1.2, fontface = "italic")
    }

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
