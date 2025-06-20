# Create LaTeX table for publication with oracle coverage and width
create_latex_table <- function(summary_table, caption = "Plateau Selector Results") {
  # Format numeric columns
  summary_table$Abs_Bias <- sprintf("%.6f", summary_table$Abs_Bias)
  summary_table$Std_Err <- sprintf("%.6f", summary_table$Std_Err)
  summary_table$MSE <- sprintf("%.6f", summary_table$MSE)
  summary_table$Coverage <- sprintf("%.2f", summary_table$Coverage)
  summary_table$Avg_CI_Width <- sprintf("%.4f", summary_table$Avg_CI_Width)
  summary_table$Oracle_Coverage <- sprintf("%.2f", summary_table$Oracle_Coverage)
  summary_table$Oracle_CI_Width <- sprintf("%.4f", summary_table$Oracle_CI_Width)

  # Create LaTeX table
  latex_table <- "\\begin{table}[ht]\n\\centering\n\\caption{"
  latex_table <- paste0(latex_table, caption, "}\n")
  latex_table <- paste0(latex_table, "\\begin{tabular}{lccccccccc}\n\\toprule\n")
  latex_table <- paste0(latex_table,
                        "n & Targeting & Selected j & Abs. Bias & Std. Error & MSE & Coverage (\\%) & CI Width & Oracle Cov. (\\%) & Oracle Width \\\\\n\\midrule\n")

  # Add rows grouped by sample size
  current_n <- -1
  for (i in 1:nrow(summary_table)) {
    row <- summary_table[i, ]

    # Add midrule between sample sizes
    if (current_n != row$n && current_n != -1) {
      latex_table <- paste0(latex_table, "\\midrule\n")
    }
    current_n <- row$n

    # Add data row
    latex_table <- paste0(latex_table,
                          row$n, " & ",
                          row$Targeting, " & ",
                          row$Selected_j, " & ",
                          row$Abs_Bias, " & ",
                          row$Std_Err, " & ",
                          row$MSE, " & ",
                          row$Coverage, " & ",
                          row$Avg_CI_Width, " & ",
                          row$Oracle_Coverage, " & ",
                          row$Oracle_CI_Width,
                          " \\\\\n")
  }

  # Close table
  latex_table <- paste0(latex_table, "\\bottomrule\n\\end{tabular}\n\\end{table}")

  return(latex_table)
}# Fixed Plateau Selector with NA Handling


# Basic plateau selector algorithm with NA handling
# plateauSelector <- function(data, method = c("relax", "proj", "delta")) {
#   method <- match.arg(method)
#
#   if (nrow(data) == 0) return(NULL)
#
#   # Extract sequences
#   estimates <- data[[paste0("psi_", method)]]
#   lowerBounds <- data[[paste0("lower_", method, "_np")]]
#   upperBounds <- data[[paste0("upper_", method, "_np")]]
#   jValues <- data$j
#
#   # Handle NAs by removing corresponding elements
#   valid_indices <- which(!is.na(estimates) & !is.na(lowerBounds) & !is.na(upperBounds))
#   if (length(valid_indices) == 0) return(NULL) # All values are NA
#
#   estimates <- estimates[valid_indices]
#   lowerBounds <- lowerBounds[valid_indices]
#   upperBounds <- upperBounds[valid_indices]
#   jValues <- jValues[valid_indices]
#
#   # If we have only one valid data point, return that j
#   if (length(jValues) == 1) return(jValues[1])
#
#   # Determine if estimates are predominantly increasing or decreasing
#   increases <- sum(diff(estimates) > 0, na.rm = TRUE)
#   decreases <- sum(diff(estimates) < 0, na.rm = TRUE)
#   isIncreasing <- increases >= decreases
#
#   # Find optimal j
#   optimalJ <- jValues[1]  # Default to first j
#
#   if (isIncreasing) {
#     # For increasing estimates, find where lower bound starts decreasing
#     for (i in 2:length(jValues)) {
#       if (!is.na(lowerBounds[i]) && !is.na(lowerBounds[i-1]) &&
#           lowerBounds[i] < lowerBounds[i-1]) {
#         optimalJ <- jValues[i-1]
#         break
#       }
#
#       # If we reach the end without finding a decrease, use the last j
#       if (i == length(jValues)) {
#         optimalJ <- jValues[i]
#       }
#     }
#   } else {
#     # For decreasing estimates, find where upper bound starts increasing
#     for (i in 2:length(jValues)) {
#       if (!is.na(upperBounds[i]) && !is.na(upperBounds[i-1]) &&
#           upperBounds[i] > upperBounds[i-1]) {
#         optimalJ <- jValues[i-1]
#         break
#       }
#
#       # If we reach the end without finding an increase, use the last j
#       if (i == length(jValues)) {
#         optimalJ <- jValues[i]
#       }
#     }
#   }
#
#   return(optimalJ)
# }

# Improved plateau selector with NA handling
plateauSelector <- function(data, method = c("relax", "proj", "delta")) {
  method <- match.arg(method)
  if (nrow(data) == 0) return(NULL)

  # Extract sequences
  estimates <- data[[paste0("psi_", method)]]
  lowerBounds <- data[[paste0("lower_", method, "_np")]]
  upperBounds <- data[[paste0("upper_", method, "_np")]]
  jValues <- data$j

  # Handle NAs by removing corresponding elements
  valid_indices <- which(!is.na(estimates) & !is.na(lowerBounds) & !is.na(upperBounds))
  if (length(valid_indices) == 0) return(NULL) # All values are NA

  estimates <- estimates[valid_indices]
  lowerBounds <- lowerBounds[valid_indices]
  upperBounds <- upperBounds[valid_indices]
  jValues <- jValues[valid_indices]

  # If we have only one valid data point, return that j
  if (length(jValues) == 1) return(jValues[1])

  # Check each pair of adjacent points
  for (i in 2:length(jValues)) {
    # Check if estimate is increasing or decreasing at this step
    if (estimates[i] > estimates[i-1]) {
      # Estimate is increasing, lower bound should be increasing
      if (lowerBounds[i] < lowerBounds[i-1]) {
        # Rule broken: lower bound decreased when estimate increased
        return(jValues[i-1])
      }
    } else if (estimates[i] < estimates[i-1]) {
      # Estimate is decreasing, upper bound should be decreasing
      if (upperBounds[i] > upperBounds[i-1]) {
        # Rule broken: upper bound increased when estimate decreased
        return(jValues[i-1])
      }
    }
    # If estimates are exactly equal, we don't check bounds as there's no clear direction
  }

  # If we made it through all comparisons without breaking the rule
  # return the last j value
  return(jValues[length(jValues)])
}

# Enhanced plateau selector with NA handling
enhancedPlateauSelector <- function(data, method = c("relax", "proj", "delta"), truthValue = 1.5) {
  method <- match.arg(method)

  if (nrow(data) == 0) return(NULL)

  # Extract sequences
  estimates <- data[[paste0("psi_", method)]]
  lowerBounds <- data[[paste0("lower_", method, "_np")]]
  upperBounds <- data[[paste0("upper_", method, "_np")]]
  jValues <- data$j

  # Handle NAs by removing corresponding elements
  valid_indices <- which(!is.na(estimates) & !is.na(lowerBounds) & !is.na(upperBounds))
  if (length(valid_indices) == 0) return(NULL) # All values are NA

  estimates <- estimates[valid_indices]
  lowerBounds <- lowerBounds[valid_indices]
  upperBounds <- upperBounds[valid_indices]
  jValues <- jValues[valid_indices]

  # If we have only one valid data point, return that j
  if (length(jValues) == 1) return(jValues[1])

  # Calculate CI widths
  ciWidths <- upperBounds - lowerBounds

  # Calculate relative changes in estimates and CI widths
  estChanges <- numeric(length(estimates) - 1)
  widthChanges <- numeric(length(ciWidths) - 1)

  for (i in 2:length(estimates)) {
    # Relative change in estimate
    estChanges[i-1] <- abs((estimates[i] - estimates[i-1]) / estimates[i-1])
    # Relative change in CI width
    widthChanges[i-1] <- (ciWidths[i] - ciWidths[i-1]) / ciWidths[i-1]
  }

  # Determine if estimates are predominantly increasing or decreasing
  increases <- sum(diff(estimates) > 0, na.rm = TRUE)
  decreases <- sum(diff(estimates) < 0, na.rm = TRUE)
  isIncreasing <- increases >= decreases

  # 1. First criterion: Find where estimates stabilize (<0.5% change) but CI width grows (>1%)
  optimalJ <- NULL

  for (i in 1:length(estChanges)) {
    if (!is.na(estChanges[i]) && !is.na(widthChanges[i]) &&
        estChanges[i] < 0.005 && widthChanges[i] > 0.01) {
      optimalJ <- jValues[i]
      break
    }
  }

  # 2. If no match, use standard criterion from basic selector
  if (is.null(optimalJ)) {
    if (isIncreasing) {
      # For increasing estimates, find where lower bound starts decreasing
      for (i in 2:length(jValues)) {
        if (!is.na(lowerBounds[i]) && !is.na(lowerBounds[i-1]) &&
            lowerBounds[i] < lowerBounds[i-1]) {
          optimalJ <- jValues[i-1]
          break
        }

        # If we reach the end without finding a decrease, use the last j
        if (i == length(jValues)) {
          optimalJ <- jValues[i]
        }
      }
    } else {
      # For decreasing estimates, find where upper bound starts increasing
      for (i in 2:length(jValues)) {
        if (!is.na(upperBounds[i]) && !is.na(upperBounds[i-1]) &&
            upperBounds[i] > upperBounds[i-1]) {
          optimalJ <- jValues[i-1]
          break
        }

        # If we reach the end without finding an increase, use the last j
        if (i == length(jValues)) {
          optimalJ <- jValues[i]
        }
      }
    }
  }

  return(optimalJ)
}

# Apply plateau selector to all data subsets with NA handling
applyPlateauSelector <- function(data, uniqueN, uniqueB, useEnhanced = FALSE, truthValue = 1.5) {
  results <- list()

  for (n in uniqueN) {
    results[[as.character(n)]] <- list()

    for (b in uniqueB) {
      # Get data subset for this n and B
      subset <- data[data$n == n & data$B == b, ]

      # Skip if no data
      if (nrow(subset) == 0) next

      # Sort by j
      subset <- subset[order(subset$j), ]

      # Apply selector for each method
      if (useEnhanced) {
        results[[as.character(n)]][[as.character(b)]] <- list(
          relax = tryCatch(enhancedPlateauSelector(subset, "relax", truthValue), error = function(e) NULL),
          proj = tryCatch(enhancedPlateauSelector(subset, "proj", truthValue), error = function(e) NULL),
          delta = tryCatch(enhancedPlateauSelector(subset, "delta", truthValue), error = function(e) NULL)
        )
      } else {
        results[[as.character(n)]][[as.character(b)]] <- list(
          relax = tryCatch(plateauSelector(subset, "relax"), error = function(e) NULL),
          proj = tryCatch(plateauSelector(subset, "proj"), error = function(e) NULL),
          delta = tryCatch(plateauSelector(subset, "delta"), error = function(e) NULL)
        )
      }
    }
  }

  return(results)
}

# Create summary table based on plateau selections with NA handling
createSummaryTable <- function(data, selections, truthValue) {
  # Initialize result table
  result <- data.frame()

  # Get unique sample sizes
  uniqueN <- sort(unique(data$n))

  for (n in uniqueN) {
    # For each method
    for (method in c("relax", "proj", "delta")) {
      # Get method name for table
      methodName <- switch(method,
                           relax = "Relaxed",
                           proj = "Projection",
                           delta = "Delta-method")

      # Get all selected j values for this n and method
      selectedJs <- list()
      for (b in names(selections[[as.character(n)]])) {
        j <- selections[[as.character(n)]][[b]][[method]]
        if (!is.null(j) && !is.na(j)) {
          selectedJs <- c(selectedJs, j)
        }
      }

      if (length(selectedJs) == 0) next

      # Find most common j (mode) - this is just for reporting, not for statistics
      jTable <- table(unlist(selectedJs))
      modeJ <- as.numeric(names(jTable)[which.max(jTable)])

      # Get rows with selected j values for each B
      selectedRows <- data.frame()
      for (b in names(selections[[as.character(n)]])) {
        j <- selections[[as.character(n)]][[b]][[method]]
        if (!is.null(j) && !is.na(j)) {
          row <- data[data$n == n & data$B == as.numeric(b) & data$j == j, ]
          if (nrow(row) > 0) {
            selectedRows <- rbind(selectedRows, row)
          }
        }
      }

      # Calculate statistics
      if (nrow(selectedRows) > 0) {
        estimates <- selectedRows[[paste0("psi_", method)]]
        # Handle NAs in estimates
        valid_estimates <- estimates[!is.na(estimates)]
        if (length(valid_estimates) == 0) next

        mean_est <- mean(valid_estimates)
        bias <- mean_est - truthValue
        abs_bias <- abs(bias)

        # Standard error
        std_err <- sd(valid_estimates)

        # MSE
        mse <- mean((valid_estimates - truthValue)^2)

        # Coverage using nonparametric CI
        lower_bound <- selectedRows[[paste0("lower_", method, "_np")]]
        upper_bound <- selectedRows[[paste0("upper_", method, "_np")]]
        valid_ci <- !is.na(lower_bound) & !is.na(upper_bound)
        if (sum(valid_ci) == 0) next

        coverage <- mean(truthValue >= lower_bound[valid_ci] & truthValue <= upper_bound[valid_ci], na.rm = TRUE) * 100

        # Average CI width
        avg_ci_width <- mean(upper_bound[valid_ci] - lower_bound[valid_ci], na.rm = TRUE)

        # Oracle coverage (based on empirical standard error of plateau-selected estimates)
        z_quantile <- 1.96 # For 95% coverage
        oracle_lower <- valid_estimates - z_quantile * std_err
        oracle_upper <- valid_estimates + z_quantile * std_err
        oracle_coverage <- mean(truthValue >= oracle_lower & truthValue <= oracle_upper, na.rm = TRUE) * 100

        # Oracle CI width
        oracle_ci_width <- mean(oracle_upper - oracle_lower, na.rm = TRUE)

        # Add to result table
        row <- data.frame(
          n = n,
          Targeting = methodName,
          Selected_j = modeJ,
          Abs_Bias = abs_bias,
          Std_Err = std_err,
          MSE = mse,
          Coverage = coverage,
          Avg_CI_Width = avg_ci_width,
          Oracle_Coverage = oracle_coverage,
          Oracle_CI_Width = oracle_ci_width
        )

        result <- rbind(result, row)
      }
    }
  }

  return(result)
}

# Main function to run the plateau selector analysis with NA handling
runPlateauAnalysis <- function(data, truthValue = 1.5, useEnhanced = FALSE, estimator_name = "A-TMLE") {
  # Get unique sample sizes and simulation rounds
  uniqueN <- unique(data$n)
  uniqueB <- unique(data$B)

  # Apply plateau selector to all data subsets
  selections <- applyPlateauSelector(data, uniqueN, uniqueB, useEnhanced, truthValue)

  # Create summary table
  summaryTable <- createSummaryTable(data, selections, truthValue)

  # Create LaTeX tables - both standard and formatted styles
  standard_latex <- create_latex_table(summaryTable,
                                       "Plateau Selector Results")

  formatted_latex <- create_formatted_latex_table(summaryTable,
                                                  paste0("Simulation results for ", estimator_name,
                                                         " with plateau-selected working model."),
                                                  estimator_name)

  return(list(
    selections = selections,
    summaryTable = summaryTable,
    standard_latex = standard_latex,
    formatted_latex = formatted_latex
  ))
}
