# Enhanced analysis function to create table matching the desired format
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

# Function to create formatted table similar to the example
create_tmle_table <- function(results_df, truth, estimator_name = "HAL-TMLE") {

  # Compute summary statistics
  summary_stats <- results_df %>%
    group_by(n) %>%
    summarize(
      # Basic performance metrics
      abs_bias = abs(mean(psi_tmle - truth)),
      std_err = sd(psi_tmle),
      mse = mean((psi_tmle - truth)^2),

      # EIC-based coverage (Non-parametric)
      cov_np_pct = mean(truth >= ci_lower & truth <= ci_upper) * 100,
      cov_np_width = mean(ci_upper - ci_lower),

      # Oracle coverage
      oracle_cov_pct = mean(truth >= psi_tmle - 1.96*sd(psi_tmle) &
                              truth <= psi_tmle + 1.96*sd(psi_tmle)) * 100,
      oracle_cov_width = 2 * 1.96 * sd(psi_tmle),

      .groups = 'drop'
    ) %>%
    mutate(
      # Format coverage columns as "percentage (width)"
      cov_np = sprintf("%.2f (%.3f)", cov_np_pct, cov_np_width),
      oracle_cov = sprintf("%.2f (%.3f)", oracle_cov_pct, oracle_cov_width),

      # Set projection and delta coverage to NA/empty
      cov_proj = NA_character_,
      cov_proj_cv = NA_character_,
      cov_delta = NA_character_,

      # Add estimator and targeting columns
      estimator = estimator_name,
      targeting = "Standard" # You can modify this based on your specific setup
    ) %>%
    select(estimator, n, targeting, abs_bias, std_err, mse,
           cov_np, cov_proj, cov_proj_cv, cov_delta, oracle_cov)

  return(summary_stats)
}

# Create LaTeX table - manually constructed to match desired format
create_latex_table <- function(table_data, caption = "Simulation results for HAL-TMLE") {

  # Round numeric columns for display
  table_formatted <- table_data %>%
    mutate(
      abs_bias = sprintf("%.4f", abs_bias),
      std_err = sprintf("%.4f", std_err),
      mse = sprintf("%.6f", mse)
    )

  n_rows <- nrow(table_formatted)
  estimator_name <- table_formatted$estimator[1]

  # Start building the LaTeX manually
  latex_lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\centering",
    "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{",
    "\\fontsize{8}{10}\\selectfont",
    "\\begin{tabular}[t]{cccrrrrrrrr}",
    "\\toprule",
    "Estimator & n & Targeting & Abs\\_Bias & Std\\_Err & MSE & Cov NP (\\%) & Cov Proj (\\%) & Cov Proj CV (\\%) & Cov Delta (\\%) & Oracle Cov (\\%)\\\\",
    "\\midrule"
  )

  # Add data rows
  for (i in 1:n_rows) {
    row_data <- table_formatted[i, ]

    if (i < n_rows) {
      # Empty estimator cell for first n-1 rows
      estimator_cell <- ""
    } else {
      # Multirow cell for last row
      estimator_cell <- paste0("\\multirow{-", n_rows, "}{*}{\\centering\\arraybackslash ", estimator_name, "}")
    }

    data_row <- paste(
      estimator_cell,
      row_data$n,
      row_data$targeting,
      row_data$abs_bias,
      row_data$std_err,
      row_data$mse,
      row_data$cov_np,
      ifelse(is.na(row_data$cov_proj), "NA", row_data$cov_proj),
      ifelse(is.na(row_data$cov_proj_cv), "NA", row_data$cov_proj_cv),
      ifelse(is.na(row_data$cov_delta), "NA", row_data$cov_delta),
      row_data$oracle_cov,
      sep = " & "
    )

    latex_lines <- c(latex_lines, paste0(data_row, "\\\\"))
  }

  # Close the table
  latex_lines <- c(
    latex_lines,
    "\\bottomrule",
    "\\end{tabular}}",
    "\\end{table}"
  )

  return(paste(latex_lines, collapse = "\n"))
}

# Alternative: If you want to match the exact multi-row format from your example
create_multirow_table <- function(results_df, truth, estimator_name = "HAL-TMLE") {

  # If you have multiple targeting methods, you would expand this
  # For now, creating a single method but showing how to structure for multiple

  summary_stats <- results_df %>%
    group_by(n) %>%
    summarize(
      abs_bias = abs(mean(psi_tmle - truth)),
      std_err = sd(psi_tmle),
      mse = mean((psi_tmle - truth)^2),
      cov_np_pct = mean(truth >= ci_lower & truth <= ci_upper) * 100,
      cov_np_width = mean(ci_upper - ci_lower),
      oracle_cov_pct = mean(truth >= psi_tmle - 1.96*sd(psi_tmle) &
                              truth <= psi_tmle + 1.96*sd(psi_tmle)) * 100,
      oracle_cov_width = 2 * 1.96 * sd(psi_tmle),
      .groups = 'drop'
    ) %>%
    # Create multiple rows per n if you have different targeting methods
    crossing(targeting = "Standard") %>%  # Replace with actual methods if you have them
    mutate(
      estimator = estimator_name,
      cov_np = sprintf("%.2f (%.3f)", cov_np_pct, cov_np_width),
      oracle_cov = sprintf("%.2f (%.3f)", oracle_cov_pct, oracle_cov_width),
      cov_proj = "",  # Empty string instead of NA for LaTeX
      cov_proj_cv = "",
      cov_delta = ""
    ) %>%
    arrange(n) %>%
    select(estimator, n, targeting, abs_bias, std_err, mse,
           cov_np, cov_proj, cov_proj_cv, cov_delta, oracle_cov)

  return(summary_stats)
}

# Usage with your data
source("dgp/sim_data_1.R")
truth <- get_truth_1()
res_df <- read.csv("out/tmle_linear_results_dgp1_parallel_0528.csv")

# Create the formatted table
formatted_table <- create_tmle_table(res_df, truth)
print(formatted_table)

# Generate LaTeX table with the desired format
latex_output <- create_latex_table(formatted_table,
                                   caption = "Simulation results for HAL-TMLE under DGP 1.")
cat(latex_output)

# Create multirow version (alternative format)
multirow_table <- create_multirow_table(res_df, truth)
print("Multi-row format:")
print(multirow_table)
