# Create LaTeX table in the requested format
create_formatted_latex_table <- function(summary_table, caption = "Plateau Selector Results", estimator_name = "A-TMLE") {
  # Format numeric columns
  summary_table$Abs_Bias <- sprintf("%.4f", summary_table$Abs_Bias)
  summary_table$Std_Err <- sprintf("%.4f", summary_table$Std_Err)
  summary_table$MSE <- sprintf("%.6f", summary_table$MSE)
  summary_table$Coverage <- sprintf("%.2f", summary_table$Coverage)
  summary_table$Avg_CI_Width <- sprintf("%.3f", summary_table$Avg_CI_Width)
  summary_table$Oracle_Coverage <- sprintf("%.2f", summary_table$Oracle_Coverage)
  summary_table$Oracle_CI_Width <- sprintf("%.3f", summary_table$Oracle_CI_Width)

  # Get unique sample sizes for separating with cmidrule
  unique_n <- unique(summary_table$n)
  unique_n <- sort(unique_n)

  # Start building the LaTeX table
  latex_table <- "\\begin{table}\n\\centering\n\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{\n"
  latex_table <- paste0(latex_table, "\\fontsize{8}{10}\\selectfont\n")
  latex_table <- paste0(latex_table, "\\begin{tabular}[t]{cccrrrrrr}\n\\toprule\n")

  # Headers
  latex_table <- paste0(latex_table,
                        "Estimator & n & Targeting & Abs\\_Bias & Std\\_Err & MSE & ",
                        "Cov NP (\\%) & Selected j & Oracle Cov (\\%) \\\\\n\\midrule\n")

  # Group rows by sample size
  for (i in 1:length(unique_n)) {
    n_value <- unique_n[i]
    n_rows <- summary_table[summary_table$n == n_value, ]

    # Sort by targeting method to ensure consistent order
    n_rows <- n_rows[order(n_rows$Targeting), ]

    # Add rows for this sample size
    for (j in 1:nrow(n_rows)) {
      row <- n_rows[j, ]

      # For first row, leave Estimator blank; for last row of all rows, add the estimator name with multirow
      if (i == 1 && j == 1) {
        estimator_cell <- " "
      } else if (i == length(unique_n) && j == nrow(n_rows)) {
        total_rows <- nrow(summary_table)
        estimator_cell <- paste0("\\multirow{-", total_rows, "}{*}{\\centering\\arraybackslash ", estimator_name, "}")
      } else {
        estimator_cell <- " "
      }

      # Format CI width as "Coverage (Width)"
      coverage_with_width <- paste0(row$Coverage, " (", row$Avg_CI_Width, ")")
      oracle_with_width <- paste0(row$Oracle_Coverage, " (", row$Oracle_CI_Width, ")")

      # Add data row
      latex_table <- paste0(latex_table,
                            estimator_cell, " & ",
                            row$n, " & ",
                            row$Targeting, " & ",
                            row$Abs_Bias, " & ",
                            row$Std_Err, " & ",
                            row$MSE, " & ",
                            coverage_with_width, " & ",
                            row$Selected_j, " & ",
                            oracle_with_width,
                            " \\\\\n")
    }

    # Add cmidrule between sample sizes, except after the last group
    if (i < length(unique_n)) {
      latex_table <- paste0(latex_table, "\\cmidrule(lr){2-9}\n")
    }
  }

  # Close table
  latex_table <- paste0(latex_table, "\\bottomrule\n\\end{tabular}}\n")
  latex_table <- paste0(latex_table, "\\caption{", caption, "}\n\\end{table}")

  return(latex_table)
}

