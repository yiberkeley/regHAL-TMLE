library(dplyr)
source("dgp/sim_data_2.R")
truth <- get_truth_2()

res_df <- read.csv("out/res_dgp_2_0613_140411.csv")

names(res_df)<-c(
  "n", "B", "j",
  "proj_lambda",
  "psi_proj",
  # Non-parametric EIC-based inference
  "se_proj_np",
  "lower_proj_np", "upper_proj_np",
  # Projected EIC-based inference(weak)
  "se_proj_proj",
  "lower_proj_proj", "upper_proj_proj",
  # Projected EIC-based inference(cv)
  "se_proj_proj_cv",
  "lower_proj_proj_cv", "upper_proj_proj_cv",
  # Delta EIC-based inference
  "se_proj_delta",
  "lower_proj_delta", "upper_proj_delta"
)

# Quick summary for projection method only
res_df %>% filter(j==2) %>%
  summarize(abs_bias_proj = abs(mean(psi_proj - truth)),
            se_proj = sd(psi_proj),
            mse_proj = mean((psi_proj - truth)^2),
            cover_proj_np = mean(truth >= lower_proj_np & truth <= upper_proj_np),
            cover_proj_proj = mean(truth >= lower_proj_proj & truth <= upper_proj_proj),
            cover_proj_proj_cv = mean(truth >= lower_proj_proj_cv & truth <= upper_proj_proj_cv),
            cover_proj_delta = mean(truth >= lower_proj_delta & truth <= upper_proj_delta),
            cover_oracle_proj = mean(truth >= psi_proj-1.96*sd(psi_proj) & truth <= psi_proj+1.96*sd(psi_proj)),
            .by = n)

###################################################################################################################

# Load required packages
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)

# Compute summary statistics for projection method only
temp <- res_df %>%
  filter(j == 2) %>%
  group_by(n) %>%
  summarize(
    abs_bias_proj  = abs(mean(psi_proj - truth)),
    se_proj        = sd(psi_proj),
    mse_proj       = mean((psi_proj - truth)^2),
    # Non-parametric EIC-based inference (coverage and CI length)
    cov_proj_np    = sprintf("%.2f (%.3f)",
                             mean(truth >= lower_proj_np & truth <= upper_proj_np) * 100,
                             mean(upper_proj_np - lower_proj_np)),
    # Projected EIC-based inference (weak)
    cov_proj_proj  = sprintf("%.2f (%.3f)",
                             mean(truth >= lower_proj_proj & truth <= upper_proj_proj) * 100,
                             mean(upper_proj_proj - lower_proj_proj)),
    # Projected EIC-based inference (cv)
    cov_proj_proj_cv  = sprintf("%.2f (%.3f)",
                                mean(truth >= lower_proj_proj_cv & truth <= upper_proj_proj_cv) * 100,
                                mean(upper_proj_proj_cv - lower_proj_proj_cv)),
    # Delta EIC-based inference
    cov_proj_delta  = sprintf("%.2f (%.3f)",
                              mean(truth >= lower_proj_delta & truth <= upper_proj_delta) * 100,
                              mean(upper_proj_delta - lower_proj_delta)),
    # Oracle coverage (CI length computed from 1.96*sd)
    cov_oracle_proj  = sprintf("%.2f (%.3f)",
                               mean(truth >= psi_proj - 1.96*sd(psi_proj) & truth <= psi_proj + 1.96*sd(psi_proj)) * 100,
                               2 * 1.96 * sd(psi_proj))
  )

# Create table for projection method
temp_table <- temp %>% transmute(
  n,
  Targeting = "Projection",
  Abs_Bias = sprintf("%.4f", abs_bias_proj),
  Std_Err = sprintf("%.4f", se_proj),
  MSE = sprintf("%.6f", mse_proj),
  `Cov NP (%)` = cov_proj_np,
  `Cov Proj (%)` = cov_proj_proj,
  `Cov Proj CV (%)` = cov_proj_proj_cv,
  `Cov Delta (%)` = cov_proj_delta,
  `Oracle Cov (%)` = cov_oracle_proj
)

# Add an Estimator column so that "A-TMLE" appears in every row
temp_table <- temp_table %>%
  mutate(Estimator = "A-TMLE") %>%
  select(Estimator, everything())

# Generate the LaTeX table using kable
latex_table <- kable(temp_table, format = "latex", booktabs = TRUE,
                     caption = "Simulation results for A-TMLE. DGP 1. Undersmoothed (Projection Method Only)",
                     align = c("c", "c", "c", "r", "r", "r", "r", "r", "r", "r", "r")) %>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle")

# Make the table font smaller using kable_styling()
latex_table <- latex_table %>%
  kable_styling(latex_options = c("scale_down"), font_size = 8)

# Print the LaTeX table code
cat(latex_table)

