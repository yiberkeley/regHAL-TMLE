library(dplyr)
source("dgp/sim_data_1.R")
truth <- get_truth_1()

res_df1 <- read.csv("out/res_dgp_500_1000_Part1_1_0408_195103.csv")
res_df2 <- read.csv("out/res_dgp_500_1000_Part2_1_0408_195147.csv")
res_df3 <- read.csv("out/res_dgp_1500_Part1_1_0408_194632.csv")
res_df4 <- read.csv("out/res_dgp_1500_Part2_1_0408_194717.csv")
res_df5 <- read.csv("out/res_dgp_2000_Part1_1_0408_194907.csv")
res_df6 <- read.csv("out/res_dgp_2000_Part2_1_0408_194837.csv")

# res_df1 <- read.csv("out/res_dgp_500_1000_Part1_2_0408_195103.csv")
# res_df2 <- read.csv("out/res_dgp_500_1000_Part2_2_0408_195147.csv")
# res_df3 <- read.csv("out/res_dgp_1500_Part1_2_0408_194632.csv")
# res_df4 <- read.csv("out/res_dgp_1500_Part2_2_0408_194717.csv")
# res_df5 <- read.csv("out/res_dgp_2000_Part1_2_0408_194907.csv")
# res_df6 <- read.csv("out/res_dgp_2000_Part2_2_0408_194837.csv")

all_dfs<-list(res_df1,
              res_df2,
              res_df3,
              res_df4,
              res_df5,
              res_df6)

# Combine all dataframes into one
res_df <- do.call(rbind, all_dfs)

names(res_df)<-c(
  "n", "B", "j",
  "relax_lambda", "proj_lambda", "delta_lambda",
  "psi_relax", "psi_proj", "psi_delta",
  # Non-parametric EIC-based inference
  "se_relax_np", "se_proj_np", "se_delta_np",
  "lower_relax_np", "upper_relax_np",
  "lower_proj_np", "upper_proj_np",
  "lower_delta_np", "upper_delta_np",
  # Projected EIC-based inference(weak)
  "se_relax_proj", "se_proj_proj", "se_delta_proj",
  "lower_relax_proj", "upper_relax_proj",
  "lower_proj_proj", "upper_proj_proj",
  "lower_delta_proj", "upper_delta_proj",
  # Projected EIC-based inference(cv)
  "se_relax_proj_cv", "se_proj_proj_cv", "se_delta_proj_cv",
  "lower_relax_proj_cv", "upper_relax_proj_cv",
  "lower_proj_proj_cv", "upper_proj_proj_cv",
  "lower_delta_proj_cv", "upper_delta_proj_cv",
  # Delta EIC-based inference
  "se_relax_delta", "se_proj_delta", "se_delta_delta",
  "lower_relax_delta", "upper_relax_delta",
  "lower_proj_delta", "upper_proj_delta",
  "lower_delta_delta", "upper_delta_delta"
)

results <- runPlateauAnalysis(res_df, truthValue = truth, useEnhanced = FALSE, estimator_name = "A-TMLE")

# View the standard LaTeX table
cat(results$standard_latex)

# View the formatted LaTeX table
cat(results$formatted_latex)
