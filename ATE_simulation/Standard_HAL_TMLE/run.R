# Note: Using undersmoothed TMLE approach with LINEAR updates and propensity score truncation
# Modified lambda sequence to only use first and last values
# Targeting step uses linear regression instead of logistic regression
# Propensity scores are truncated to improve stability
# PARALLELIZED across simulation replications

.libPaths(c("/global/home/users/yili/R/x86_64-pc-linux-gnu-library/4.4",
            .libPaths()))

`%+%` <- function(a, b) paste0(a, b)
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")
library(purrr)
library(torch)
library(origami)
library(hal9001)
library(glmnet)
library(furrr)
library(doMC)
library(devtools)
load_all()

# Set up parallel processing
library(future)
# Option 1: Use multicore (recommended for Unix/Linux/Mac)
plan(multicore, workers = availableCores() - 1)
# Option 2: Alternative for Windows or if multicore doesn't work
# plan(multisession, workers = availableCores() - 1)

source("dgp/sim_data.R")
set.seed(123)

run <- function(dgp_num) {
  B <- 500
  n_seq <- seq(500, 2000, 500)

  res_df <- map_dfr(n_seq, function(.n) {
    # PARALLEL across replications using future_map_dfr
    future_map_dfr(seq(B), function(.b) {
      print("n = " %+% .n %+% ", run: " %+% .b)

      # Generate data
      data <- do.call(eval("sim_data_" %+% dgp_num), list(n = .n))
      W <- data[, grep("W", colnames(data))]
      A <- data$A
      Y <- data$Y
      folds <- make_folds(n = .n, V = 5)

      # Estimate P(A=1|W) - propensity score
      g1W <- learn_g(W = W,
                     A = A,
                     method = "glm",
                     folds = folds,
                     g_bounds = c(0.01, 0.99),
                     cross_fit_nuisance = TRUE)

      # Additional propensity score truncation for stability
      truncation_level <- 0.01  # Can be adjusted
      g1W_truncated <- pmax(pmin(g1W, 1 - truncation_level), truncation_level)

      print(paste0("Propensity score range before truncation: [",
                   round(min(g1W), 4), ", ", round(max(g1W), 4), "]"))
      print(paste0("Propensity score range after truncation: [",
                   round(min(g1W_truncated), 4), ", ", round(max(g1W_truncated), 4), "]"))

      # Estimate Q(A,W) = E(Y|A,W) - outcome regression
      QAW <- learn_Q(W = W,
                     A = A,
                     Y = Y,
                     delta = rep(1, .n),
                     v_folds = 5,
                     weights = rep(1, .n),
                     enumerate_basis_args = list(max_degree = 2,
                                                 smoothness_orders = 1,
                                                 num_knots = length(Y)/20),
                     browse = FALSE)

      print("QAW initial finished")

      # ========== TARGETING STEP (LINEAR UPDATE) ==========

      # Get initial Q estimates for A=1 and A=0 using HAL predictions
      # QAW is a list with design matrices and fitted model
      Q1W_init <- predict(QAW$fit, newx = QAW$X1_hal, s = "lambda.min")
      Q0W_init <- predict(QAW$fit, newx = QAW$X0_hal, s = "lambda.min")

      # Convert to vectors (predict returns matrix)
      Q1W_init <- as.vector(Q1W_init)
      Q0W_init <- as.vector(Q0W_init)
      QAW_init <- ifelse(A == 1, Q1W_init, Q0W_init)

      # Compute clever covariate using TRUNCATED propensity scores
      # For ATE: H_n(A,W) = A/g1W - (1-A)/(1-g1W)
      H1W <- 1 / g1W_truncated  # Clever covariate for A=1
      H0W <- -1 / (1 - g1W_truncated)  # Clever covariate for A=0
      HAW <- A * H1W + (1 - A) * H0W  # Observed clever covariate

      # LINEAR Targeting step: Fit linear regression with clever covariate
      # Q^* = Q_init + epsilon * H
      targeting_data <- data.frame(
        Y = Y,
        Q_init = QAW_init,
        H = HAW
      )

      # Fit linear targeting model: Y ~ Q_init + epsilon * H (without intercept for Q_init)
      targeting_fit <- lm(Y ~ Q_init + H - 1, data = targeting_data)

      # Extract epsilon (targeting parameter) and alpha (coefficient for Q_init)
      epsilon <- coef(targeting_fit)["H"]
      alpha <- coef(targeting_fit)["Q_init"]

      # Handle convergence issues
      if(is.na(epsilon)) epsilon <- 0
      if(is.na(alpha)) alpha <- 1

      # Update Q estimates using LINEAR update
      Q1W_updated <- alpha * Q1W_init + epsilon * H1W
      Q0W_updated <- alpha * Q0W_init + epsilon * H0W
      QAW_updated <- ifelse(A == 1, Q1W_updated, Q0W_updated)

      print("Linear targeting step completed")

      # ========== POINT ESTIMATION ==========

      # Compute TMLE estimate of ATE
      psi_tmle <- mean(Q1W_updated) - mean(Q0W_updated)

      # Also compute initial (untargeted) estimate for comparison
      psi_initial <- mean(Q1W_init) - mean(Q0W_init)

      # ========== EFFICIENT INFLUENCE CURVE AND INFERENCE ==========

      # Compute efficient influence curve for each observation
      # EIC_i = H1W_i * (Y_i - Q1W_updated_i) * A_i + H0W_i * (Y_i - Q0W_updated_i) * (1-A_i) +
      #         Q1W_updated_i - Q0W_updated_i - psi_tmle

      eic <- rep(NA, .n)
      for(i in 1:.n) {
        if(A[i] == 1) {
          eic_term1 <- H1W[i] * (Y[i] - Q1W_updated[i])
          eic_term2 <- 0
        } else {
          eic_term1 <- 0
          eic_term2 <- H0W[i] * (Y[i] - Q0W_updated[i])
        }
        eic_term3 <- Q1W_updated[i] - Q0W_updated[i] - psi_tmle
        eic[i] <- eic_term1 + eic_term2 + eic_term3
      }

      # Alternative vectorized EIC calculation
      eic_vec <- A * H1W * (Y - Q1W_updated) +
        (1 - A) * H0W * (Y - Q0W_updated) +
        (Q1W_updated - Q0W_updated) - psi_tmle

      # Standard error based on EIC
      se_tmle <- sqrt(var(eic_vec, na.rm = TRUE) / .n)

      # 95% Confidence interval
      ci_lower <- psi_tmle - 1.96 * se_tmle
      ci_upper <- psi_tmle + 1.96 * se_tmle

      # P-value for testing H0: psi = 0
      z_stat <- psi_tmle / se_tmle
      p_value <- 2 * (1 - pnorm(abs(z_stat)))

      # ========== DIAGNOSTICS ==========

      # Check EIC mean (should be close to 0)
      eic_mean <- mean(eic_vec, na.rm = TRUE)

      # Epsilon convergence check
      epsilon_converged <- abs(epsilon) < 1e-4

      # Range of propensity scores (before and after truncation)
      g1W_range_orig <- range(g1W)
      g1W_range_trunc <- range(g1W_truncated)

      # Proportion of observations affected by truncation
      prop_truncated <- mean(g1W != g1W_truncated)

      # Range of Q estimates
      Q_range_init <- range(c(Q1W_init, Q0W_init))
      Q_range_updated <- range(c(Q1W_updated, Q0W_updated))

      print("Inference completed")

      # ========== RETURN RESULTS ==========

      return(data.frame(
        n = .n,
        B = .b,
        # Point estimates
        psi_initial = psi_initial,
        psi_tmle = psi_tmle,
        # Targeting parameters (linear model)
        epsilon = epsilon,
        alpha = alpha,
        epsilon_converged = epsilon_converged,
        # EIC-based inference
        se_tmle = se_tmle,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        z_stat = z_stat,
        p_value = p_value,
        # Diagnostics
        eic_mean = eic_mean,
        # Propensity score info (original and truncated)
        g1W_min_orig = g1W_range_orig[1],
        g1W_max_orig = g1W_range_orig[2],
        g1W_min_trunc = g1W_range_trunc[1],
        g1W_max_trunc = g1W_range_trunc[2],
        prop_truncated = prop_truncated,
        truncation_level = truncation_level,
        # Q estimate ranges
        Q_range_init_min = Q_range_init[1],
        Q_range_init_max = Q_range_init[2],
        Q_range_updated_min = Q_range_updated[1],
        Q_range_updated_max = Q_range_updated[2]
      ))
    }, .options = furrr_options(seed = TRUE))  # Ensures reproducible random seeds across workers
  })

  return(res_df)
}

# Example usage:
#results_dgp1 <- run(dgp_num = 1)
results_dgp2 <- run(dgp_num = 2)

# Save results as CSV with information about linear updates and truncation
timestamp <- format(Sys.Date(), "%m%d")
#write.csv(results_dgp1, file = paste0("out/tmle_linear_results_dgp1_parallel_", timestamp, ".csv"), row.names = FALSE)
write.csv(results_dgp2, file = paste0("out/tmle_linear_results_dgp2_parallel_", timestamp, ".csv"), row.names = FALSE)

# Clean up parallel workers when done
plan(sequential)
