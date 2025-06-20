#Note to only consider the undersmoothed TMLE
#We need to change line 29 in target_TMLE.R to be:
#lambda_seq <- lambda_seq[1:min(nlambda_max, length(lambda_seq))]


.libPaths(c("/global/home/users/yili/R/x86_64-pc-linux-gnu-library/4.4",
            .libPaths()))

# Modify the beginning of your run.R script:
`%+%` <- function(a, b) paste0(a, b)
timestamp <- format(Sys.Date(), "%m%d") %+% "_" %+% format(Sys.time(), "%H%M%S")

library(purrr)
library(torch)
library(origami)
library(hal9001)
library(glmnet)
library(foreach)
library(doMC)
library(devtools)

# Load your package
load_all()

# Set up parallel backend with doMC
registerDoMC(cores = availableCores()-1)

source("dgp/sim_data.R")
#set.seed(123)
set.seed(234)

run <- function(dgp_num) {
  B <- 100
  n_seq <- c(500,1000)

  # Create an empty dataframe to store results
  results_all <- data.frame()

  for (n in n_seq) {
    # Use foreach for parallelization over the B iterations
    results_n <- foreach(b = 1:B, .combine = rbind) %dopar% {
      print(paste0("n = ", n, ", run: ", b, ", worker: ", Sys.getpid()))

      data <- do.call(eval(parse(text = paste0("sim_data_", dgp_num))), list(n = n))
      W <- data[, grep("W", colnames(data))]
      A <- data$A
      Y <- data$Y
      folds <- make_folds(n = n, V = 5)

      # Your existing code continues here...
      # estimate P(A=1|W)
      g1W <- learn_g(W = W,
                     A = A,
                     method = "glm",
                     folds = folds,
                     g_bounds = c(0.01, 0.99),
                     cross_fit_nuisance = TRUE)

      # estimate Q(A,W) = E(Y|A,W)
      QAW <- learn_Q(W = W,
                     A = A,
                     Y = Y,
                     delta = rep(1, n),
                     v_folds = 5,
                     weights = rep(1, n),
                     enumerate_basis_args = list(max_degree = 2,
                                                 smoothness_orders = 1),
                     browse = FALSE)

      print("QAW initial finished")

      # Rest of your code...
      # I'll keep this part as is since it looks correct
      target_args <- list(W = W,
                          A = A,
                          Y = Y,
                          g1W = g1W,
                          QAW = QAW,
                          delta = rep(1, n),
                          pseudo_outcome = QAW$pseudo_outcome,
                          pseudo_weights = QAW$pesudo_weights,
                          X = QAW$X,
                          basis_list = QAW$basis_list,
                          X_hal = QAW$X_hal,
                          fit = QAW$fit,
                          dx = 1e-4,
                          max_iter = 2000,
                          seq = TRUE,
                          nlambda_max = 20,
                          verbose = FALSE,
                          browse = FALSE)

      target_args_relax <- target_args
      target_args_relax$grad_proj <- TRUE
      target_args_relax$method <- "relax analytic"

      target_args_proj <- target_args
      target_args_proj$grad_proj <- TRUE
      target_args_proj$method <- "weak reg tmle"

      target_args_delta <- target_args
      target_args_delta$grad_proj <- FALSE
      target_args_delta$method <- "weak reg tmle"

      QAW_relax <- do.call(target_TMLE, target_args_relax)
      print("QAW relax finished")

      QAW_proj <- do.call(target_TMLE, target_args_proj)
      print("QAW proj finished")

      QAW_delta <- do.call(target_TMLE, target_args_delta)
      print("QAW delta finished")

      # point estimate and inference
      inner_res <- map_dfr(seq(length(QAW_proj)), function(.j) {
        .relax <- QAW_relax[[.j]]
        .proj <- QAW_proj[[.j]]
        .delta <- QAW_delta[[.j]]

        psi_relax <- .relax$pred
        psi_proj <-  .proj$pred
        psi_delta <- .delta$pred

        #Using the non parametric eic
        eic_relax <- .relax$eic
        eic_proj <- .proj$eic
        eic_delta <- .delta$eic

        se_relax_np <- sqrt(var(eic_relax, na.rm = TRUE) / n)
        se_proj_np <- sqrt(var(eic_proj, na.rm = TRUE) / n)
        se_delta_np<- sqrt(var(eic_delta, na.rm = TRUE) / n)

        lower_relax_np <- psi_relax - 1.96 * se_relax_np
        upper_relax_np <- psi_relax + 1.96 * se_relax_np
        lower_proj_np <- psi_proj - 1.96 * se_proj_np
        upper_proj_np <- psi_proj + 1.96 * se_proj_np
        lower_delta_np <- psi_delta - 1.96 * se_delta_np
        upper_delta_np <- psi_delta + 1.96 * se_delta_np

        #Using the projected eic (weak)
        eic_relax <- .relax$eic_proj
        eic_proj <- .proj$eic_proj
        eic_delta <- .delta$eic_proj

        se_relax_proj <- sqrt(var(eic_relax, na.rm = TRUE) / n)
        se_proj_proj <- sqrt(var(eic_proj, na.rm = TRUE) / n)
        se_delta_proj <- sqrt(var(eic_delta, na.rm = TRUE) / n)

        lower_relax_proj <- psi_relax - 1.96 * se_relax_proj
        upper_relax_proj <- psi_relax + 1.96 * se_relax_proj
        lower_proj_proj <- psi_proj - 1.96 * se_proj_proj
        upper_proj_proj <- psi_proj + 1.96 * se_proj_proj
        lower_delta_proj <- psi_delta - 1.96 * se_delta_proj
        upper_delta_proj <- psi_delta + 1.96 * se_delta_proj

        #Using the projected eic (cv)
        eic_relax <- .relax$eic_proj_cv
        eic_proj <- .proj$eic_proj_cv
        eic_delta <- .delta$eic_proj_cv

        se_relax_proj_cv <- sqrt(var(eic_relax, na.rm = TRUE) / n)
        se_proj_proj_cv <- sqrt(var(eic_proj, na.rm = TRUE) / n)
        se_delta_proj_cv <- sqrt(var(eic_delta, na.rm = TRUE) / n)

        lower_relax_proj_cv <- psi_relax - 1.96 * se_relax_proj_cv
        upper_relax_proj_cv <- psi_relax + 1.96 * se_relax_proj_cv
        lower_proj_proj_cv <- psi_proj - 1.96 * se_proj_proj_cv
        upper_proj_proj_cv <- psi_proj + 1.96 * se_proj_proj_cv
        lower_delta_proj_cv <- psi_delta - 1.96 * se_delta_proj_cv
        upper_delta_proj_cv <- psi_delta + 1.96 * se_delta_proj_cv

        #Using the delta eic
        eic_relax <- .relax$eic_delta
        eic_proj <- .proj$eic_delta
        eic_delta <- .delta$eic_delta

        se_relax_delta <- sqrt(var(eic_relax, na.rm = TRUE) / n)
        se_proj_delta <- sqrt(var(eic_proj, na.rm = TRUE) / n)
        se_delta_delta <- sqrt(var(eic_delta, na.rm = TRUE) / n)

        lower_relax_delta <- psi_relax - 1.96 * se_relax_delta
        upper_relax_delta <- psi_relax + 1.96 * se_relax_delta
        lower_proj_delta <- psi_proj - 1.96 * se_proj_delta
        upper_proj_delta <- psi_proj + 1.96 * se_proj_delta
        lower_delta_delta <- psi_delta - 1.96 * se_delta_delta
        upper_delta_delta <- psi_delta + 1.96 * se_delta_delta

        return(data.frame(n = n,
                          B = b,
                          j = .j,
                          relax_lambda = as.numeric(.relax$lambda[1]),
                          proj_lambda = as.numeric(.proj$lambda[1]),
                          delta_lambda = as.numeric(.delta$lambda[1]),
                          psi_relax = psi_relax,
                          psi_proj = psi_proj,
                          psi_delta = psi_delta,
                          ###non parametric eic based inference
                          se_relax_np = se_relax_np,
                          se_proj_np = se_proj_np,
                          se_delta_np = se_delta_np,
                          lower_relax_np = lower_relax_np,
                          upper_relax_np = upper_relax_np,
                          lower_proj_np = lower_proj_np,
                          upper_proj_np = upper_proj_np,
                          lower_delta_np = lower_delta_np,
                          upper_delta_np = upper_delta_np,
                          ###projected eic based inference (weak)
                          se_relax_proj = se_relax_proj,
                          se_proj_proj = se_proj_proj,
                          se_delta_proj = se_delta_proj,
                          lower_relax_proj = lower_relax_proj,
                          upper_relax_proj = upper_relax_proj,
                          lower_proj_proj = lower_proj_proj,
                          upper_proj_proj = upper_proj_proj,
                          lower_delta_proj = lower_delta_proj,
                          upper_delta_proj = upper_delta_proj,
                          ###projected eic based inference (cv)
                          se_relax_proj_cv = se_relax_proj_cv,
                          se_proj_proj_cv = se_proj_proj_cv,
                          se_delta_proj_cv = se_delta_proj_cv,
                          lower_relax_proj_cv = lower_relax_proj_cv,
                          upper_relax_proj_cv = upper_relax_proj_cv,
                          lower_proj_proj_cv = lower_proj_proj_cv,
                          upper_proj_proj_cv = upper_proj_proj_cv,
                          lower_delta_proj_cv = lower_delta_proj_cv,
                          upper_delta_proj_cv = upper_delta_proj_cv,
                          ###Delta eic based inference
                          se_relax_delta = se_relax_delta,
                          se_proj_delta = se_proj_delta,
                          se_delta_delta = se_delta_delta,
                          lower_relax_delta = lower_relax_delta,
                          upper_relax_delta = upper_relax_delta,
                          lower_proj_delta = lower_proj_delta,
                          upper_proj_delta = upper_proj_delta,
                          lower_delta_delta = lower_delta_delta,
                          upper_delta_delta = upper_delta_delta))
      })
      return(inner_res)
    }

    # Combine results for this n
    results_all <- rbind(results_all, results_n)
  }

  return(results_all)
}
