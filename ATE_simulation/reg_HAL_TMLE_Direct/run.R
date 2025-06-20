#Note to only consider the undersmoothed TMLE
#We need to change line 29 in target_TMLE.R to be:
#lambda_seq <- lambda_seq[c(1,min(nlambda_max, length(lambda_seq)))]



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
#registerDoMC(cores = 8)
registerDoMC(cores = availableCores()-1)
source("dgp/sim_data.R")
set.seed(123)

run <- function(dgp_num) {
    B <- 500
    n_seq <- seq(500, 2000, 500)


  res_df <- map_dfr(n_seq, function(.n) {
    map_dfr(seq(B), function(.b) {
      print("n = " %+% .n %+% ", run: " %+% .b)
      data <- do.call(eval("sim_data_" %+% dgp_num), list(n = .n))
      W <- data[, grep("W", colnames(data))]
      A <- data$A
      Y <- data$Y
      folds <- make_folds(n = .n, V = 5)
      # estimate P(A=1|W)
      g1W <- learn_g(W = W,
                     A = A,
                     method = "glm",
                     folds = folds,
                     g_bounds = c(0.01, 0.99),
                     cross_fit_nuisance = TRUE)


      # estimate Q(A,W) = E(Y|A,W)
      QAW<- learn_Q(W = W,
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


      target_args <- list(W = W,
                          A = A,
                          Y = Y,
                          g1W = g1W,
                          QAW = QAW,
                          delta = rep(1, .n),
                          pseudo_outcome = QAW$pseudo_outcome,
                          pseudo_weights = QAW$pesudo_weights,
                          X = QAW$X,
                          basis_list = QAW$basis_list,
                          X_hal = QAW$X_hal,
                          fit = QAW$fit,
                          dx = 1e-4,
                          max_iter = 2000,
                          seq = TRUE,
                          nlambda_max = 10,
                          verbose = FALSE,
                          browse = FALSE)

      target_args_proj <- target_args; target_args_proj$grad_proj <- TRUE; target_args_proj$method <- "tmle direct"
      QAW_proj <- do.call(target_TMLE, target_args_proj)



      # point estimate and inference
      res_df <- map_dfr(seq(length(QAW_proj)), function(.j) {

        .proj <- QAW_proj[[.j]]

        psi_proj <-  .proj$pred

        #Using the non parametric eic

        eic_proj <- .proj$eic

        se_proj_np <- sqrt(var(eic_proj, na.rm = TRUE) / .n)

        lower_proj_np <- psi_proj - 1.96 * se_proj_np
        upper_proj_np <- psi_proj + 1.96 * se_proj_np

        #Using the projected eic (weak)
        eic_proj <- .proj$eic_proj

        se_proj_proj <- sqrt(var(eic_proj, na.rm = TRUE) / .n)

        lower_proj_proj <- psi_proj - 1.96 * se_proj_proj
        upper_proj_proj <- psi_proj + 1.96 * se_proj_proj

        #Using the projected eic (cv)
        eic_proj <- .proj$eic_proj_cv

        se_proj_proj_cv <- sqrt(var(eic_proj, na.rm = TRUE) / .n)

        lower_proj_proj_cv <- psi_proj - 1.96 * se_proj_proj_cv
        upper_proj_proj_cv <- psi_proj + 1.96 * se_proj_proj_cv

        #Using the delta eic
        eic_proj <- .proj$eic_delta

        se_proj_delta <- sqrt(var(eic_proj, na.rm = TRUE) / .n)

        lower_proj_delta <- psi_proj - 1.96 * se_proj_delta
        upper_proj_delta <- psi_proj + 1.96 * se_proj_delta

        return(data.frame(n = .n,
                          B = .b,
                          j = .j,
                          proj_lambda = as.numeric(.proj$lambda[1]),
                          psi_proj = psi_proj,
                          ###non parametric eic based inference
                          se_proj_np = se_proj_np,
                          lower_proj_np = lower_proj_np,
                          upper_proj_np = upper_proj_np,
                          ###projected eic based inference (weak)

                          se_proj_proj = se_proj_proj,


                          lower_proj_proj = lower_proj_proj,
                          upper_proj_proj = upper_proj_proj,

                          ###projected eic based inference (cv)

                          se_proj_proj_cv = se_proj_proj_cv,


                          lower_proj_proj_cv = lower_proj_proj_cv,
                          upper_proj_proj_cv = upper_proj_proj_cv,

                          ###Delta eic based inference

                          se_proj_delta = se_proj_delta,


                          lower_proj_delta = lower_proj_delta,
                          upper_proj_delta = upper_proj_delta
                          ))
      })
      return(res_df)
    })
  })

  return(res_df)
}
