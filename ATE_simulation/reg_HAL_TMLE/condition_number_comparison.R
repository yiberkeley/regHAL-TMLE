# .libPaths(c("/global/home/users/skyqiu/R/x86_64-pc-linux-gnu-library/4.2",
#             .libPaths()))
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
#registerDoMC(cores = 3)
registerDoMC(cores = availableCores()-1)
source("dgp/sim_data.R")
set.seed(123)

run <- function(dgp_num) {
  B <- 200
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
                                                smoothness_orders = 1),
                    browse = FALSE)

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
                          seq = FALSE,
                          nlambda_max = 1,
                          verbose = FALSE,
                          browse = FALSE)
      target_args_delta <- target_args; target_args_delta$grad_proj <- FALSE; target_args_delta$method <- "weak reg tmle"
      QAW_delta <- do.call(target_TMLE, target_args_delta)

      # point estimate and inference
      res_df <- map_dfr(seq(length(QAW_delta)), function(.j) {
        .delta <- QAW_delta[[.j]]
        psi_delta <- mean(.delta$pred)
        eic_delta <- .delta$eic_delta
        kappa_IM <- .delta$eic_delta_kappa
        se_delta <- sqrt(var(eic_delta, na.rm = TRUE) / .n)
        lower_delta <- psi_delta - 1.96 * se_delta
        upper_delta <- psi_delta + 1.96 * se_delta

        return(data.frame(n = .n,
                          B = .b,
                          j = .j,
                          delta_lambda = .delta$lambda,
                          kappa_IM = kappa_IM,
                          psi_delta = psi_delta,
                          se_delta = se_delta,
                          lower_delta = lower_delta,
                          upper_delta = upper_delta))
      })
      return(res_df)
    })
  })

  return(res_df)
}

res_df <- run(2)
write.csv(res_df, file = "out/res_kappa_dgp_" %+% 2 %+% "_" %+% timestamp %+% ".csv", row.names = FALSE)
