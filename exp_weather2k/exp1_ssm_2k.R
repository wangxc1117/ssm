library(Rcpp)
Rcpp::sourceCpp("em_use.cpp")
library(Matrix)
source("function_fast.R")

run_one <- function(Z_all, gg_all, T, t0, m, m_t, k0, seed,
                    error = 1e-6, max_iter = 200,
                    cache_dir = "C:/Users/user/Desktop/ssm/cache") {

  T_all <- ncol(Z_all)
  stopifnot(t0 >= 1)
  stopifnot(t0 + T - 1 <= T_all)
  stopifnot(T_all >= T + 1)

  set.seed(seed)

  pickm  <- sample(seq_len(nrow(Z_all)), m)
  pool   <- setdiff(seq_len(nrow(Z_all)), pickm)
  pickn2 <- sample(pool, m_t)

  Z     <- Z_all[pickm,  t0:(t0 + T - 1), drop = FALSE]
  ytrue <- Z_all[pickn2, t0:(t0 + T - 1), drop = FALSE]

  gg    <- gg_all[pickm,  , drop = FALSE]
  grids <- gg_all[pickn2, , drop = FALSE]

  mu0 <- mean(Z, na.rm = TRUE)
  sd0 <- sd(as.vector(Z), na.rm = TRUE)
  if (!is.finite(sd0) || sd0 <= 0) sd0 <- 1

  Z_sc     <- (Z - mu0) / sd0
  ytrue_sc <- (ytrue - mu0) / sd0

  obs_idx <- lapply(seq_len(T), function(tt) which(!is.na(Z_sc[, tt])))
  z_list  <- lapply(seq_len(T), function(tt) {
    idxt <- obs_idx[[tt]]
    Z_sc[idxt, tt]
  })

  make_z_adj_list <- function(z_list, beta_t) {
    lapply(seq_along(z_list), function(tt) z_list[[tt]] - beta_t[tt])
  }

  make_d_list <- function(sig2_delta, sigma2_eps, obs_idx) {
    lapply(seq_along(obs_idx), function(tt) {
      rep(sigma2_eps + sig2_delta[tt], length(obs_idx[[tt]]))
    })
  }

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(
    cache_dir,
    sprintf("basis_w2k_m%d_k%d_T%d_seed%d_t0%d.rds", m, k0, T, seed, t0)
  )

  if (file.exists(cache_file)) {
    basis_obj <- readRDS(cache_file)
    B      <- basis_obj$B
    B_test <- basis_obj$B_test
  } else {
    n  <- nrow(Z_sc)
    Qc <- diag(1, n) - matrix(1/n, n, n)

    K_mat <- cpp_K(gg[, 1], gg[, 2], n)
    eiK   <- eigs_sym_cpp(Qc %*% K_mat %*% Qc, k0)

    V_mat <- sweep(eiK$vectors, 2, eiK$values, "/")
    K_one <- K_mat %*% rep(1/n, n)

    B      <- cpp_Kmatrix(k0, gg, gg,    K_one, V_mat, n, n)
    B_test <- cpp_Kmatrix(k0, gg, grids, K_one, V_mat, n, nrow(grids))

    saveRDS(list(B = B, B_test = B_test), cache_file)
  }

  B_list <- lapply(obs_idx, function(idxt) B[idxt, , drop = FALSE])

  r      <- ncol(B)
  eta_00 <- rep(0, r)

  var_Z <- var(as.vector(Z_sc), na.rm = TRUE)
  P_00  <- diag(c(var_Z, rep(var_Z/10, r - 1)))

  H  <- diag(1, r)
  U  <- diag(c(1, rep(0.1, r - 1)))
  K0 <- P_00

  sigma2_eps <- var_Z
  sig2_delta <- rep(var_Z / T, T)

  beta_t     <- rep(0, T)
  z_adj_list <- make_z_adj_list(z_list, beta_t)

  ll_old <- -Inf
  n_iter <- 0

  for (iter in 1:max_iter) {

    d_list <- make_d_list(sig2_delta, sigma2_eps, obs_idx)

    res_f <- forward_filter_cpp(H, U, B_list, d_list, z_adj_list, eta_00, K0)
    res_s <- backward_smooth_cross_cpp(H, res_f$eta_f, res_f$P_f, U)

    ll_new <- -0.5 * sum(res_f$log_parts + res_f$quad_parts)

    rd_res     <- mstep_rd_cpp(res_s$eta_s, res_s$P_s, B_list, z_adj_list, obs_idx, time_varying = TRUE)
    sig2_delta <- pmax(rd_res$rd2, 1e-12)

    beta_t <- sapply(seq_len(T), function(tt) {
      gt <- res_s$eta_s[[tt]]
      Bt <- B_list[[tt]]
      zt <- z_list[[tt]]
      mean(zt - as.vector(Bt %*% gt))
    })

    num <- 0
    den <- 0
    for (tt in 1:T) {
      gt <- res_s$eta_s[[tt]]
      Pt <- res_s$P_s[[tt]]
      Bt <- B_list[[tt]]
      zt <- z_list[[tt]]

      resid <- zt - beta_t[tt] - as.vector(Bt %*% gt)
      num <- num + sum(resid^2) + sum(diag(Bt %*% Pt %*% t(Bt)))
      den <- den + length(resid)
    }
    sigma2_eps <- max(num / den, 1e-12)

    z_adj_list <- make_z_adj_list(z_list, beta_t)

    HU_res <- mstep_HU_cpp(res_s$eta_s, res_s$P_s, res_s$P_cs)
    H_tmp  <- (HU_res$H + t(HU_res$H)) / 2
    U_tmp  <- (HU_res$U + t(HU_res$U)) / 2

    eU <- eigen(U_tmp, symmetric = TRUE)$values
    if (min(eU) < 1e-8) U_tmp <- U_tmp + (abs(min(eU)) + 1e-8) * diag(r)
    H <- H_tmp
    U <- U_tmp

    K0 <- update_K0_cpp(res_s$eta_s[[1]], res_s$P_s[[1]])
    K0 <- (K0 + t(K0)) / 2
    eK <- eigen(K0, symmetric = TRUE)$values
    if (min(eK) < 1e-8) K0 <- K0 + (abs(min(eK)) + 1e-8) * diag(r)

    diff_rel <- if (is.finite(ll_old) && ll_old != 0) abs((ll_new - ll_old) / ll_old) else Inf
    ll_old <- ll_new
    n_iter <- iter

    if (diff_rel < error) break
  }

  z_adj_list <- make_z_adj_list(z_list, beta_t)
  d_list     <- make_d_list(sig2_delta, sigma2_eps, obs_idx)

  res_f <- forward_filter_cpp(H, U, B_list, d_list, z_adj_list, eta_00, K0)
  res_s <- backward_smooth_cross_cpp(H, res_f$eta_f, res_f$P_f, U)

  pred_train_sc <- sapply(seq_len(T), function(tt) beta_t[tt] + as.vector(B %*% res_s$eta_s[[tt]]))
  pred_test_sc  <- sapply(seq_len(T), function(tt) beta_t[tt] + as.vector(B_test %*% res_s$eta_s[[tt]]))

  pred_train_in <- mu0 + sd0 * pred_train_sc
  pred_test_in  <- mu0 + sd0 * pred_test_sc

  mask_tr <- !is.na(Z)
  rmse_tr <- sqrt(mean((pred_train_in[mask_tr] - Z[mask_tr])^2))

  mask_te  <- !is.na(ytrue)
  rmspe_te <- sqrt(mean((pred_test_in[mask_te] - ytrue[mask_te])^2))

  pred_next_tr_sc <- matrix(NA_real_, nrow = m, ncol = T)
  for (tt in 1:(T - 1)) {
    eta_tp1_t  <- H %*% res_f$eta_f[[tt]]
    beta_tp1_t <- beta_t[tt]
    pred_next_tr_sc[, tt + 1] <- beta_tp1_t + as.vector(B %*% eta_tp1_t)
  }
  pred_next_tr_in <- mu0 + sd0 * pred_next_tr_sc
  mask_next_tr    <- !is.na(Z[, 2:T, drop = FALSE])
  rmspe_next_tr   <- sqrt(mean((pred_next_tr_in[, 2:T][mask_next_tr] - Z[, 2:T][mask_next_tr])^2))

  pred_next_te_sc <- matrix(NA_real_, nrow = m_t, ncol = T)
  for (tt in 1:(T - 1)) {
    eta_tp1_t  <- H %*% res_f$eta_f[[tt]]
    beta_tp1_t <- beta_t[tt]
    pred_next_te_sc[, tt + 1] <- beta_tp1_t + as.vector(B_test %*% eta_tp1_t)
  }
  pred_next_te_in <- mu0 + sd0 * pred_next_te_sc
  mask_next_te    <- !is.na(ytrue[, 2:T, drop = FALSE])
  rmspe_next_te   <- sqrt(mean((pred_next_te_in[, 2:T][mask_next_te] - ytrue[, 2:T][mask_next_te])^2))

  pred_next_tr_s_sc <- matrix(NA_real_, nrow = m, ncol = T)
  for (tt in 1:(T - 1)) {
    eta_tp1_t  <- H %*% res_s$eta_s[[tt]]
    beta_tp1_t <- beta_t[tt]
    pred_next_tr_s_sc[, tt + 1] <- beta_tp1_t + as.vector(B %*% eta_tp1_t)
  }
  pred_next_tr_s_in <- mu0 + sd0 * pred_next_tr_s_sc
  rmspe_next_tr_s   <- sqrt(mean((pred_next_tr_s_in[, 2:T][mask_next_tr] - Z[, 2:T][mask_next_tr])^2))

  data.frame(
    mu0 = mu0, sd0 = sd0,
    Train_RMSE = rmse_tr,
    Test_RMSPE_same_time = rmspe_te,
    Next_RMSPE_train_space = rmspe_next_tr,
    Next_RMSPE_test_space = rmspe_next_te,
    Sanity_Smoother_Next_train = rmspe_next_tr_s
  )
}

w_t <- readRDS("C:/Users/user/Desktop/ssm/data/Weather2K/weather2k_t.rds")
Z_all  <- w_t$Z
gg_all <- w_t$coords
rm(w_t); gc()

error    <- 1e-6
k0       <- 10
max_iter <- 200

# ----------------------------
# Exp1: sweep T (m fixed)
# ----------------------------
m_exp1      <- 1000
m_t_exp1    <- 300
T_set_exp1  <- c(200, 1000, 2920)
t0_set_exp1 <- c(1, 4000, 9000)
seed_set    <- c(1, 2, 3)

# ----------------------------
# Exp2: sweep m (T fixed)
# ----------------------------
T_exp2      <- 1000
m_set_exp2  <- c(300, 500, 1000, 2000)
t0_set_exp2 <- c(1, 9000)

# ----------------------------
# Build job grid
# ----------------------------
jobs <- list()

for (T in T_set_exp1) {
  for (t0 in t0_set_exp1) {
    for (seed in seed_set) {
      jobs[[length(jobs) + 1]] <- data.frame(
        exp_id = "exp1_T",
        T = T, t0 = t0, seed = seed,
        m = m_exp1, m_t = m_t_exp1, k0 = k0,
        stringsAsFactors = FALSE
      )
    }
  }
}

for (m in m_set_exp2) {
  m_t <- as.integer(round(0.3 * m))
  for (t0 in t0_set_exp2) {
    for (seed in seed_set) {
      jobs[[length(jobs) + 1]] <- data.frame(
        exp_id = "exp2_m",
        T = T_exp2, t0 = t0, seed = seed,
        m = m, m_t = m_t, k0 = k0,
        stringsAsFactors = FALSE
      )
    }
  }
}

jobs <- do.call(rbind, jobs)

# ----------------------------
# Output folders (NEW)
# ----------------------------
root_dir  <- "C:/Users/user/Desktop/ssm"
cache_dir <- file.path(root_dir, "cache")
out_dir   <- file.path(root_dir, "exp_runs")

dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir,   showWarnings = FALSE, recursive = TRUE)

out_rds <- file.path(out_dir, "results_exp1_exp2_raw.rds")
out_csv <- file.path(out_dir, "results_exp1_exp2_raw.csv")
log_file <- file.path(out_dir, "run_log.txt")

log_msg <- function(msg) {
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

# ----------------------------
# Resume if exists
# ----------------------------
results <- if (file.exists(out_rds)) readRDS(out_rds) else data.frame()

key_exists <- function(df, exp_id, T, t0, seed, m, m_t, k0) {
  if (nrow(df) == 0) return(FALSE)
  any(df$exp_id == exp_id &
        df$T == T & df$t0 == t0 & df$seed == seed &
        df$m == m & df$m_t == m_t & df$k0 == k0,
      na.rm = TRUE)
}

# ----------------------------
# Run all jobs (skip done, continue on error, checkpoint)
# ----------------------------
for (i in seq_len(nrow(jobs))) {

  job <- jobs[i, ]

  if (key_exists(results, job$exp_id, job$T, job$t0, job$seed, job$m, job$m_t, job$k0)) {
    log_msg(sprintf("Skip (done): %s | T=%d | t0=%d | seed=%d | m=%d | m_t=%d | k0=%d",
                    job$exp_id, job$T, job$t0, job$seed, job$m, job$m_t, job$k0))
    next
  }

  log_msg(sprintf("Running: %s | T=%d | t0=%d | seed=%d | m=%d | m_t=%d | k0=%d",
                  job$exp_id, job$T, job$t0, job$seed, job$m, job$m_t, job$k0))

  res <- tryCatch(
    run_one(
      Z_all = Z_all, gg_all = gg_all,
      T = job$T, t0 = job$t0, m = job$m, m_t = job$m_t, k0 = job$k0, seed = job$seed,
      error = error, max_iter = max_iter,
      cache_dir = cache_dir
    ),
    error = function(e) {
      data.frame(
        mu0 = NA_real_, sd0 = NA_real_,
        Train_RMSE = NA_real_,
        Test_RMSPE_same_time = NA_real_,
        Next_RMSPE_train_space = NA_real_,
        Next_RMSPE_test_space = NA_real_,
        Sanity_Smoother_Next_train = NA_real_,
        err_msg = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )

  if (!("err_msg" %in% names(res))) res$err_msg <- NA_character_

  res$exp_id <- job$exp_id
  res$T      <- job$T
  res$t0     <- job$t0
  res$seed   <- job$seed
  res$m      <- job$m
  res$m_t    <- job$m_t
  res$k0     <- job$k0

  results <- rbind(results, res)

  saveRDS(results, out_rds)
  write.csv(results, out_csv, row.names = FALSE)

  gc()
}

print(results)
