library(Rcpp)
Rcpp::sourceCpp("em_use.cpp")
library(Matrix)
source("function_fast.R")

set.seed(123)

# -----------------------------
# Global settings
# -----------------------------
error    <- 1e-6
k0       <- 10
T        <- 200
m        <- 1000
m_t      <- 300
max_iter <- 200

rds_path <- "C:/Users/user/Desktop/ssm/data/Weather2K/weather2k_t.rds"

cache_dir <- "C:/Users/user/Desktop/ssm/cache"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

cache_file <- file.path(
  cache_dir,
  sprintf("basis_weather2k_m%d_k%d_T%d_seed123.rds", m, k0, T)
)

# -----------------------------
# Read Weather2K (light version)
# -----------------------------
w_t <- readRDS(rds_path)
Z_all  <- w_t$Z
gg_all <- w_t$coords
rm(w_t); gc()

stopifnot(nrow(Z_all) == nrow(gg_all))
stopifnot(ncol(Z_all) >= (T + 1))  # Next-step needs T+1

# -----------------------------
# Train / test split (space split)
# -----------------------------
pickm  <- sample(seq_len(nrow(Z_all)), m)
pool   <- setdiff(seq_len(nrow(Z_all)), pickm)
pickn2 <- sample(pool, m_t)

Z     <- Z_all[pickm,  1:T, drop = FALSE]
ytrue <- Z_all[pickn2, 1:T, drop = FALSE]

gg    <- gg_all[pickm,  , drop = FALSE]
grids <- gg_all[pickn2, , drop = FALSE]

obs_idx <- lapply(seq_len(T), function(t) which(!is.na(Z[, t])))
z_list  <- lapply(seq_len(T), function(t) {
  idxt <- obs_idx[[t]]
  Z[idxt, t]
})

make_z_adj_list <- function(z_list, beta_t) {
  lapply(seq_along(z_list), function(t) z_list[[t]] - beta_t[t])
}

make_d_list <- function(sig2_delta, sigma2_eps, obs_idx) {
  lapply(seq_along(obs_idx), function(t) {
    rep(sigma2_eps + sig2_delta[t], length(obs_idx[[t]]))
  })
}

# -----------------------------
# Basis (cacheable)
# -----------------------------
if (file.exists(cache_file)) {

  basis_obj <- readRDS(cache_file)
  stopifnot(all(c("B","B_test","k0","pickm","pickn2") %in% names(basis_obj)))
  B      <- basis_obj$B
  B_test <- basis_obj$B_test

} else {

  tm_basis <- system.time({
    n  <- nrow(Z)
    Qc <- diag(1, n) - matrix(1/n, n, n)

    # kernel eigenspace on training coords
    K_mat <- cpp_K(gg[, 1], gg[, 2], n)
    eiK   <- eigs_sym_cpp(Qc %*% K_mat %*% Qc, k0)

    V_mat <- sweep(eiK$vectors, 2, eiK$values, "/")
    K_one <- K_mat %*% rep(1/n, n)

    # Basis on train + on test space
    B      <- cpp_Kmatrix(k0, gg, gg,    K_one, V_mat, n, n)
    B_test <- cpp_Kmatrix(k0, gg, grids, K_one, V_mat, n, nrow(grids))
  })

  cat(sprintf("Basis 計算耗時：%.2f 秒\n", tm_basis["elapsed"]))

  saveRDS(
    list(
      B = B, B_test = B_test, k0 = k0,
      pickm = pickm, pickn2 = pickn2
    ),
    cache_file
  )
  cat("已存到：", normalizePath(cache_file), "\n")
}

stopifnot(!isTRUE(all(B[, 1] == 1)))

B_list <- lapply(obs_idx, function(idxt) B[idxt, , drop = FALSE])

# -----------------------------
# Init parameters
# -----------------------------
r      <- ncol(B)
eta_00 <- rep(0, r)

var_Z <- var(as.vector(Z), na.rm = TRUE)
P_00  <- diag(c(var_Z, rep(var_Z/10, r - 1)))

H  <- diag(1, r)
U  <- diag(c(1, rep(0.1, r - 1)))
K0 <- P_00

sigma2_eps <- var_Z
sig2_delta <- rep(var_Z / T, T)

beta_t     <- rep(0, T)
z_adj_list <- make_z_adj_list(z_list, beta_t)

# -----------------------------
# EM loop
# -----------------------------
tol_rel <- error
ll_old  <- -Inf

for (iter in 1:max_iter) {

  d_list <- make_d_list(sig2_delta, sigma2_eps, obs_idx)

  res_f <- forward_filter_cpp(H, U, B_list, d_list, z_adj_list, eta_00, K0)
  res_s <- backward_smooth_cross_cpp(H, res_f$eta_f, res_f$P_f, U)

  ll_new <- -0.5 * sum(res_f$log_parts + res_f$quad_parts)

  rd_res     <- mstep_rd_cpp(res_s$eta_s, res_s$P_s, B_list, z_adj_list, obs_idx, time_varying = TRUE)
  sig2_delta <- pmax(rd_res$rd2, 1e-12)

  beta_t <- sapply(seq_len(T), function(t) {
    gt <- res_s$eta_s[[t]]
    Bt <- B_list[[t]]
    zt <- z_list[[t]]
    mean(zt - as.vector(Bt %*% gt))
  })

  num <- 0
  den <- 0
  for (t in 1:T) {
    gt <- res_s$eta_s[[t]]
    Pt <- res_s$P_s[[t]]
    Bt <- B_list[[t]]
    zt <- z_list[[t]]

    resid <- zt - beta_t[t] - as.vector(Bt %*% gt)
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
  cat(sprintf(
    "Iter %03d | logLik = %.6f | Δrel = %.3e | sigma2_eps = %.6e | beta_t range = [%.4f, %.4f]\n",
    iter, ll_new, diff_rel, sigma2_eps, min(beta_t), max(beta_t)
  ))

  if (diff_rel < tol_rel) {
    cat("EM 收斂於第", iter, "次迭代\n")
    break
  }
  ll_old <- ll_new
}

# -----------------------------
# Final filter + smoother
# -----------------------------
z_adj_list <- make_z_adj_list(z_list, beta_t)
d_list     <- make_d_list(sig2_delta, sigma2_eps, obs_idx)

res_f <- forward_filter_cpp(H, U, B_list, d_list, z_adj_list, eta_00, K0)
res_s <- backward_smooth_cross_cpp(H, res_f$eta_f, res_f$P_f, U)

# -----------------------------
# Metrics (4 targets)
# -----------------------------

# 1) Train RMSE (in-sample, train space, use observed mask)
pred_train_in <- sapply(seq_len(T), function(t) beta_t[t] + as.vector(B %*% res_s$eta_s[[t]]))
mask_tr       <- !is.na(Z)
rmse_tr       <- sqrt(mean((pred_train_in[mask_tr] - Z[mask_tr])^2))
cat(sprintf("Train RMSE = %.6f\n", rmse_tr))

# 2) Test RMSPE (same time, different space)
pred_test_in <- sapply(seq_len(T), function(t) beta_t[t] + as.vector(B_test %*% res_s$eta_s[[t]]))
mask_te      <- !is.na(ytrue)
rmspe_te     <- sqrt(mean((pred_test_in[mask_te] - ytrue[mask_te])^2))
cat(sprintf("Test RMSPE (same time, test space) = %.6f\n", rmspe_te))

# 3) Next-step RMSPE (same space, next time) — one-step-ahead forecast from FILTER
pred_next_tr <- matrix(NA_real_, nrow = m, ncol = T)
for (t in 1:(T - 1)) {
  eta_tp1_t  <- H %*% res_f$eta_f[[t]]
  beta_tp1_t <- beta_t[t]  # persistence
  pred_next_tr[, t + 1] <- beta_tp1_t + as.vector(B %*% eta_tp1_t)
}
mask_next_tr <- !is.na(Z[, 2:T, drop = FALSE])
rmspe_next_tr <- sqrt(mean((pred_next_tr[, 2:T][mask_next_tr] - Z[, 2:T][mask_next_tr])^2))
cat(sprintf("Next-step RMSPE (train space) = %.6f\n", rmspe_next_tr))

# 4) Next-step + test space RMSPE — one-step-ahead forecast projected to test space
pred_next_te <- matrix(NA_real_, nrow = m_t, ncol = T)
for (t in 1:(T - 1)) {
  eta_tp1_t  <- H %*% res_f$eta_f[[t]]
  beta_tp1_t <- beta_t[t]
  pred_next_te[, t + 1] <- beta_tp1_t + as.vector(B_test %*% eta_tp1_t)
}
mask_next_te <- !is.na(ytrue[, 2:T, drop = FALSE])
rmspe_next_te <- sqrt(mean((pred_next_te[, 2:T][mask_next_te] - ytrue[, 2:T][mask_next_te])^2))
cat(sprintf("Next-step + test space RMSPE = %.6f\n", rmspe_next_te))
