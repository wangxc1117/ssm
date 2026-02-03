library(Rcpp)
Rcpp::sourceCpp("em_use.cpp")
library(terra)
library(stringr)
library(Matrix)
source("function_fast.R")
library(astsa)


# read data
nc_dir = "C:/Users/user/Desktop/MERRA_2_temp/nc4/"
start_date = as.Date("2024-01-01")
end_date = as.Date("2024-01-14")

nc_files = list.files(nc_dir, pattern="\\.nc4$", full.names=TRUE)
fn_dates = as.Date(str_extract(basename(nc_files), "\\d{8}"), format="%Y%m%d")
files_sel = nc_files[fn_dates >= start_date & fn_dates <= end_date]
if (length(files_sel) == 0) stop("找不到符合日期的 .nc4 檔案")

r_stack = rast(files_sel)             
y1 = values(r_stack)                     # cell × time
gg = as.matrix(crds(r_stack, df=TRUE))   # cell × (lon, lat)


# parameter
error = 1e-06
r = bigK = 10
T = 100
m = 1000
m_t = 50
max_iter = 200
k0 = bigK

# train/test
set.seed(123)
pickm = sample(seq_len(nrow(y1)), m)
pickn2 = sample(seq_len(nrow(y1)), m_t)

y_train = y1[pickm, ]
y_test = y1[pickn2, ]
ytrue = y_test[, 1:T]
grids = gg[pickn2, ]

Z_raw = y_train[, 1:T]
Z = Z_raw

obs_idx = lapply(seq_len(T), function(t) which(!is.na(Z[, t])))
z_list = lapply(seq_len(T), function(t) {
  idxt = obs_idx[[t]]
  Z[idxt, t]
})

# fixed effect: time-varying intercept beta_t
z_var = var(as.vector(Z), na.rm = TRUE)
sigma2_eps = z_var
sig2_delta = rep(z_var / T, T)

beta_t = rep(0, T)

make_z_adj_list = function(z_list, beta_t) {
  lapply(seq_along(z_list), function(t) z_list[[t]] - beta_t[t])
}

z_adj_list = make_z_adj_list(z_list, beta_t)


# 有cache可以不用再算
# basis functions
tm_basis = system.time({
  n = nrow(Z)
  Qc = diag(1, n) - matrix(1/n, n, n)
  K_mat = cpp_K(gg[pickm, 1], gg[pickm, 2], n)
  eiK = eigs_sym_cpp(Qc %*% K_mat %*% Qc, k0)
  V_mat = sweep(eiK$vectors, 2, eiK$values, "/")
  B_full = cpp_Kmatrix(k0, gg[pickm, ], gg[pickm, ], K_mat %*% rep(1/n, n), V_mat, n, n)
  B_test_full = cpp_Kmatrix(k0, gg[pickm, ], grids, K_mat %*% rep(1/n, n), V_mat, n, nrow(grids))
  B      = B_full
  B_test = B_test_full
  r      = ncol(B)
  B_list = lapply(obs_idx, function(idxt) B[idxt, , drop=FALSE])
})
cat(sprintf("Basis 計算耗時：%.2f 秒\n", tm_basis["elapsed"]))


tm_basis = system.time({
  n  = nrow(Z)
  Qc = diag(1, n) - matrix(1/n, n, n)

  K_mat = cpp_K(gg[pickm, 1], gg[pickm, 2], n)
  eiK   = eigs_sym_cpp(Qc %*% K_mat %*% Qc, k0)
  V_mat = sweep(eiK$vectors, 2, eiK$values, "/")
  K_one = K_mat %*% rep(1/n, n)

  B_full = cpp_Kmatrix(k0, gg[pickm, ], gg[pickm, ], K_one, V_mat, n, n)
  B_test_full = cpp_Kmatrix(k0, gg[pickm, ], grids, K_one, V_mat, n, nrow(grids))

  B      = B_full
  B_test = B_test_full
  r      = ncol(B)

  B_list = lapply(obs_idx, function(idxt) B[idxt, , drop = FALSE])
})
cat(sprintf("Basis 計算耗時：%.2f 秒\n", tm_basis["elapsed"]))

# ---------- 打包 + 存檔 ----------
cache_dir  = "C:/Users/user/Desktop/ssm/cache"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

cache_file = file.path(
  cache_dir,
  sprintf("basis_m%d_k%d_T%d_seed123_%s_%s.rds",
          m, k0, T, format(start_date, "%Y%m%d"), format(end_date, "%Y%m%d"))
)

vars = c(
  "B","B_test","r","k0","V_mat","K_one",
  "pickm","pickn2","Z","obs_idx","z_list"
)
if (exists("B_list")) vars = c(vars, "B_list")

basis_obj = mget(vars, inherits = TRUE)
saveRDS(basis_obj, cache_file)
cat("已存到：", normalizePath(cache_file), "\n")


cache_dir  = "C:/Users/user/Desktop/ssm/cache"
cache_file = file.path(cache_dir,
  "basis_m1000_k10_T100_seed123_20240101_20240114.rds"  # 改成你的實際檔名
)

stopifnot(file.exists(cache_file))

basis_obj = readRDS(cache_file)
list2env(basis_obj, envir = .GlobalEnv)
cat("B first col all ones? ", isTRUE(all(B[,1] == 1)), "\n")
stopifnot(!isTRUE(all(B[,1] == 1)))

B_list = lapply(obs_idx, function(idxt) B[idxt, , drop = FALSE])


# 參數初始化
r        = ncol(B)
eta_00   = rep(0, r)
var_Z    = var(as.vector(Z), na.rm = TRUE)
P_00     = diag(c(var_Z, rep(var_Z/10, r-1)))

H        = diag(1, r)
U        = diag(c(1, rep(0.1, r-1)))
K0       = P_00

z_var      = var_Z
sigma2_eps = z_var
sig2_delta = rep(z_var / T, T)

beta_t = rep(0, T)

make_z_adj_list = function(z_list, beta_t) {
  lapply(seq_along(z_list), function(t) z_list[[t]] - beta_t[t])
}

z_adj_list = make_z_adj_list(z_list, beta_t)

make_d_list = function(sig2_delta, sigma2_eps, obs_idx) {
  lapply(seq_along(obs_idx), function(t)
    rep(sigma2_eps + sig2_delta[t],
        length(obs_idx[[t]]))
  )
}

# --- EM 迴圈（fixed effect: time-varying intercept；不鎖 state intercept） ---
tol_rel = error
ll_old  = -Inf

for (iter in 1:max_iter) {

  ## — E-step —
  d_list = make_d_list(sig2_delta, sigma2_eps, obs_idx)
  res_f  = forward_filter_cpp(H, U, B_list, d_list, z_adj_list, eta_00, K0)
  res_s  = backward_smooth_cross_cpp(H, res_f$eta_f, res_f$P_f, U)

  ll_new = -0.5 * sum(res_f$log_parts + res_f$quad_parts)

  ## — M-step —

  # 1) 更新 σ²_{δ,t}（用扣掉固定效應後的 z_adj_list）
  rd_res     = mstep_rd_cpp(res_s$eta_s, res_s$P_s,
                            B_list, z_adj_list, obs_idx,
                            time_varying = TRUE)
  sig2_delta = pmax(rd_res$rd2, 1e-12)

  # 2) 更新 beta_t（time-varying intercept）
  beta_t = sapply(seq_len(T), function(t) {
    gt = res_s$eta_s[[t]]
    Bt = B_list[[t]]
    zt = z_list[[t]]
    mean(zt - as.vector(Bt %*% gt))
  })

  # 3) 更新 measurement error variance sigma2_eps
  num = 0
  den = 0
  for (t in 1:T) {
    gt = res_s$eta_s[[t]]
    Pt = res_s$P_s[[t]]
    Bt = B_list[[t]]
    zt = z_list[[t]]

    resid = zt - beta_t[t] - as.vector(Bt %*% gt)
    num = num + sum(resid^2) + sum(diag(Bt %*% Pt %*% t(Bt)))
    den = den + length(resid)
  }
  sigma2_eps = max(num / den, 1e-12)

  # 4) 更新 z_adj_list（供下一輪 E-step）
  z_adj_list = make_z_adj_list(z_list, beta_t)

  # 5) 更新 H, U 並對稱化（不再鎖 intercept 維度）
  HU_res = mstep_HU_cpp(res_s$eta_s, res_s$P_s, res_s$P_cs)
  H_tmp  = (HU_res$H + t(HU_res$H)) / 2
  U_tmp  = (HU_res$U + t(HU_res$U)) / 2

  # 確保 U_tmp 正定
  eU = eigen(U_tmp, symmetric = TRUE)$values
  if (min(eU) < 1e-8) {
    U_tmp = U_tmp + (abs(min(eU)) + 1e-8) * diag(r)
  }

  H = H_tmp
  U = U_tmp

  # 6) 更新 K0
  K0 = update_K0_cpp(res_s$eta_s[[1]], res_s$P_s[[1]])
  K0 = (K0 + t(K0)) / 2
  eK = eigen(K0, symmetric = TRUE)$values
  if (min(eK) < 1e-8) {
    K0 = K0 + (abs(min(eK)) + 1e-8) * diag(r)
  }

  ## — 收斂判定 —
  diff_rel = if (is.finite(ll_old) && ll_old != 0) abs((ll_new - ll_old) / ll_old) else Inf
  cat(sprintf("Iter %03d | logLik = %.6f | Δrel = %.3e | sigma2_eps = %.6e | beta_t range = [%.4f, %.4f]\n",
              iter, ll_new, diff_rel, sigma2_eps, min(beta_t), max(beta_t)))
  if (diff_rel < tol_rel) {
    cat("EM 收斂於第", iter, "次迭代\n")
    break
  }
  ll_old = ll_new
}

# 最終平滑（用扣掉固定效應的 z_adj_list）
z_adj_list = make_z_adj_list(z_list, beta_t)
d_list = make_d_list(sig2_delta, sigma2_eps, obs_idx)
res_f  = forward_filter_cpp(H, U, B_list, d_list, z_adj_list, eta_00, K0)
res_s  = backward_smooth_cross_cpp(H, res_f$eta_f, res_f$P_f, U)


# Training
y_train_pred = sapply(seq_len(T), function(t) beta_t[t] + as.vector(B %*% res_s$eta_s[[t]]))
mask_tr      = !is.na(Z)
mse_tr       = mean((y_train_pred[mask_tr] - Z[mask_tr])^2)
cat(sprintf("Training MSE: %.6f | RMSE: %.6f\n", mse_tr, sqrt(mse_tr)))

# Test
y_test_pred  = sapply(seq_len(T), function(t) beta_t[t] + as.vector(B_test %*% res_s$eta_s[[t]]))
mask_te      = !is.na(ytrue)
mspe         = mean((y_test_pred[mask_te] - ytrue[mask_te])^2)
cat(sprintf("Test MSPE: %.6f | RMSPE: %.6f\n", mspe, sqrt(mspe)))


eta_next   = H %*% res_s$eta_s[[T]]        # η_{T+1|T}
beta_next  = beta_t[T]                     # 最小改動：β_{T+1|T} := β_T（random-walk / persistence）
y_em_pred  = beta_next + as.vector(B %*% eta_next)

y_em_true  = y1[pickm, T + 1]

mask       = !is.na(y_em_true)
mspe_em    = mean((y_em_pred[mask] - y_em_true[mask])^2)
rmspe_em   = sqrt(mspe_em)
cat(sprintf("EM 全點 Next-step MSPE = %.6f | RMSPE = %.6f\n", mspe_em, rmspe_em))



