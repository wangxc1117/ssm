// em_use.cpp
// RcppArmadillo：SMW、截斷特徵分解、前向濾波、後向平滑(含交叉共變)、M-step(估 H/U/K0/variance)
// 加入 safe_chol + pinv 後援，避免 chol fail

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;

// ---------- Utilities ----------

// 安全 chol：對稱化 + nugget + retry
static inline void safe_chol(const arma::mat& A_in,
                             arma::mat& C_out,
                             double base_jitter = 1e-8,
                             int max_try = 6) {
  arma::mat A = 0.5 * (A_in + A_in.t());
  double jitter = base_jitter;
  bool ok = false;
  for (int i = 0; i < max_try; ++i) {
    try {
      C_out = arma::chol(A);
      ok = true;
      break;
    } catch (...) {
      A.diag() += jitter;
      jitter *= 10.0;
    }
  }
  if (!ok) Rcpp::stop("safe_chol: decomposition failed.");
}

// pseudo-inverse & pseudo-logdet for symmetric matrix
static inline void pinv_logdet_sym(const arma::mat& A_in,
                                   arma::mat& A_inv,
                                   double& logdet,
                                   double tol = 1e-8) {
  arma::mat A = 0.5 * (A_in + A_in.t());
  arma::vec eigval; arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, A);
  arma::uvec keep = arma::find(eigval > tol);
  if (keep.n_elem == 0) {
    // last resort
    arma::mat A2 = A + tol * arma::eye(A.n_rows, A.n_cols);
    arma::eig_sym(eigval, eigvec, A2);
    keep = arma::find(eigval > tol);
    if (keep.n_elem == 0) Rcpp::stop("pinv_logdet_sym: singular matrix.");
  }
  arma::vec lam = eigval(keep);
  arma::mat V   = eigvec.cols(keep);
  A_inv  = V * arma::diagmat(1.0 / lam) * V.t();
  logdet = arma::sum(arma::log(lam));
}

// 用 safe_chol 的逆
static inline arma::mat safe_chol_inv(const arma::mat& A_in) {
  arma::mat C;
  safe_chol(A_in, C);
  return arma::solve(arma::trimatu(C),
                     arma::solve(arma::trimatl(C.t()),
                                 arma::eye<arma::mat>(C.n_rows, C.n_cols)));
}

// ===== 0. SMW 函式（保留） =====
static arma::vec smw_vec_cpp(const arma::vec& d,
                             const arma::mat& B,
                             const arma::mat& P,
                             const arma::vec& v) {
  arma::vec d_inv = 1.0 / d;
  arma::mat BtDinv = B.t();
  BtDinv.each_row() %= d_inv.t();
  arma::mat M = arma::inv(P) + BtDinv * B;
  arma::mat M_inv = arma::inv(M);
  arma::vec v1 = d_inv % v;
  arma::vec v2 = BtDinv * v;
  return v1 - d_inv % (B * (M_inv * v2));
}

// ===== 1. 截斷特徵分解 =====
// [[Rcpp::export]]
List eigs_sym_cpp(const arma::mat& M, arma::uword k) {
  arma::vec all_vals; arma::mat all_vecs;
  arma::eig_sym(all_vals, all_vecs, M);
  arma::uword p = all_vals.n_elem;
  return List::create(
    Named("values")  = all_vals.subvec(p - k, p - 1),
    Named("vectors") = all_vecs.cols(p - k, p - 1)
  );
}

// ===== 2. 前向濾波 =====
// [[Rcpp::export]]
List forward_filter_cpp(
    const arma::mat& H,
    const arma::mat& U,
    const List& B_list,
    const List& d_list,
    const List& z_list,
    const arma::vec& eta0,
    const arma::mat& P0) {
  
  int T = B_list.size();
  List eta_f(T), P_f(T);
  arma::vec log_parts(T), quad_parts(T);
  
  arma::vec eta_prev = eta0;
  arma::mat P_prev   = P0;
  
  for (int t = 0; t < T; ++t) {
    arma::mat B_t = as<arma::mat>(B_list[t]);
    arma::vec d_t = as<arma::vec>(d_list[t]);
    arma::vec z_t = as<arma::vec>(z_list[t]);
    
    // 1) Prediction
    arma::vec eta_pr = H * eta_prev;
    arma::mat P_pr   = H * P_prev * H.t() + U;
    
    // 2) P_pr_inv & logdetP  (safe chol -> pinv fallback)
    arma::mat P_pr_inv; double logdetP = 0.0;
    try {
      arma::mat chol_P;
      safe_chol(P_pr, chol_P);
      logdetP   = 2.0 * arma::sum(arma::log(chol_P.diag()));
      P_pr_inv  = arma::solve(arma::trimatu(chol_P),
                              arma::solve(arma::trimatl(chol_P.t()),
                                          arma::eye<arma::mat>(chol_P.n_rows, chol_P.n_cols)));
    } catch (...) {
      pinv_logdet_sym(P_pr, P_pr_inv, logdetP, 1e-8);
    }
    
    // 3) D^{-1}B 與 G
    arma::vec d_inv = 1.0 / d_t;
    arma::mat DinvB = B_t.each_col() % d_inv;
    arma::mat G     = B_t.t() * DinvB;
    
    // 4) M = P_pr_inv + G
    arma::mat Mmat   = P_pr_inv + G;
    arma::mat Mmat_inv; double logdetM = 0.0;
    try {
      arma::mat chol_M;
      safe_chol(Mmat, chol_M);
      logdetM   = 2.0 * arma::sum(arma::log(chol_M.diag()));
      Mmat_inv  = arma::solve(arma::trimatu(chol_M),
                              arma::solve(arma::trimatl(chol_M.t()),
                                          arma::eye<arma::mat>(chol_M.n_rows, chol_M.n_cols)));
    } catch (...) {
      pinv_logdet_sym(Mmat, Mmat_inv, logdetM, 1e-8);
    }
    
    // 5) log|R_t|
    log_parts[t] = arma::sum(arma::log(d_t)) + logdetP + logdetM;
    
    // 6) quadratic term
    arma::vec resid = z_t - B_t * eta_pr;
    arma::vec v1    = d_inv % resid;
    arma::vec tmp   = B_t.t() * v1;
    arma::vec x     = Mmat_inv * tmp;
    arma::vec v2    = d_inv % (B_t * x);
    arma::vec Sinv_r= v1 - v2;
    quad_parts[t]   = arma::dot(resid, Sinv_r);
    
    // 7) Posterior mean/cov
    arma::vec eta_up = eta_pr + P_pr * B_t.t() * Sinv_r;
    
    arma::mat E      = Mmat_inv * G;
    arma::mat BE     = B_t * E;
    arma::mat SinvB  = DinvB - (BE.each_col() % d_inv);
    arma::mat P_up   = P_pr - P_pr * (B_t.t() * SinvB) * P_pr;
    
    eta_f[t] = eta_up;
    P_f[t]   = P_up;
    eta_prev = eta_up;
    P_prev   = P_up;
  }
  
  return List::create(
    Named("eta_f")      = eta_f,
    Named("P_f")        = P_f,
    Named("log_parts")  = log_parts,
    Named("quad_parts") = quad_parts
  );
}

// ===== 3. 後向平滑 (含 cross-cov) =====
// [[Rcpp::export]]
List backward_smooth_cross_cpp(
    const arma::mat& H,
    const List& eta_f,
    const List& P_f,
    const arma::mat& U) {
  
  int T = eta_f.size();
  List eta_s(T), P_s(T), P_cs(T);
  std::vector<arma::mat> J_list(T);
  
  arma::mat P0tmp = as<arma::mat>(P_f[0]);
  arma::uword r = P0tmp.n_rows, c = P0tmp.n_cols;
  
  eta_s[T-1] = eta_f[T-1];
  P_s[T-1]   = P_f[T-1];
  P_cs[T-1]  = arma::mat(r, c, arma::fill::zeros);
  
  for (int t = T-2; t >= 0; --t) {
    arma::vec eta_ft = as<arma::vec>(eta_f[t]);
    arma::mat P_ft   = as<arma::mat>(P_f[t]);
    
    arma::mat P_pr   = H * P_ft * H.t() + U;
    arma::mat P_pr_inv;
    try {
      P_pr_inv = safe_chol_inv(P_pr);
    } catch (...) {
      double dummy;
      pinv_logdet_sym(P_pr, P_pr_inv, dummy, 1e-8);
    }
    arma::mat J = P_ft * H.t() * P_pr_inv;
    J_list[t]   = J;
    
    arma::vec eta_sp = as<arma::vec>(eta_s[t+1]);
    arma::mat P_sp   = as<arma::mat>(P_s[t+1]);
    
    arma::vec eta_st = eta_ft + J * (eta_sp - H * eta_ft);
    arma::mat P_st   = P_ft + J * (P_sp - P_pr) * J.t();
    
    eta_s[t] = eta_st;
    P_s[t]   = P_st;
  }
  
  // cross-cov
  P_cs[0] = arma::mat(r, c, arma::fill::zeros);
  for (int t = 1; t < T; ++t) {
    arma::mat J_tm1 = J_list[t-1];
    arma::mat P_t   = as<arma::mat>(P_s[t]);
    P_cs[t]         = J_tm1 * P_t;
  }
  
  return List::create(
    Named("eta_s") = eta_s,
    Named("P_s")   = P_s,
    Named("P_cs")  = P_cs
  );
}

// ===== 4a. M-step：fine-scale / measurement 變異 =====
// [[Rcpp::export]]
List mstep_variance_cpp(
    const List& eta_s,
    const List& P_s,
    const List& B_list,
    const List& z_list) {
  
  int T = eta_s.size();
  std::vector<arma::vec> eta_vec(T);
  std::vector<arma::mat> P_mat(T);
  std::vector<arma::mat> B_mat(T);
  std::vector<arma::vec> z_vec(T);
  
  for (int t = 0; t < T; ++t) {
    eta_vec[t] = as<arma::vec>(eta_s[t]);
    P_mat[t]   = as<arma::mat>(P_s[t]);
    B_mat[t]   = as<arma::mat>(B_list[t]);
    z_vec[t]   = as<arma::vec>(z_list[t]);
  }
  
  arma::vec sig2_delta(T);
  double sum_eps = 0.0;
  double total_n = 0.0;
  
#pragma omp parallel for reduction(+:sum_eps,total_n)
  for (int t = 0; t < T; ++t) {
    const arma::vec& eta_st = eta_vec[t];
    const arma::mat& P_st   = P_mat[t];
    arma::mat Mmat = P_st + eta_st * eta_st.t();
    sig2_delta[t]  = arma::mean(arma::diagvec(Mmat));
    
    const arma::mat& B_t = B_mat[t];
    const arma::vec& z_t = z_vec[t];
    arma::vec resid      = z_t - B_t * eta_st;
    double eps_part      = arma::dot(resid, resid) + arma::trace(B_t * P_st * B_t.t());
    int    n_part        = resid.n_elem;
    
    sum_eps  += eps_part;
    total_n  += n_part;
  }
  
  double sigma2_eps = sum_eps / total_n;
  
  return List::create(
    Named("sigma2_eps") = sigma2_eps,
    Named("sig2_delta") = sig2_delta
  );
}

// ===== 4b. M-step：估 H、U =====
// [[Rcpp::export]]
List mstep_HU_cpp(
    const List& eta_s,
    const List& P_s,
    const List& P_cs) {
  
  int T = eta_s.size();
  arma::uword r2 = as<arma::vec>(eta_s[0]).n_elem;
  
  arma::mat SumK = arma::zeros(r2, r2);  // t=0..T-2
  arma::mat SumL = arma::zeros(r2, r2);  // t=1..T-1
  
  for (int t = 0; t < T; ++t) {
    arma::vec gt = as<arma::vec>(eta_s[t]);
    arma::mat Pt = as<arma::mat>(P_s[t]);
    arma::mat Kt = Pt + gt * gt.t();
    if (t < T-1) SumK += Kt;
    
    if (t > 0) {
      arma::vec gt1 = as<arma::vec>(eta_s[t-1]);
      arma::mat Pcs = as<arma::mat>(P_cs[t]);
      arma::mat Lt  = Pcs + gt * gt1.t();
      SumL += Lt;
    }
  }
  
  // 正則化防奇異
  arma::mat H_new = SumL * arma::inv(SumK + 1e-10 * arma::eye(r2, r2));
  
  arma::mat U_new = arma::zeros(r2, r2);
  for (int t = 1; t < T; ++t) {
    arma::vec gt = as<arma::vec>(eta_s[t]);
    arma::mat Pt = as<arma::mat>(P_s[t]);
    arma::mat Kt = Pt + gt * gt.t();
    
    arma::vec gt1 = as<arma::vec>(eta_s[t-1]);
    arma::mat Pcs = as<arma::mat>(P_cs[t]);
    arma::mat Lt  = Pcs + gt * gt1.t();
    
    U_new += (Kt - H_new * Lt.t());
  }
  U_new /= static_cast<double>(T - 1);
  
  return List::create(
    Named("H") = H_new,
    Named("U") = U_new
  );
}

// ===== 4c. (可選) M-step：估 K0 =====
// [[Rcpp::export]]
arma::mat update_K0_cpp(const arma::vec& g0_s, const arma::mat& P0_s) {
  return P0_s + g0_s * g0_s.t();
}

// ===== 4d. M-step：論文式 (19)/(21) 的 r_{d,t}^2 =====
// 假設 V_{d,t} = I_{n_t}；若要一般化，可再加一個 Vd_list 輸入
// [[Rcpp::export]]
List mstep_rd_cpp(const List& eta_s,
                  const List& P_s,
                  const List& B_list,
                  const List& z_list,
                  const List& obs_idx_list,
                  const bool time_varying = true) {
  int T = eta_s.size();
  arma::vec rd2(T, arma::fill::zeros);
  double num_sum = 0.0;
  double den_sum = 0.0;
  
  for (int t = 0; t < T; ++t) {
    arma::vec gt = as<arma::vec>(eta_s[t]);
    arma::mat Pt = as<arma::mat>(P_s[t]);
    arma::mat Bt = as<arma::mat>(B_list[t]);
    
    // 取出該期觀測的索引長度 (n_t)
    IntegerVector idv = obs_idx_list[t];
    int nt = idv.size();
    
    // d_t|T = z_t - B_t g_t|T   (假設 X b_t 無)
    arma::vec zt   = as<arma::vec>(z_list[t]);
    arma::vec dt_s = zt - Bt * gt;
    
    // R_t|T ≈ B_t P_t|T B_t'    (論文式中 R_t|T)
    arma::mat Rt_s = Bt * Pt * Bt.t();
    
    double trace_part = arma::trace(Rt_s + dt_s * dt_s.t());
    if (time_varying) {
      rd2[t] = trace_part / (double)nt;
    } else {
      num_sum += trace_part;
      den_sum += (double)nt;
    }
  }
  
  if (!time_varying) {
    double rd_const = num_sum / den_sum;
    rd2.fill(rd_const);
  }
  
  return List::create(
    Named("rd2") = rd2
  );
}
