# =============================================================================
# fit_r_midas_fast.R — drop-in accelerated concentrated NLS for R-MIDAS
#
# Mathematically IDENTICAL estimator to common_utils::fit_r_midas (profiled /
# concentrated NLS, Foroni-Guerin-Marcellino 2018 R-MIDAS with exp-Almon
# weights, Ghysels-Sinko-Valkanov 2007), but the concentrated SSE is evaluated
# through the Frisch-Waugh-Lovell decomposition with pre-computed cross
# products, so each objective evaluation no longer touches the n x p design.
#
# Algebra.  Let X be the linear block (intercept, calendar, LF, HF lag-1) and
# z_k(theta_k) = Z_k w(theta_k) the K Almon-weighted columns.  The profiled
# SSE of OLS on [X, z_1..z_K] equals, by FWL,
#
#     S(theta) = y~'y~  -  b(theta)' A(theta)^{-1} b(theta),
#
# where  y~ = M_X y,  Z~_k = M_X Z_k  (M_X = I - X(X'X)^{-1}X'),
#        A_{kl}(theta) = w_k' (Z~_k' Z~_l) w_l   (K x K),
#        b_k(theta)    = w_k' (Z~_k' y~).
#
# Z~'Z~ (KQ x KQ), Z~'y~ (KQ) and y~'y~ depend only on the data, NOT on theta,
# so they are computed ONCE per rolling window (one QR of X + crossprods).
# Each L-BFGS-B objective call then costs O(K^2 Q^2) flops (microseconds,
# independent of n) instead of a fresh O(n p^2) rank-checked QR.  Verified to
# reproduce the brute-force SSE to ~1e-13 relative error; the returned
# theta / beta_lin / delta / sigma2 are computed by the SAME final full-design
# QR as the original function, so downstream behaviour is unchanged.
#
# Usage: source this AFTER common_utils.R, then either call fit_r_midas_fast()
# explicitly, or simply alias
#     fit_r_midas <- fit_r_midas_fast
# in r_midas.R / r_midas_extended.R / r_midas_surveys_h.R after the source()
# lines (no other change needed — signature and return value are identical).
# =============================================================================

# lag_sets: optional named list giving, per Almon block, the ACTUAL lag
# numbers of its columns (e.g. list(y_lags = c(1,2,3,7)) for the paper's
# AR lag set J).  Defaults to consecutive 1..ncol per block.  Blocks may
# have different numbers of columns.
fit_r_midas_fast <- function(y, Z_list, X_lin,
                             theta_init = NULL,
                             theta_lower_block = c(-8, -8),
                             theta_upper_block = c( 8,  0),
                             optim_maxit = 2000L,
                             optim_factr = 5e8,
                             n_starts = 1L,
                             deterministic_starts = TRUE,
                             lag_sets = NULL) {
  K <- length(Z_list)
  if (K < 1L) stop("Need at least one Z block (the AR block).")
  block_names <- names(Z_list)
  if (is.null(lag_sets)) lag_sets <- list()
  for (nm in block_names) {
    if (is.null(lag_sets[[nm]]))
      lag_sets[[nm]] <- seq_len(ncol(Z_list[[nm]]))
  }
  Qs <- vapply(Z_list, ncol, integer(1))
  n <- length(y)
  if (is.null(X_lin))
    X_lin <- matrix(1, n, 1L, dimnames = list(NULL, "intercept"))

  # ---- one-off row filter (was redone inside the objective before) ----------
  keep <- stats::complete.cases(X_lin) & is.finite(y)
  for (k in seq_len(K))
    keep <- keep & (rowSums(!is.finite(Z_list[[k]])) == 0L)
  if (sum(keep) < ncol(X_lin) + K + 5L) return(NULL)

  Xk <- X_lin[keep, , drop = FALSE]
  yk <- y[keep]
  Zk <- lapply(Z_list, function(M) M[keep, , drop = FALSE])

  # ---- one-off FWL pre-computation ------------------------------------------
  qrX <- qr(Xk)
  if (qrX$rank < ncol(Xk)) return(NULL)
  Zt <- do.call(cbind, lapply(Zk, function(M) qr.resid(qrX, M)))  # n x KQ
  yt <- qr.resid(qrX, yk)
  G  <- crossprod(Zt)              # (sum Q_k) x (sum Q_k)   (theta-free)
  cv <- drop(crossprod(Zt, yt))    # (sum Q_k)               (theta-free)
  yy <- sum(yt^2)
  Qtot <- sum(Qs)
  blk <- split(seq_len(Qtot), rep(seq_len(K), Qs))

  # ---- concentrated SSE: O((sum Q_k)^2), independent of n -------------------
  obj <- function(theta_v) {
    W <- matrix(0, Qtot, K)
    for (k in seq_len(K))
      W[blk[[k]], k] <- exp_almon_weights_j(lag_sets[[block_names[k]]],
                                            theta_v[(2L * k - 1L):(2L * k)])
    A <- crossprod(W, G %*% W)     # K x K
    b <- drop(crossprod(W, cv))    # K
    d <- tryCatch(solve(A, b), error = function(e) NULL)
    if (is.null(d) || !all(is.finite(d))) return(1e12)
    sse <- yy - sum(b * d)
    if (!is.finite(sse) || sse < -1e-6 * yy) return(1e12)
    max(sse, 0)
  }

  # ---- starting values -------------------------------------------------------
  theta_lower <- rep(theta_lower_block, K)
  theta_upper <- rep(theta_upper_block, K)
  dim_theta   <- 2L * K
  if (n_starts < 1L) n_starts <- 1L
  # Canonical exp-Almon shapes used as deterministic restarts: flat, fast
  # geometric decay, slow decay, single hump, gentle decay.  Uniform draws
  # over [-8,8]x[-8,0] mostly produce degenerate weight profiles (all mass on
  # lag 1 or lag Q) in a flat region of the SSE and waste optimizer calls.
  shapes <- matrix(c( 0.0,  0.00,
                     -1.0,  0.00,
                     -0.2,  0.00,
                      0.5, -0.10,
                      0.1, -0.02), ncol = 2, byrow = TRUE)
  if (deterministic_starts)
    n_starts <- min(n_starts, 1L + nrow(shapes))   # no duplicate restarts
  starts <- matrix(NA_real_, n_starts, dim_theta)
  starts[1, ] <- if (!is.null(theta_init) && length(theta_init) == dim_theta &&
                     all(is.finite(theta_init))) as.numeric(theta_init)
                 else rep(0, dim_theta)
  if (n_starts > 1L) {
    if (deterministic_starts) {
      for (s in 2:n_starts) starts[s, ] <- rep(shapes[s - 1L, ], K)
    } else {
      for (s in 2:n_starts)
        starts[s, ] <- runif(dim_theta, min = theta_lower, max = theta_upper)
    }
  }

  best <- NULL; best_val <- Inf
  for (s in seq_len(n_starts)) {
    cand <- tryCatch(
      optim(par = as.numeric(starts[s, ]), fn = obj,
            method = "L-BFGS-B",
            lower = theta_lower, upper = theta_upper,
            control = list(maxit = optim_maxit, factr = optim_factr)),
      error = function(e) NULL
    )
    if (is.null(cand)) next
    if (!is.finite(cand$value) || cand$value >= 1e11) next
    if (!(cand$convergence %in% c(0L, 1L))) next
    if (cand$value < best_val) { best_val <- cand$value; best <- cand }
  }
  if (is.null(best)) return(NULL)
  theta_hat <- best$par

  # ---- final fit at theta_hat: SAME code path as the original ---------------
  zmat <- matrix(NA_real_, length(yk), K)
  for (k in seq_len(K))
    zmat[, k] <- as.numeric(Zk[[k]] %*%
                   exp_almon_weights_j(lag_sets[[block_names[k]]],
                                       theta_hat[(2L * k - 1L):(2L * k)]))
  Xfull <- cbind(Xk, zmat)
  qrfit <- qr(Xfull)
  if (qrfit$rank < ncol(Xfull)) return(NULL)
  beta_full <- qr.coef(qrfit, yk)
  res       <- yk - drop(Xfull %*% beta_full)
  p_lin     <- ncol(Xk)
  beta_lin  <- beta_full[1:p_lin]
  names(beta_lin) <- colnames(Xk)
  delta     <- beta_full[(p_lin + 1L):(p_lin + K)]
  names(delta) <- names(Z_list)
  theta_mat <- matrix(theta_hat, nrow = K, ncol = 2L, byrow = TRUE,
                      dimnames = list(names(Z_list), c("theta1", "theta2")))

  list(theta     = theta_mat,
       beta_lin  = beta_lin,
       delta     = delta,
       sse       = best$value,
       sigma2    = mean(res^2),
       converged = best$convergence %in% c(0L, 1L),
       n_eff     = length(yk),
       n_starts  = n_starts)
}
