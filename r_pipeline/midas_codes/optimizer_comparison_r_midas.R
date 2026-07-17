# =============================================================================
# optimizer_comparison_r_midas.R
# Diagnostic comparison for the existing R-MIDAS optimizer.
#
# PURPOSE
#   Compare the CURRENT estimator (warm start + deterministic canonical starts)
#   with an ADDITIONAL reproducible random-multistart diagnostic, without
#   changing or replacing fit_r_midas_fast().
#
# WORKS FOR
#   * Standard R-MIDAS: K = 1 Almon block (theta dimension = 2)
#   * Extended R-MIDAS: K >= 2 Almon blocks (theta dimension = 2K)
#
# OUTPUT PER FIT
#   best_sse, second_best_sse, sse_gap, theta estimates, convergence code,
#   number of distinct local solutions, number of successful starts, and a
#   comparison with the current deterministic fit.
#
# USAGE
#   source("common_utils.R")
#   source("fit_r_midas_fast.R")
#   source("optimizer_comparison_r_midas.R")
#
#   cmp <- compare_r_midas_optimizers(
#     y          = yf,
#     Z_list     = list(y_lags = Zf),
#     X_lin      = Xf,
#     theta_init = theta_warm,
#     lag_sets   = list(y_lags = c(1, 2, 3, 7)),
#     n_random   = 30L,
#     seed       = 20240501L + si * 100000L + oi
#   )
#   append_optimizer_diagnostic(cmp$summary,
#                               "results/optimizer_diagnostics_r_midas.csv")
#
# NOTE
#   This file intentionally does not overwrite fit_r_midas_fast().
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

# Cluster optimized theta vectors into numerically distinct solutions.
# Two solutions are treated as the same when both their SSE and parameter
# vectors are sufficiently close.  Tolerances can be adjusted by the caller.
.count_distinct_solutions <- function(runs,
                                      theta_tol = 1e-4,
                                      sse_rel_tol = 1e-8) {
  good <- Filter(function(z) isTRUE(z$success), runs)
  if (!length(good)) return(0L)

  ord <- order(vapply(good, `[[`, numeric(1), "sse"))
  good <- good[ord]
  reps <- list()

  for (g in good) {
    is_new <- TRUE
    if (length(reps)) {
      for (r in reps) {
        sse_scale <- max(1, abs(r$sse), abs(g$sse))
        close_sse <- abs(g$sse - r$sse) <= sse_rel_tol * sse_scale
        close_th  <- max(abs(g$theta - r$theta)) <= theta_tol
        if (close_sse && close_th) {
          is_new <- FALSE
          break
        }
      }
    }
    if (is_new) reps[[length(reps) + 1L]] <- g
  }
  length(reps)
}

# Build the same FWL concentrated objective used by fit_r_midas_fast().
# Kept local to this diagnostic so the production estimator remains untouched.
.make_fwl_objective <- function(y, Z_list, X_lin, lag_sets = NULL) {
  K <- length(Z_list)
  if (K < 1L) stop("Need at least one Almon block.")
  block_names <- names(Z_list)
  if (is.null(block_names) || any(!nzchar(block_names)))
    stop("Z_list must be a named list.")

  lag_sets <- lag_sets %||% list()
  for (nm in block_names) {
    if (is.null(lag_sets[[nm]])) lag_sets[[nm]] <- seq_len(ncol(Z_list[[nm]]))
    if (length(lag_sets[[nm]]) != ncol(Z_list[[nm]]))
      stop("lag_sets[[", nm, "]] has the wrong length.")
  }

  n <- length(y)
  if (is.null(X_lin))
    X_lin <- matrix(1, n, 1L, dimnames = list(NULL, "intercept"))

  keep <- stats::complete.cases(X_lin) & is.finite(y)
  for (k in seq_len(K))
    keep <- keep & (rowSums(!is.finite(Z_list[[k]])) == 0L)

  if (sum(keep) < ncol(X_lin) + K + 5L)
    return(NULL)

  Xk <- X_lin[keep, , drop = FALSE]
  yk <- y[keep]
  Zk <- lapply(Z_list, function(M) M[keep, , drop = FALSE])

  qrX <- qr(Xk)
  if (qrX$rank < ncol(Xk)) return(NULL)

  Zt <- do.call(cbind, lapply(Zk, function(M) qr.resid(qrX, M)))
  yt <- qr.resid(qrX, yk)
  G  <- crossprod(Zt)
  cv <- drop(crossprod(Zt, yt))
  yy <- sum(yt^2)

  Qs   <- vapply(Z_list, ncol, integer(1))
  Qtot <- sum(Qs)
  blk  <- split(seq_len(Qtot), rep(seq_len(K), Qs))

  obj <- function(theta_v) {
    W <- matrix(0, Qtot, K)
    for (k in seq_len(K)) {
      th <- theta_v[(2L * k - 1L):(2L * k)]
      W[blk[[k]], k] <- exp_almon_weights_j(lag_sets[[block_names[k]]], th)
    }
    A <- crossprod(W, G %*% W)
    b <- drop(crossprod(W, cv))
    d <- tryCatch(solve(A, b), error = function(e) NULL)
    if (is.null(d) || !all(is.finite(d))) return(1e12)
    sse <- yy - sum(b * d)
    if (!is.finite(sse) || sse < -1e-6 * max(yy, 1)) return(1e12)
    max(sse, 0)
  }

  list(obj = obj, K = K, block_names = block_names, n_eff = sum(keep))
}

.run_one_start <- function(start, obj, lower, upper,
                           optim_maxit = 2000L,
                           optim_factr = 5e8,
                           start_type = "random",
                           start_id = NA_integer_) {
  ans <- tryCatch(
    optim(par = as.numeric(start), fn = obj, method = "L-BFGS-B",
          lower = lower, upper = upper,
          control = list(maxit = optim_maxit, factr = optim_factr)),
    error = function(e) NULL
  )

  success <- !is.null(ans) && is.finite(ans$value) && ans$value < 1e11 &&
    ans$convergence %in% c(0L, 1L)

  list(
    start_id = start_id,
    start_type = start_type,
    start = as.numeric(start),
    theta = if (success) as.numeric(ans$par) else rep(NA_real_, length(start)),
    sse = if (success) as.numeric(ans$value) else NA_real_,
    convergence = if (is.null(ans)) NA_integer_ else as.integer(ans$convergence),
    message = if (is.null(ans)) "optim error" else as.character(ans$message %||% ""),
    success = success
  )
}

compare_r_midas_optimizers <- function(
    y, Z_list, X_lin,
    theta_init = NULL,
    lag_sets = NULL,
    n_current_starts = 5L,
    n_random = 30L,
    seed = 20240501L,
    theta_lower_block = c(-8, -8),
    theta_upper_block = c( 8,  0),
    optim_maxit = 2000L,
    optim_factr = 5e8,
    theta_tol = 1e-4,
    sse_rel_tol = 1e-8,
    metadata = list()) {

  prep <- .make_fwl_objective(y, Z_list, X_lin, lag_sets)
  if (is.null(prep)) return(NULL)

  K <- prep$K
  dim_theta <- 2L * K
  lower <- rep(theta_lower_block, K)
  upper <- rep(theta_upper_block, K)

  # Existing production implementation: unchanged.
  current_fit <- fit_r_midas_fast(
    y = y, Z_list = Z_list, X_lin = X_lin,
    theta_init = theta_init,
    theta_lower_block = theta_lower_block,
    theta_upper_block = theta_upper_block,
    optim_maxit = optim_maxit,
    optim_factr = optim_factr,
    n_starts = n_current_starts,
    deterministic_starts = TRUE,
    lag_sets = lag_sets
  )

  # Diagnostic random protocol: warm start (or zero) plus n_random independent
  # reproducible starts.  For K > 1 every block receives its own random pair,
  # unlike the canonical protocol that repeats the same shape across blocks.
  set.seed(seed)
  start0 <- if (!is.null(theta_init) && length(theta_init) == dim_theta &&
                all(is.finite(theta_init))) as.numeric(theta_init) else rep(0, dim_theta)

  runs <- list(.run_one_start(
    start0, prep$obj, lower, upper,
    optim_maxit, optim_factr,
    start_type = if (all(start0 == 0)) "zero" else "warm",
    start_id = 0L
  ))

  if (n_random > 0L) {
    for (s in seq_len(n_random)) {
      st <- stats::runif(dim_theta, min = lower, max = upper)
      runs[[length(runs) + 1L]] <- .run_one_start(
        st, prep$obj, lower, upper,
        optim_maxit, optim_factr,
        start_type = "random", start_id = s
      )
    }
  }

  good <- Filter(function(z) isTRUE(z$success), runs)
  if (!length(good)) {
    random_best <- NULL
    best_sse <- second_best_sse <- NA_real_
    convergence_code <- NA_integer_
    theta_best <- rep(NA_real_, dim_theta)
  } else {
    ord <- order(vapply(good, `[[`, numeric(1), "sse"))
    good <- good[ord]
    random_best <- good[[1L]]
    best_sse <- random_best$sse
    second_best_sse <- if (length(good) >= 2L) good[[2L]]$sse else NA_real_
    convergence_code <- random_best$convergence
    theta_best <- random_best$theta
  }

  n_distinct <- .count_distinct_solutions(runs, theta_tol, sse_rel_tol)
  current_sse <- if (is.null(current_fit)) NA_real_ else current_fit$sse
  current_theta <- if (is.null(current_fit)) rep(NA_real_, dim_theta) else
    as.numeric(t(current_fit$theta))

  # as.numeric(t(theta matrix)) gives block-wise theta1,theta2 ordering.
  theta_names <- unlist(lapply(prep$block_names,
                               function(nm) paste0(c("theta1_", "theta2_"), nm)))

  summary <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!ncol(summary)) summary <- data.frame(row_id = 1L)[, FALSE, drop = FALSE]
  summary$n_almon_blocks <- K
  summary$theta_dimension <- dim_theta
  summary$n_effective <- prep$n_eff
  summary$seed <- seed
  summary$n_random_requested <- n_random
  summary$n_successful_random_protocol <- length(good)
  summary$n_distinct_local_solutions <- n_distinct
  summary$best_sse <- best_sse
  summary$second_best_sse <- second_best_sse
  summary$sse_gap <- second_best_sse - best_sse
  summary$convergence_code <- convergence_code
  summary$current_deterministic_sse <- current_sse
  summary$random_minus_current_sse <- best_sse - current_sse
  summary$random_better_than_current <- is.finite(best_sse) && is.finite(current_sse) &&
    best_sse < current_sse - sse_rel_tol * max(1, abs(current_sse))

  for (j in seq_len(dim_theta)) {
    summary[[paste0("best_", theta_names[j])]] <- theta_best[j]
    summary[[paste0("current_", theta_names[j])]] <- current_theta[j]
  }

  run_table <- do.call(rbind, lapply(runs, function(z) {
    row <- data.frame(
      start_id = z$start_id,
      start_type = z$start_type,
      success = z$success,
      sse = z$sse,
      convergence_code = z$convergence,
      message = z$message,
      stringsAsFactors = FALSE
    )
    for (j in seq_len(dim_theta)) {
      row[[paste0("start_", theta_names[j])]] <- z$start[j]
      row[[paste0("hat_", theta_names[j])]] <- z$theta[j]
    }
    row
  }))

  list(summary = summary,
       runs = run_table,
       current_fit = current_fit,
       random_best = random_best)
}

append_optimizer_diagnostic <- function(x, file) {
  if (is.null(x) || !nrow(x)) return(invisible(FALSE))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.table(x, file = file, sep = ",", row.names = FALSE,
              col.names = !file.exists(file), append = file.exists(file),
              quote = TRUE, na = "")
  invisible(TRUE)
}
