# =============================================================================
# Simulation code for:
#
#   Wu, L. (2025). Restricted Mean Survival Time Under Non-Proportional
#   Hazards: How Much Power Can Covariate Adjustment Recover?
#   Submitted to Statistics in Medicine.
#
# Author : Longyang Wu  <Longyang.wu@gmail.com>
# Repo   : https://github.com/longyangw/rmst-nph-adjustment
#
# Methods compared
#   LR            : ordinary log-rank
#   FH(0,1)       : Fleming-Harrington weighted log-rank, weight on late events
#   FH(1,1)       : Fleming-Harrington weighted log-rank, weight on mid events
#   KM-RMST       : Kaplan-Meier RMST with Greenwood-type variance
#   Plug-in Cox   : retained as a diagnostic comparator only (T1E inflation)
#   PO-AIPW+Cox   : jackknife pseudo-observation augmentation with pooled Cox
#
# Data-generating mechanisms (see Table 1 of the paper)
#   DGM0  Null                       HR ≡ 1.0
#   DGM1  Delayed benefit            HR 1.0 / 0.55, changepoint t = 4
#   DGM2  Crossing hazards           HR 1.3 / 0.50, changepoint t = 4
#   DGM3  Early harm, later benefit  HR 1.8 / 0.45, changepoint t = 3
#   DGM4  Gamma frailty              theta = 0.5, conditional HR_trt = 0.70
#   DGM5  Covariate-dependent cens.  DGM1 event process; dropout ∝ exp(0.4·η)
#   DGM6  Treatment switching        DGM1 + 30% control rescue at t ≥ 6m
#                                    (treatment-policy estimand)
#
# True Δ(τ) values are computed at runtime via Monte Carlo (n = 500,000 per DGM)
# under the actual covariate distribution.
#
# Design
#   N ∈ {250, 500, 1000}
#   τ = 24 months
#   α = 0.025 (one-sided)
#   B = 2000 Monte Carlo replications per (DGM, N) cell
#   12-month uniform accrual; minimum follow-up 12 months
#   Covariates: X1 ~ N(0,1), X2 ~ Bernoulli(0.5)
#   Cox log-hazard coefficients: β_X1 = 0.3, β_X2 = 0.2
#
# Required packages: survival, pseudo, parallel
# Approximate runtime: 2-4 hours on 8 cores for the full grid.
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(pseudo)
  library(parallel)
})

# ── 0. Global parameters ──────────────────────────────────────────────────────

TAU       <- 24          # restriction time (months)
ALPHA     <- 0.025       # one-sided significance level
ACCRUAL   <- 12          # uniform accrual window (months)
LAM0      <- log(2)/12   # control-arm baseline hazard (median = 12 months)
BETA_X1   <- 0.3         # Cox log-hazard coefficient for X1
BETA_X2   <- 0.2         # Cox log-hazard coefficient for X2
B         <- 2000        # Monte Carlo replications per (DGM, N) cell
N_VEC     <- c(250, 500, 1000)
DGM_VEC   <- 0:6
SEED_BASE <- 20250101
RUN_PLUG_COX <- TRUE     # diagnostic comparator; excluded from main power table

# ── 1. Cox-correct event-time generator ──────────────────────────────────────
#
# For a Cox PH model  λ(t|X) = λ_0(t) exp(β'X), generate T by inverting the
# cumulative hazard:  H_0(T) = V,  where V = -log(U) / exp(β'X), U ~ U(0,1).
# This keeps changepoints at fixed calendar times for all subjects.
#
# For piecewise-exponential baseline with K pieces:
#   λ_0(t) = lam_k  for  t ∈ [cp_{k-1}, cp_k)
#   H_0(t) = sum_{j<k} lam_j*(cp_j - cp_{j-1}) + lam_k*(t - cp_{k-1})

H0_inv <- function(v, lams, cps) {
  cum_H <- cumsum(lams * diff(c(0, cps)))
  k <- findInterval(v, c(0, cum_H))
  k <- pmax(pmin(k, length(lams)), 1L)
  prev_H  <- c(0, cum_H)[k]
  prev_cp <- c(0, cps)[k]
  t_raw   <- prev_cp + (v - prev_H) / lams[k]
  pmax(t_raw, prev_cp)
}

gen_cox_time <- function(lams, cps, eta) {
  U <- runif(1)
  V <- -log(U) / exp(eta)
  H0_inv(V, lams, cps)
}

# ── 2. Dataset generator ──────────────────────────────────────────────────────

sim_dataset <- function(n, dgm, tau = TAU, accrual = ACCRUAL) {

  arm <- rep(0:1, each = n/2)
  X1  <- rnorm(n)
  X2  <- rbinom(n, 1, 0.5)
  eta <- BETA_X1 * X1 + BETA_X2 * X2

  # Baseline hazard specifications per DGM
  lams0 <- rep(LAM0, 1); cps0 <- Inf

  if (dgm %in% c(0, 1, 5, 6)) {
    if (dgm == 0) {
      lams1 <- LAM0; cps1 <- Inf
    } else {
      lams1 <- c(LAM0, 0.55*LAM0); cps1 <- c(4, Inf)
    }
  } else if (dgm == 2) {
    lams1 <- c(1.3*LAM0, 0.50*LAM0); cps1 <- c(4, Inf)
  } else if (dgm == 3) {
    lams1 <- c(1.8*LAM0, 0.45*LAM0); cps1 <- c(3, Inf)
  } else if (dgm == 4) {
    lams1 <- 0.70*LAM0; cps1 <- Inf
  }

  # Generate event times
  T_true <- numeric(n)

  if (dgm == 4) {
    # Gamma frailty with conditional proportional hazards
    theta   <- 0.5
    frailty <- rgamma(n, shape = 1/theta, rate = 1/theta)
    for (i in seq_len(n)) {
      lam_i <- (if (arm[i] == 0) LAM0 else 0.70*LAM0) * frailty[i]
      T_true[i] <- gen_cox_time(lam_i, Inf, eta[i])
    }
  } else {
    for (i in seq_len(n)) {
      lams_i <- if (arm[i] == 0) lams0 else lams1
      cps_i  <- if (arm[i] == 0) cps0  else cps1
      T_true[i] <- gen_cox_time(lams_i, cps_i, eta[i])
    }
  }

  # Treatment switching (DGM6): 30% of control with T > 6 switch at month 6
  switch_flag <- rep(FALSE, n)
  if (dgm == 6) {
    ctrl_idx  <- which(arm == 0)
    eligible  <- ctrl_idx[T_true[ctrl_idx] > 6]
    switchers <- eligible[runif(length(eligible)) < 0.30]
    switch_flag[switchers] <- TRUE
    for (i in switchers) {
      T_true[i] <- 6 + gen_cox_time(0.55*LAM0, Inf, eta[i])
    }
  }

  # Censoring
  entry   <- runif(n, 0, accrual)
  C_admin <- tau - entry

  if (dgm == 5) {
    # Covariate-dependent dropout
    C_drop <- rexp(n, rate = 0.1 * exp(0.4 * eta))
  } else {
    # ~10% uniform dropout over [0, tau]
    C_drop <- rexp(n, rate = 0.10 / tau)
  }

  C_time <- pmin(C_admin, C_drop)

  # Observed time and event indicator
  time  <- pmin(T_true, C_time, tau)
  event <- as.integer(T_true <= C_time & T_true <= tau)
  time  <- pmax(time, 1e-6)

  data.frame(
    id          = seq_len(n),
    arm         = arm,
    X1          = X1,
    X2          = X2,
    eta         = eta,
    T_true      = T_true,
    C_time      = C_time,
    time        = time,
    event       = event,
    entry       = entry,
    switch_flag = switch_flag
  )
}


# ── 3. KM-RMST with Greenwood variance ───────────────────────────────────────

km_rmst <- function(time, event, tau) {

  km <- survfit(Surv(time, event) ~ 1)

  t_ev <- km$time[km$time <= tau]
  d_ev <- km$n.event[km$time <= tau]
  y_ev <- km$n.risk[km$time <= tau]
  s_ev <- km$surv[km$time <= tau]

  # RMST = integral_0^tau S(t) dt via the KM step function
  breaks <- c(0, t_ev, tau)
  s_step <- c(1, s_ev)
  if (length(s_step) < length(breaks))
    s_step <- c(s_step, tail(s_step, 1))
  rmst <- sum(s_step[-length(s_step)] * diff(breaks))

  # Greenwood variance: Var(RMST) = sum_k A_k^2 * d_k / (y_k*(y_k - d_k))
  # A_k = int_{t_k}^{tau} S(t) dt
  K <- length(t_ev)
  if (K == 0) return(list(rmst = rmst, var = 0))

  A_k <- numeric(K)
  for (k in seq_len(K)) {
    t_bounds <- c(t_ev[k:K], tau)
    s_vals   <- s_ev[k:K]
    widths   <- diff(t_bounds)
    A_k[k]  <- sum(s_vals * widths)
  }

  denom  <- y_ev * (y_ev - d_ev)
  safe   <- denom > 0
  var_rmst <- sum(A_k[safe]^2 * d_ev[safe] / denom[safe])

  list(rmst = rmst, var = var_rmst)
}


# ── 4. Test statistics ────────────────────────────────────────────────────────

# 4a. Log-rank
lr_test <- function(dat) {
  fit <- survdiff(Surv(time, event) ~ arm, data = dat)
  O <- fit$obs; E <- fit$exp
  sign(E[2] - O[2]) * sqrt(fit$chisq)
}

# 4b. Fleming-Harrington weighted log-rank
fh_test <- function(dat, rho, gamma) {
  km_pool  <- survfit(Surv(time, event) ~ 1, data = dat)
  s_pool   <- stepfun(km_pool$time, c(1, km_pool$surv))
  ev_times <- sort(unique(dat$time[dat$event == 1]))

  num <- 0; den_sq <- 0
  for (tj in ev_times) {
    at_risk <- dat$time >= tj
    n1j <- sum(at_risk & dat$arm == 1)
    n0j <- sum(at_risk & dat$arm == 0)
    nj  <- n1j + n0j
    if (nj < 2) next

    d1j <- sum(dat$time == tj & dat$event == 1 & dat$arm == 1)
    dj  <- sum(dat$time == tj & dat$event == 1)

    S_tj <- s_pool(tj - 1e-8)
    wj   <- S_tj^rho * (1 - S_tj)^gamma

    ej1    <- n1j * dj / nj
    vj_raw <- n1j * n0j * dj * (nj - dj) / (nj^2 * (nj - 1))
    vj     <- if (!is.finite(vj_raw)) 0 else vj_raw

    num    <- num    + wj * (ej1 - d1j)
    den_sq <- den_sq + wj^2 * vj
  }
  if (den_sq <= 0) return(NA_real_)
  num / sqrt(den_sq)
}

# 4c. KM-RMST test
km_rmst_test <- function(dat, tau = TAU) {
  d1 <- dat[dat$arm == 1, ]; d0 <- dat[dat$arm == 0, ]
  r1 <- km_rmst(d1$time, d1$event, tau)
  r0 <- km_rmst(d0$time, d0$event, tau)
  delta <- r1$rmst - r0$rmst
  se    <- sqrt(r1$var + r0$var)
  list(z = delta / se, delta = delta, se = se)
}

# 4d. Plug-in Cox (diagnostic only; severe T1E inflation)
plug_cox_test <- function(dat, tau = TAU) {
  tryCatch({
    cox_fit <- coxph(Surv(time, event) ~ arm + X1 + X2, data = dat)

    pred_rmst_arm <- function(arm_val) {
      nd   <- dat; nd$arm <- arm_val
      sf   <- survfit(cox_fit, newdata = nd)
      apply(sf$surv, 2, function(sv) {
        tt <- sf$time
        if (max(tt) < tau) { tt <- c(tt, tau); sv <- c(sv, tail(sv, 1)) }
        keep <- tt <= tau
        tt   <- c(0, tt[keep]); sv <- c(1, sv[keep])
        sum(diff(tt) * head(sv, -1))
      })
    }

    r1 <- pred_rmst_arm(1); r0 <- pred_rmst_arm(0)
    delta <- mean(r1) - mean(r0)
    se    <- sd(r1 - r0) / sqrt(nrow(dat))
    list(z = delta / se, delta = delta, se = se)
  }, error = function(e) list(z = NA_real_, delta = NA_real_, se = NA_real_))
}


# 4e. PO-AIPW+Cox: jackknife pseudo-observations + pooled Cox working model
augmentation_test <- function(dat, tau = TAU) {
  tryCatch({
    n <- nrow(dat)

    # Step 1: jackknife pseudo-observations for RMST
    theta <- pseudomean(dat$time, dat$event, tmax = tau)

    # Step 2: pooled Cox working model
    cox_fit <- coxph(Surv(time, event) ~ arm + X1 + X2, data = dat)

    # Step 3: predict arm-specific RMST for each subject under arm = 0 and arm = 1
    pred_rmst_cox <- function(arm_val) {
      nd <- dat; nd$arm <- arm_val
      sf <- survfit(cox_fit, newdata = nd)
      apply(sf$surv, 2, function(sv) {
        tt <- sf$time
        if (max(tt) < tau) { tt <- c(tt, tau); sv <- c(sv, tail(sv, 1)) }
        keep <- tt <= tau
        tt   <- c(0, tt[keep]); sv <- c(1, sv[keep])
        sum(diff(tt) * head(sv, -1))
      })
    }
    mu1 <- pred_rmst_cox(1)
    mu0 <- pred_rmst_cox(0)

    # Step 4: AIPW influence functions with known propensity pi = 0.5
    arm_i   <- dat$arm
    phi1_if <- 2*arm_i*(theta - mu1) + mu1
    phi0_if <- 2*(1-arm_i)*(theta - mu0) + mu0

    delta_aug <- mean(phi1_if) - mean(phi0_if)
    psi       <- (phi1_if - phi0_if) - delta_aug
    se_aug    <- sqrt(mean(psi^2) / n)
    z_aug     <- delta_aug / se_aug

    list(z = z_aug, delta = delta_aug, se = se_aug)
  }, error = function(e)
    list(z = NA_real_, delta = NA_real_, se = NA_real_))
}


# ── 5. True RMST contrasts via Monte Carlo ───────────────────────────────────
#
# Computed at runtime under the actual covariate distribution
# (n = 500,000 per DGM). DGM6 uses the treatment-policy estimand.

compute_true_delta <- function(dgm, tau = TAU, n_mc = 500000,
                                seed_mc = 99999999) {
  set.seed(seed_mc + dgm)
  X1  <- rnorm(n_mc)
  X2  <- rbinom(n_mc, 1, 0.5)
  eta <- BETA_X1*X1 + BETA_X2*X2

  # Control arm: simple exponential
  V0 <- rexp(n_mc, 1)
  T0 <- V0 / (LAM0 * exp(eta))
  rmst0 <- mean(pmin(T0, tau))

  if (dgm == 0) return(0)

  if (dgm %in% c(1, 5)) {
    lams1 <- c(LAM0, 0.55*LAM0); cps1 <- c(4, Inf)
  } else if (dgm == 2) {
    lams1 <- c(1.3*LAM0, 0.50*LAM0); cps1 <- c(4, Inf)
  } else if (dgm == 3) {
    lams1 <- c(1.8*LAM0, 0.45*LAM0); cps1 <- c(3, Inf)
  } else if (dgm == 4) {
    # Frailty must be drawn once per subject and applied to BOTH arms.
    theta <- 0.5
    fr    <- rgamma(n_mc, 1/theta, 1/theta)
    T0_4  <- rexp(n_mc, 1) / (LAM0      * fr * exp(eta))
    T1_4  <- rexp(n_mc, 1) / (0.70*LAM0 * fr * exp(eta))
    return(mean(pmin(T1_4, tau)) - mean(pmin(T0_4, tau)))
  } else if (dgm == 6) {
    lams1 <- c(LAM0, 0.55*LAM0); cps1 <- c(4, Inf)
  }

  # Invert piecewise-exp Cox for treatment
  V1  <- rexp(n_mc, 1) / exp(eta)
  H_cp <- lams1[1] * cps1[1]
  if (length(lams1) == 1) {
    T1 <- V1 / lams1[1]
  } else {
    T1 <- ifelse(V1 <= H_cp,
                 V1 / lams1[1],
                 cps1[1] + (V1 - H_cp) / lams1[2])
  }

  # DGM6: treatment-policy estimand with 30% control switching at t = 6
  if (dgm == 6) {
    switchers <- which(T0 > 6 & runif(n_mc) < 0.30)
    V_sw <- rexp(length(switchers), 1) / exp(eta[switchers])
    T0[switchers] <- 6 + V_sw / (0.55*LAM0)
    rmst0 <- mean(pmin(T0, tau))
  }

  mean(pmin(T1, tau)) - rmst0
}

cat("Computing TRUE_DELTA via Monte Carlo (n = 500,000 per DGM)...\n")
TRUE_DELTA <- setNames(
  vapply(0:6, function(d) compute_true_delta(d), numeric(1)),
  as.character(0:6)
)
cat("TRUE_DELTA:\n")
print(round(TRUE_DELTA, 4))
cat("\n")


# ── 6. Single replication ─────────────────────────────────────────────────────

one_rep <- function(n, dgm, tau = TAU, alpha = ALPHA) {
  dat    <- sim_dataset(n, dgm, tau)
  z_crit <- qnorm(1 - alpha)
  reject <- function(z) as.integer(!is.na(z) && z > z_crit)

  z_lr   <- lr_test(dat)
  z_fh01 <- fh_test(dat, rho = 0, gamma = 1)
  z_fh11 <- fh_test(dat, rho = 1, gamma = 1)

  km_res  <- km_rmst_test(dat, tau)
  aug_res <- augmentation_test(dat, tau)
  z_cox   <- if (RUN_PLUG_COX) plug_cox_test(dat, tau)$z else NA_real_

  true_d <- TRUE_DELTA[as.character(dgm)]

  ci_cover <- function(delta_val, se_val) {
    if (is.na(delta_val) || is.na(se_val)) return(NA_integer_)
    lo <- delta_val - qnorm(0.975) * se_val
    hi <- delta_val + qnorm(0.975) * se_val
    as.integer(lo <= true_d && true_d <= hi)
  }

  data.frame(
    rej_lr   = reject(z_lr),
    rej_fh01 = reject(z_fh01),
    rej_fh11 = reject(z_fh11),
    rej_km   = reject(km_res$z),
    rej_cox  = reject(z_cox),
    rej_aug  = reject(aug_res$z),

    delta_km  = km_res$delta,
    delta_aug = aug_res$delta,
    se_km     = km_res$se,
    se_aug    = aug_res$se,

    cov_km  = ci_cover(km_res$delta,  km_res$se),
    cov_aug = ci_cover(aug_res$delta, aug_res$se),

    true_delta = as.numeric(true_d)
  )
}


# ── 7. Main simulation loop ───────────────────────────────────────────────────

run_simulation <- function(n_vec   = N_VEC,
                           dgm_vec = DGM_VEC,
                           B       = B,
                           ncores  = max(1L, detectCores() - 1L),
                           seed    = SEED_BASE) {
  results_list <- list()

  for (dgm in dgm_vec) {
    for (n in n_vec) {
      cat(sprintf("\n=== DGM%d  N=%4d  B=%d  cores=%d ===\n",
                  dgm, n, B, ncores))
      t0 <- proc.time()

      if (ncores > 1) {
        cl <- makeCluster(ncores)
        on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
        needed <- c("sim_dataset","one_rep","km_rmst","km_rmst_test",
                    "lr_test","fh_test","plug_cox_test",
                    "augmentation_test","H0_inv","gen_cox_time",
                    "TAU","ALPHA","ACCRUAL","LAM0","BETA_X1","BETA_X2",
                    "TRUE_DELTA","RUN_PLUG_COX","SEED_BASE")
        clusterExport(cl, varlist = needed, envir = .GlobalEnv)
        clusterEvalQ(cl, { library(survival); library(pseudo) })
        reps <- parLapply(cl, seq_len(B), function(b) {
          set.seed(seed + dgm * 10000 + n * 10 + b)
          tryCatch(one_rep(n, dgm), error = function(e) NULL)
        })
        stopCluster(cl)
      } else {
        reps <- lapply(seq_len(B), function(b) {
          set.seed(seed + dgm * 10000 + n * 10 + b)
          tryCatch(one_rep(n, dgm), error = function(e) NULL)
        })
      }

      reps  <- Filter(Negate(is.null), reps)
      n_ok  <- length(reps)
      skip  <- B - n_ok
      res   <- do.call(rbind, reps)

      mu <- function(x) mean(x, na.rm = TRUE)

      smry <- data.frame(
        dgm  = dgm, n = n, B = B, n_ok = n_ok, skip = skip,
        power_lr   = mu(res$rej_lr),
        power_fh01 = mu(res$rej_fh01),
        power_fh11 = mu(res$rej_fh11),
        power_km   = mu(res$rej_km),
        power_cox  = mu(res$rej_cox),
        power_aug  = mu(res$rej_aug),
        bias_km    = mu(res$delta_km)  - res$true_delta[1],
        bias_aug   = mu(res$delta_aug) - res$true_delta[1],
        var_km     = var(res$delta_km,  na.rm = TRUE),
        var_aug    = var(res$delta_aug, na.rm = TRUE),
        cov_km     = mu(res$cov_km),
        cov_aug    = mu(res$cov_aug),
        true_delta = res$true_delta[1]
      )
      smry$re_aug  <- smry$var_aug / smry$var_km
      smry$prr_km  <- smry$power_km  / smry$power_lr
      smry$prr_aug <- smry$power_aug / smry$power_lr

      elapsed <- proc.time() - t0
      cat(sprintf("  Done: %.1f s  skip=%d\n", elapsed["elapsed"], skip))
      cat(sprintf("  LR=%.1f  FH01=%.1f  FH11=%.1f  KM=%.1f  AUG=%.1f  RE_aug=%.3f\n",
                  smry$power_lr*100, smry$power_fh01*100, smry$power_fh11*100,
                  smry$power_km*100, smry$power_aug*100, smry$re_aug))

      results_list[[paste0("DGM", dgm, "_N", n)]] <- list(summary = smry, reps = res)
    }
  }
  results_list
}


# ── 8. Summary tables ─────────────────────────────────────────────────────────

print_power_table <- function(results, n_target = 500) {
  dgms    <- 0:6
  methods <- c("lr","fh01","fh11","km","cox","aug")
  labels  <- c("LR","WLR FH(0,1)","WLR FH(1,1)","KM-RMST",
                "Plug-in Cox*","PO-AIPW+Cox")

  cat(sprintf("\n=== Power (%%) and PRR at N=%d ===\n", n_target))
  cat(sprintf("  * Plug-in Cox retained for diagnostics; severe T1E inflation\n\n"))
  cat(sprintf("%-20s", "Method"))
  for (d in dgms) cat(sprintf("  DGM%-4d", d))
  cat("\n", strrep("-", 20 + 8*length(dgms)), "\n", sep="")

  for (j in seq_along(methods)) {
    m <- methods[j]
    cat(sprintf("%-20s", labels[j]))
    for (d in dgms) {
      key  <- paste0("DGM", d, "_N", n_target)
      if (!key %in% names(results)) { cat("      NA"); next }
      s    <- results[[key]]$summary
      pow  <- s[[paste0("power_", m)]] * 100
      if (m == "lr" || d == 0) {
        cat(sprintf("  %6.1f", pow))
      } else {
        prr <- s[[paste0("power_", m)]] / s$power_lr
        cat(sprintf("  %5.1f[%.2f]", pow, prr))
      }
    }
    cat("\n")
  }
}

print_bias_table <- function(results, n_target = 500) {
  cat(sprintf("\n=== Bias (months) and 95%% Coverage at N=%d ===\n\n", n_target))
  cat(sprintf("%-8s  %-22s  %-22s\n", "DGM", "KM-RMST", "PO-AIPW+Cox"))
  for (d in 0:6) {
    key <- paste0("DGM", d, "_N", n_target)
    if (!key %in% names(results)) next
    s <- results[[key]]$summary
    cat(sprintf("DGM%-5d  %+.3f (%.3f)         %+.3f (%.3f)\n",
                d,
                s$bias_km,  s$cov_km,
                s$bias_aug, s$cov_aug))
  }
}

print_var_table <- function(results, n_target = 500) {
  cat(sprintf("\n=== Variance and RE at N=%d  (RE<1 = better than KM) ===\n\n",
              n_target))
  cat(sprintf("%-8s  %8s  %8s  %7s\n",
              "DGM","Var(KM)","Var(Aug)","RE(Aug)"))
  for (d in 0:6) {
    key <- paste0("DGM", d, "_N", n_target)
    if (!key %in% names(results)) next
    s <- results[[key]]$summary
    cat(sprintf("DGM%-5d  %8.4f  %8.4f  %7.3f\n",
                d, s$var_km, s$var_aug, s$re_aug))
  }
}


# ── 9. Entry point ────────────────────────────────────────────────────────────

main <- function() {
  cat("============================================================\n")
  cat("  RMST Covariate Adjustment Simulation\n")
  cat(sprintf("  B=%d | tau=%g m | alpha=%g | seed=%d\n",
              B, TAU, ALPHA, SEED_BASE))
  cat(sprintf("  DGMs: %s\n", paste(DGM_VEC, collapse=", ")))
  cat(sprintf("  N   : %s\n", paste(N_VEC,   collapse=", ")))
  cat("============================================================\n")

  results <- run_simulation(
    n_vec   = N_VEC,
    dgm_vec = DGM_VEC,
    B       = B,
    ncores  = max(1L, detectCores() - 1L),
    seed    = SEED_BASE
  )

  outfile <- "RMST_simulation_results.RData"
  save(results, file = outfile)
  cat("\nResults saved to", outfile, "\n")

  print_power_table(results, 500)
  print_bias_table(results,  500)
  print_var_table(results,   500)

  invisible(results)
}

if (!interactive()) main()
