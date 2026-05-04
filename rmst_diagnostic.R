#' RMST Design-Stage Diagnostic Tool
#'
#' Computes C(tau), SRR, t_bar_delta and provides a recommendation
#' on whether to use RMST or log-rank based on the PH reference range.
#'
#' Based on: Wu (2025) "Restricted Mean Survival Time Under
#' Non-Proportional Hazards: How Much Power Can Covariate Adjustment Recover?"
#'
#' Usage:
#'   rmst_diagnostic(type = "exponential", lambda0, HR, tau)
#'   rmst_diagnostic(type = "delayed",     lambda0, HR, tau, delay)
#'   rmst_diagnostic(type = "crossing",    lambda0, HR_early, HR_late, tau, changepoint)
#'   rmst_diagnostic(type = "frailty",     lambda0, HR_cond, tau, theta)

rmst_diagnostic <- function(
    type       = c("exponential", "delayed", "crossing", "frailty"),
    lambda0    = NULL,   # control arm hazard (or specify median0)
    median0    = NULL,   # control arm median survival (months)
    HR         = NULL,   # overall or post-onset HR
    HR_early   = NULL,   # HR before changepoint (crossing hazards)
    HR_late    = NULL,   # HR after changepoint  (crossing hazards)
    HR_cond    = NULL,   # conditional HR (frailty)
    tau        = 24,     # restriction time (months)
    delay      = NULL,   # onset delay d (months), for type="delayed"
    changepoint= NULL,   # changepoint (months), for type="crossing"
    theta      = 0.5,    # frailty variance, for type="frailty"
    n_grid     = 500,    # number of integration grid points
    verbose    = TRUE
) {
  type <- match.arg(type)

  # --- resolve lambda0 ---
  if (is.null(lambda0) && !is.null(median0)) {
    lambda0 <- log(2) / median0
  }
  if (is.null(lambda0)) stop("Specify either lambda0 or median0.")

  # ----------------------------------------------------------------
  # Define S0, S1, f0, f1 by type
  # ----------------------------------------------------------------

  if (type == "exponential") {
    if (is.null(HR)) stop("Specify HR for type='exponential'.")
    lambda1 <- HR * lambda0
    S0 <- function(t) exp(-lambda0 * t)
    S1 <- function(t) exp(-lambda1 * t)
    f0 <- function(t) lambda0 * exp(-lambda0 * t)
    f1 <- function(t) lambda1 * exp(-lambda1 * t)
    label <- sprintf("Exponential PH  (HR=%.2f)", HR)

  } else if (type == "delayed") {
    if (is.null(HR) || is.null(delay))
      stop("Specify HR and delay for type='delayed'.")
    d <- delay
    lambda1_post <- HR * lambda0
    S0 <- function(t) exp(-lambda0 * t)
    S1 <- function(t) {
      ifelse(t < d,
             exp(-lambda0 * t),
             exp(-lambda0 * d - lambda1_post * (t - d)))
    }
    f0 <- function(t) lambda0 * S0(t)
    f1 <- function(t) {
      ifelse(t < d,
             lambda0 * S1(t),
             lambda1_post * S1(t))
    }
    label <- sprintf("Delayed benefit (HR=%.2f, delay=%g months)", HR, d)

  } else if (type == "crossing") {
    if (is.null(HR_early) || is.null(HR_late) || is.null(changepoint))
      stop("Specify HR_early, HR_late, changepoint for type='crossing'.")
    cp <- changepoint
    lam1_e <- HR_early * lambda0
    lam1_l <- HR_late  * lambda0
    S0 <- function(t) exp(-lambda0 * t)
    S1 <- function(t) {
      ifelse(t < cp,
             exp(-lam1_e * t),
             exp(-lam1_e * cp - lam1_l * (t - cp)))
    }
    f0 <- function(t) lambda0 * S0(t)
    f1 <- function(t) {
      ifelse(t < cp,
             lam1_e * S1(t),
             lam1_l * S1(t))
    }
    label <- sprintf("Crossing hazards (HR=%.2f/%.2f, cp=%g months)",
                     HR_early, HR_late, cp)

  } else if (type == "frailty") {
    if (is.null(HR_cond))
      stop("Specify HR_cond for type='frailty'.")
    lambda1 <- HR_cond * lambda0
    # Marginal survival under Gamma frailty (no covariates for simplicity)
    # S_A(t) = (1 + theta*lambda_A*t)^(-1/theta)
    S0 <- function(t) (1 + theta * lambda0 * t)^(-1/theta)
    S1 <- function(t) (1 + theta * lambda1 * t)^(-1/theta)
    f0 <- function(t) lambda0 * (1 + theta * lambda0 * t)^(-1/theta - 1)
    f1 <- function(t) lambda1 * (1 + theta * lambda1 * t)^(-1/theta - 1)
    label <- sprintf("Gamma frailty   (cond.HR=%.2f, theta=%.2f)",
                     HR_cond, theta)
  }

  # ----------------------------------------------------------------
  # Numerical integration on grid
  # ----------------------------------------------------------------
  t_grid <- seq(1e-6, tau, length.out = n_grid)

  S0v <- S0(t_grid); S1v <- S1(t_grid)
  f0v <- f0(t_grid); f1v <- f1(t_grid)
  dv  <- f0v - f1v   # delta(t)

  C_tau  <- trapz(t_grid, t_grid * dv)
  Delta  <- trapz(t_grid, S1v - S0v)
  phi    <- trapz(t_grid, dv)
  t_bar  <- if (abs(phi) > 1e-10) C_tau / phi else NA
  SRR    <- if (!is.na(t_bar))    1 - t_bar / tau else NA

  # ----------------------------------------------------------------
  # t_star: sign-change time of delta(t) = f0(t) - f1(t)
  # Find where delta(t) changes sign within [0, tau]
  # ----------------------------------------------------------------
  sign_changes <- which(diff(sign(dv)) != 0)
  if (length(sign_changes) == 0) {
    t_star <- NA   # no sign change within [0, tau]
    t_star_note <- sprintf("> tau (delta(t) does not change sign in [0, %.0f])", tau)
  } else {
    # Linear interpolation for first sign change
    i  <- sign_changes[1]
    t_star <- t_grid[i] + (t_grid[i+1] - t_grid[i]) *
              abs(dv[i]) / (abs(dv[i]) + abs(dv[i+1]))
    t_star_note <- sprintf("%.2f months (within [0, tau])", t_star)
  }

  # ----------------------------------------------------------------
  # PH reference range (exponential, HR 0.50-0.90, tau=24, median=12)
  # Precomputed: t_bar/tau in [0.18, 0.28]
  # ----------------------------------------------------------------
  ph_low  <- 0.18   # HR=0.90
  ph_high <- 0.28   # HR=0.50
  tbar_ratio <- t_bar / tau

  # ----------------------------------------------------------------
  # Recommendation
  # ----------------------------------------------------------------
  if (is.na(SRR)) {
    recommendation <- "Cannot determine (phi ~ 0, no net treatment effect)"
    regime <- NA
  } else if (C_tau < 0) {
    recommendation <- paste(
      "C(tau) < 0: early harm precedes later benefit.",
      "RMST is structurally favored over flat benchmark.",
      "RMST recommended if estimand is appropriate."
    )
    regime <- "i"
  } else if (tbar_ratio <= ph_high + 0.05) {
    recommendation <- paste(
      sprintf("t_bar/tau = %.3f is within or below the PH reference range [%.3f, %.3f].",
              tbar_ratio, ph_low, ph_high),
      "RMST is competitive with log-rank.",
      "RMST recommended: interpretability advantage at no meaningful power cost."
    )
    regime <- "ii"
  } else if (tbar_ratio <= 0.50) {
    recommendation <- paste(
      sprintf("t_bar/tau = %.3f exceeds PH reference range [%.3f, %.3f].",
              tbar_ratio, ph_low, ph_high),
      "RMST faces moderate structural disadvantage.",
      "Consider log-rank or weighted log-rank; RMST still interpretable but power loss expected."
    )
    regime <- "iii-moderate"
  } else {
    recommendation <- paste(
      sprintf("t_bar/tau = %.3f >> PH reference range [%.3f, %.3f].",
              tbar_ratio, ph_low, ph_high),
      "RMST faces severe structural disadvantage (SRR < 0.5).",
      "Log-rank or Fleming-Harrington WLR strongly preferred.",
      "Covariate adjustment cannot repair this gap."
    )
    regime <- "iii-severe"
  }

  # ----------------------------------------------------------------
  # Proposition 3 upper bound (delayed type only)
  # ----------------------------------------------------------------
  srr_upper <- NA
  if (type == "delayed") {
    srr_upper <- 1 - delay / tau
  }

  # ----------------------------------------------------------------
  # Output
  # ----------------------------------------------------------------
  results <- list(
    type           = type,
    label          = label,
    tau            = tau,
    lambda0        = lambda0,
    median0        = log(2) / lambda0,
    C_tau          = C_tau,
    Delta_tau      = Delta,
    phi_tau        = phi,
    t_bar_delta    = t_bar,
    tbar_ratio     = tbar_ratio,
    SRR            = SRR,
    t_star         = t_star,
    t_star_note    = t_star_note,
    SRR_upper_bound= srr_upper,
    PH_range       = c(ph_low, ph_high),
    regime         = regime,
    recommendation = recommendation
  )

  if (verbose) {
    cat("\n", paste(rep("=", 60), collapse=""), "\n", sep="")
    cat("RMST Design-Stage Diagnostic\n")
    cat(paste(rep("=", 60), collapse=""), "\n", sep="")
    cat(sprintf("Scenario    : %s\n", label))
    cat(sprintf("tau         : %g months\n", tau))
    cat(sprintf("Control med : %.1f months\n", log(2)/lambda0))
    cat(paste(rep("-", 60), collapse=""), "\n", sep="")
    cat(sprintf("C(tau)      : %.4f\n", C_tau))
    cat(sprintf("Delta(tau)  : %.4f months\n", Delta))
    cat(sprintf("t_star      : %s  [sign-change of delta(t)]\n", t_star_note))
    cat(sprintf("t_bar_delta : %.2f months  [signal centroid]\n", t_bar))
    cat(sprintf("t_bar/tau   : %.3f  (PH range: [%.3f, %.3f])\n",
                tbar_ratio, ph_low, ph_high))
    cat(sprintf("SRR         : %.4f\n", SRR))
    if (!is.na(srr_upper))
      cat(sprintf("SRR upper   : %.4f  (Prop.3: 1 - d/tau = 1 - %g/%g)\n",
                  srr_upper, delay, tau))
    cat(paste(rep("-", 60), collapse=""), "\n", sep="")
    cat("RECOMMENDATION:\n")
    cat(strwrap(recommendation, width=58, prefix="  "), sep="\n")
    cat(paste(rep("=", 60), collapse=""), "\n\n", sep="")
  }

  invisible(results)
}

# Simple trapezoid integration
trapz <- function(x, y) sum(diff(x) * (y[-1] + y[-length(y)]) / 2)

# ----------------------------------------------------------------
# Demo
# ----------------------------------------------------------------
cat("========================================\n")
cat("DEMO: Four scenarios from Wu (2025)\n")
cat("========================================\n")

# DGM1: delayed benefit
rmst_diagnostic(type="delayed", median0=12, HR=0.55, tau=24, delay=4)

# DGM2: crossing hazards
rmst_diagnostic(type="crossing", median0=12,
                HR_early=1.3, HR_late=0.50,
                tau=24, changepoint=4)

# DGM4: frailty
rmst_diagnostic(type="frailty", median0=12,
                HR_cond=0.70, tau=24, theta=0.5)

# PH reference
rmst_diagnostic(type="exponential", median0=12, HR=0.70, tau=24)

# ----------------------------------------------------------------
# Note on frailty with covariates:
# The frailty type above uses marginal survival without covariates
# (r=1). If covariates are present, the true marginal survival
# requires integrating over the covariate distribution:
#   S_A(t) = E_X[(1 + theta*lambda_A*exp(eta)*t)^(-1/theta)]
# This requires numerical integration over X and is scenario-specific.
# The function above provides a conservative approximation.
# ----------------------------------------------------------------
