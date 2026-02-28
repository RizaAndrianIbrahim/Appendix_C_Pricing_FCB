# This appendix contains an R program that simulates Algorithm 1 to estimate solutions to the proposed flood catastrophe bond pricing model.
# This script implements:
# (1) ARIMAX/ARIMA-based forecasting for flood intensity, ONI, and SOI
# (2) Monte Carlo simulation of flood occurrence and loss processes
# (3) Valuation of flood catastrophe bonds under indemnity triggers
#
# NOTE:
# - The code below is a direct reformatting of the original program.
# - No algorithms, parameters, or logical steps have been altered.
# - Changes are limited to layout, spacing, and professional annotation
#   for journal supplementary material purposes.
# ============================================================== 

rm(list = ls())

# --------------------------------------------------------------
# Required libraries
# --------------------------------------------------------------
library(readxl)
library(tseries)
library(forecast)
library(stats)

# --------------------------------------------------------------
# Data input: flood intensity (lambda), ONI, and SOI
# --------------------------------------------------------------
data   <- read_excel(file.choose())
data   <- ts(data)

lambda <- data[, 1]
oni    <- data[, 2]
soi    <- data[, 3]

# --------------------------------------------------------------
# Box–Cox transformation for flood intensity
# --------------------------------------------------------------
l11_lambda  <- BoxCox.lambda(lambda)
bc1_lambda  <- (lambda^l11_lambda - 1) / l11_lambda

# --------------------------------------------------------------
# Model orders
# --------------------------------------------------------------
# ARIMAX model for flood intensity (lambda)
p_lambda <- 5
d_lambda <- 1
q_lambda <- 5

# ARIMA model for ONI
p_oni <- 4
d_oni <- 1
q_oni <- 5

# ARIMA model for SOI
p_soi <- 4
d_soi <- 1
q_soi <- 5

# --------------------------------------------------------------
# Model estimation
# --------------------------------------------------------------
m_lambda <- Arima(
  bc1_lambda,
  c(p_lambda, d_lambda, q_lambda),
  xreg  = cbind(oni, soi),
  method = "CSS-ML"
)

m_oni <- Arima(
  oni,
  c(p_oni, d_oni, q_oni),
  method = "CSS-ML"
)

m_soi <- Arima(
  soi,
  c(p_soi, d_soi, q_soi),
  method = "CSS-ML"
)

# --------------------------------------------------------------
# Standard deviations of ARIMA error 
# --------------------------------------------------------------
s_lambda <- sd(residuals(m_lambda))
s_oni    <- sd(residuals(m_oni))
s_soi    <- sd(residuals(m_soi))

# --------------------------------------------------------------
# Gamma parameters for individual flood losses
# --------------------------------------------------------------
alpha_G <- 0.62549
beta_G  <- 4.1955

# --------------------------------------------------------------
# Flood bond parameters
# --------------------------------------------------------------
T       <- 5        # Bond maturity
C       <- 0.1      # Annual coupon rate
R       <- 1        # Redemption value
kappa_1 <- 0        # Coupon protection fraction
theta_1 <- 0.5      # Redemption protection fraction
mu      <- 569.2214 # Cumulative loss trigger

# --------------------------------------------------------------
# Cox-Ingersoll-Ross model parameters
# --------------------------------------------------------------
kappa_2 <- 0.20845
theta_2 <- 0.08285
sigma   <- 0.10944

gamma_2 <- sqrt(kappa_2^2 + 2 * sigma^2)

B_0 <- numeric(T)

for (k in 1:T) {

  D_0k <- 2*(exp(gamma_2*k) - 1) /
          (2*gamma_2 +
           (kappa_2 + gamma_2)*(exp(gamma_2*k) - 1))

  C_0k <- (2*gamma_2 *
           exp((kappa_2 + gamma_2)*k/2) /
           (2*gamma_2 +
            (kappa_2 + gamma_2)*(exp(gamma_2*k) - 1))) ^
           (2*kappa_2*theta_2/sigma^2)

  B_0[k] <- C_0k * exp(-D_0k * 0.0610416666666667)
}

# --------------------------------------------------------------
# Monte Carlo configuration
# --------------------------------------------------------------
W <- 100          # Number of replications
M <- 10000        # Samples per replication

# --------------------------------------------------------------
# Utility functions
# --------------------------------------------------------------
valid_vec <- function(x) all(is.finite(x))
valid_num <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)

# --------------------------------------------------------------
# Storage matrices
# --------------------------------------------------------------
PL <- matrix(0, nrow = W, ncol = T)
VT <- matrix(0, nrow = W, ncol = T)

# --------------------------------------------------------------
# Main Monte Carlo loop
# --------------------------------------------------------------
for (w in 1:W) {

  frcst_loss <- matrix(0, nrow = M, ncol = T)

  for (m in 1:M) {

    # ----------------------------------------------------------
    # Forecast ONI with noise
    # ----------------------------------------------------------
    repeat {
      fc <- tryCatch(
        as.numeric(forecast(m_oni, h = T)$mean + rnorm(T, 0, s_oni)),
        error = function(e) rep(NA_real_, T)
      )
      if (!valid_vec(fc)) next
      raw_frcst_oni <- fc
      break
    }

    # ----------------------------------------------------------
    # Forecast SOI with noise
    # ----------------------------------------------------------
    repeat {
      fc <- tryCatch(
        as.numeric(forecast(m_soi, h = T)$mean + rnorm(T, 0, s_soi)),
        error = function(e) rep(NA_real_, T)
      )
      if (!valid_vec(fc)) next
      raw_frcst_soi <- fc
      break
    }

    # ----------------------------------------------------------
    # Truncation of climate indices
    # ----------------------------------------------------------
    frcst_oni <- pmax(pmin(raw_frcst_oni, 12), 0)
    frcst_soi <- pmax(pmin(raw_frcst_soi, 12), 0)

    xreg_future <- cbind(oni = frcst_oni, soi = frcst_soi)

    # ----------------------------------------------------------
    # Forecast flood intensity and inverse Box–Cox transform
    # ----------------------------------------------------------
    repeat {
      frcst_m_lambda <- tryCatch(
        as.numeric(
          forecast(m_lambda, xreg = xreg_future, h = T)$mean +
            rnorm(T, 0, s_lambda)
        ),
        error = function(e) rep(NA, T)
      )

      if (!is.numeric(frcst_m_lambda) || any(!is.finite(frcst_m_lambda))) next

      frcst_lambda <- tryCatch(
        (l11_lambda * frcst_m_lambda + 1)^(1 / l11_lambda),
        error = function(e) rep(NA, T)
      )

      if (is.numeric(frcst_lambda) &&
          all(is.finite(frcst_lambda)) &&
          all(frcst_lambda >= 0)) break
    }

    # ----------------------------------------------------------
    # Flood occurrence and loss simulation
    # ----------------------------------------------------------
    frcst_flood <- numeric(T)
    frcst_f <- 0
    frcst_l <- 0

    for (k in 1:T) {
      repeat {
        f_new <- frcst_f + rpois(1, frcst_lambda[k])
        if (!valid_num(f_new)) next

        l_new <- frcst_l + sum(rgamma(f_new, shape = alpha_G, scale = beta_G))
        if (!valid_num(l_new)) next

        frcst_flood[k]  <- f_new
        frcst_loss[m,k] <- l_new

        frcst_f <- f_new
        frcst_l <- l_new
        break
      }
    }
  }

  # --------------------------------------------------------------
  # Protection level and bond valuation
  # --------------------------------------------------------------
  PL[w, ] <- colMeans(frcst_loss <= mu)

  Ck <- 0
  for (k in 1:T) {
    Ck <- Ck + C * ((1 - kappa_1) * PL[w, k] + kappa_1) * B_0[k]
    RT <- R * ((1 - theta_1) * PL[w, k] + theta_1) * B_0[k]
    VT[w, k] <- Ck + RT
  }
}

# --------------------------------------------------------------
# Output (copied to clipboard)
# --------------------------------------------------------------
library(clipr)
write_clip(VT)
