## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  dpi = 96,
  out.width = "100%"
)

## ----setup--------------------------------------------------------------------
library(Rtwalk)
set.seed(2024)

## ----helper-------------------------------------------------------------------
compare <- function(fit, truth, burnin_frac = 0.2) {
  s <- summary(fit, burnin_frac = burnin_frac)
  data.frame(
    parameter = s$Parameter, true_mean = truth$mean, estimate = s$Mean,
    mcse = s$MCSE, true_sd = truth$sd, sd = s$SD, ess_bulk = s$ESS_bulk,
    rhat = s$Rhat
  )
}

## ----normal-------------------------------------------------------------------
fit <- twalk(function(x) dnorm(x, log = TRUE), n_iter = 20000,
             x0 = -2, xp0 = 2, show_progress = FALSE)
fit
compare(fit, list(mean = 0, sd = 1))

## ----bivariate----------------------------------------------------------------
mu <- c(1, -0.5)
Sigma <- matrix(c(4, 2.1, 2.1, 2.25), 2)
log_bvn <- function(x) mvtnorm::dmvnorm(x, mu, Sigma, log = TRUE)

starts <- function() matrix(rnorm(8, sd = 5), nrow = 4)
fit_bvn <- twalk(log_bvn, n_iter = 10000, x0 = starts(), xp0 = starts(),
                 n_chains = 4, n_cores = 1, show_progress = FALSE)
compare(fit_bvn, list(mean = mu, sd = sqrt(diag(Sigma))))

## ----bivariate-plot, fig.height = 6-------------------------------------------
plot(fit_bvn)

## ----gamma--------------------------------------------------------------------
y <- c(2.1, 0.7, 1.9, 3.4, 0.4, 1.2, 2.8, 1.5)
# Exponential likelihood with rate lambda and a Gamma(2, 1) prior:
# the posterior is Gamma(2 + n, 1 + sum(y)).
log_post_rate <- function(theta) {
  if (theta[["lambda"]] <= 0) return(-Inf)
  sum(dexp(y, theta[["lambda"]], log = TRUE)) + dgamma(theta[["lambda"]], 2, 1, log = TRUE)
}

fit_rate <- twalk(log_post_rate, n_iter = 20000, x0 = c(lambda = 0.5),
                  xp0 = c(lambda = 1), show_progress = FALSE)
shape <- 2 + length(y); rate <- 1 + sum(y)
compare(fit_rate, list(mean = shape / rate, sd = sqrt(shape) / rate))

## ----ten----------------------------------------------------------------------
d <- 10
Sigma10 <- 0.7^abs(outer(1:d, 1:d, "-"))
mu10 <- 1:d
log_mvn <- function(x) mvtnorm::dmvnorm(x, mu10, Sigma10, log = TRUE)

fit10 <- twalk(log_mvn, n_iter = 1e5, x0 = rep(0, d), xp0 = rep(2, d),
               show_progress = FALSE, thin = 10)
compare(fit10, list(mean = mu10, sd = rep(1, d)))

## ----logistic-----------------------------------------------------------------
n_obs <- 1000
X <- cbind(1, rnorm(n_obs), rnorm(n_obs, sd = 1.5))
beta_true <- c(0.5, -1.2, 0.8)
y_bin <- rbinom(n_obs, 1, plogis(X %*% beta_true))

log_post_logit <- function(beta, X, y) {
  eta <- drop(X %*% beta)
  sum(y * eta - pmax(eta, 0) - log1p(exp(-abs(eta)))) +
    sum(dnorm(beta, 0, 5, log = TRUE))
}

fit_logit <- twalk(log_post_logit, n_iter = 30000,
                   x0 = c(b0 = 0, b1 = 0, b2 = 0), xp0 = c(0.2, -0.2, 0.1),
                   X = X, y = y_bin, show_progress = FALSE)

glm_fit <- glm(y_bin ~ X - 1, family = binomial)
compare(fit_logit, list(mean = unname(coef(glm_fit)),
                        sd = unname(sqrt(diag(vcov(glm_fit))))))

## ----continue-----------------------------------------------------------------
fit_more <- twalk_continue(fit_logit, log_post_logit, n_iter = 30000,
                           X = X, y = y_bin, show_progress = FALSE)
fit_more$n_iter
summary(fit_more)$ESS_bulk

## ----interop------------------------------------------------------------------
coda::effectiveSize(coda::as.mcmc(fit_logit, burnin_frac = 0.2))

if (requireNamespace("posterior", quietly = TRUE)) {
  draws <- posterior::as_draws_df(fit_bvn, burnin_frac = 0.2)
  posterior::summarise_draws(draws)
}

## ----rates--------------------------------------------------------------------
fit10$kernel_statistics

