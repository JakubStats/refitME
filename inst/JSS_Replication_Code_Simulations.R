#----------------------------------------------------------------------------
# JSS_Replication_Code_Examples.R
#
# This file reproduces all the results for the simulations (results for
# Figures 1-2, Figures 5-9 (Section 5), and Figures 12-13 Appendix B).
#----------------------------------------------------------------------------

# Figure 1: Checking the effective sample size (ESS) across model dimensions and measurement error variance.

library(refitME)

setwd(system.file(package = 'refitME'))
source('MCEM_prog.R')

set.seed(2019)

B <- 50

n <- 10000

dim_vec <- c(1:10)

sigma.sq.e <- diag(c(rep(1, 10))) # True predictor covariance (sigma).

alpha <- c(rep(0, 10))

beta_vec <- c(0.5, rep(1, 10))

sigma.sq.u_big <- diag(c(rep(0.5, 10)))

x1 <- rnorm(n, alpha[1], sd = sqrt(sigma.sq.e[1, 1]))
x2 <- rnorm(n, alpha[2], sd = sqrt(sigma.sq.e[2, 2]))
x3 <- rnorm(n, alpha[3], sd = sqrt(sigma.sq.e[3, 3]))
x4 <- rnorm(n, alpha[4], sd = sqrt(sigma.sq.e[4, 4]))
x5 <- rnorm(n, alpha[5], sd = sqrt(sigma.sq.e[5, 5]))
x6 <- rnorm(n, alpha[6], sd = sqrt(sigma.sq.e[6, 6]))
x7 <- rnorm(n, alpha[7], sd = sqrt(sigma.sq.e[7, 7]))
x8 <- rnorm(n, alpha[8], sd = sqrt(sigma.sq.e[8, 8]))
x9 <- rnorm(n, alpha[9], sd = sqrt(sigma.sq.e[9, 9]))
x10 <- rnorm(n, alpha[10], sd = sqrt(sigma.sq.e[10, 10]))

X_big <- cbind(rep(1, n), x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)

w1 <- x1 + rnorm(n, 0, sd = sqrt(sigma.sq.u_big[1, 1]))
w2 <- x2 + rnorm(n, 0, sd = sqrt(sigma.sq.u_big[2, 2]))
w3 <- x3 + rnorm(n, 0, sd = sqrt(sigma.sq.u_big[3, 3]))
w4 <- x4 + rnorm(n, 0, sd = sqrt(sigma.sq.u_big[4, 4]))
w5 <- x5 + rnorm(n, 0, sd = sqrt(sigma.sq.u_big[5, 5]))
w6 <- x6 + rnorm(n, 0, sd = sqrt(sigma.sq.u_big[6, 6]))
w7 <- x7 + rnorm(n, 0, sd = sqrt(sigma.sq.u_big[7, 7]))
w8 <- x8 + rnorm(n, 0, sd = sqrt(sigma.sq.u_big[8, 8]))
w9 <- x9 + rnorm(n, 0, sd = sqrt(sigma.sq.u_big[9, 9]))
w10 <- x10 + rnorm(n, 0, sd = sqrt(sigma.sq.u_big[10, 10]))

W_big <- cbind(w1, w2 ,w3, w4, w5, w6, w7, w8, w9, w10)

eff_vec1 <- c()
eff_vec2 <- c()

for(i in 1:length(dim_vec)) {
  sigma.sq.u <- sigma.sq.u_big[1:i, 1:i]

  W <- as.matrix(W_big[, 1:i])
  colnames(W) <- colnames(W_big)[1:i]

  X <- X_big[, 1:(1 + i)]
  colnames(X) <- colnames(X_big)[1:(1 + i)]

  eta <- X%*%beta_vec[1:(1 + i)]

  mu_Y <- eta
  Y1 <- c(mu_Y) + rnorm(n, 0, 1)

  mu_Y <- exp(eta)/(1 + exp(eta))
  Y2 <- rbinom(n, 1, prob = mu_Y)

  dat1 <- data.frame(cbind(Y1, X[, -1], W))
  dat2 <- data.frame(cbind(Y2, X[, -1], W))

  colnames(dat1) <- c("Y1", colnames(X)[-1], colnames(W))
  colnames(dat2) <- c("Y2", colnames(X)[-1], colnames(W))

  if (dim_vec[i] == 1) {
    mod_naiv1 <- lm(Y1 ~ w1, x = TRUE, data = dat1)
    mod_naiv2 <- glm(Y2 ~ w1, x = TRUE, family = binomial, data = dat2)
  }

  if (dim_vec[i] == 2) {
    mod_naiv1 <- lm(Y1 ~ w1 + w2, x = TRUE, data = dat1)
    mod_naiv2 <- glm(Y2 ~ w1 + w2, x = TRUE, family = binomial, data = dat2)
  }

  if (dim_vec[i] == 3) {
    mod_naiv1 <- lm(Y1 ~ w1 + w2 + w3, x = TRUE, data = dat1)
    mod_naiv2 <- glm(Y2 ~ w1 + w2 + w3, x = TRUE, family = binomial, data = dat2)
  }

  if (dim_vec[i] == 4) {
    mod_naiv1 <- lm(Y1 ~ w1 + w2 + w3 + w4, x = TRUE, data = dat1)
    mod_naiv2 <- glm(Y2 ~ w1 + w2 + w3 + w4, x = TRUE, family = binomial, data = dat2)
  }

  if (dim_vec[i] == 5) {
    mod_naiv1 <- lm(Y1 ~ w1 + w2 + w3 + w4 + w5, x = TRUE, data = dat1)
    mod_naiv2 <- glm(Y2 ~ w1 + w2 + w3 + w4 + w5, x = TRUE, family = binomial, data = dat2)
  }

  if (dim_vec[i] == 6) {
    mod_naiv1 <- lm(Y1 ~ w1 + w2 + w3 + w4 + w5 + w6, x = TRUE, data = dat1)
    mod_naiv2 <- glm(Y2 ~ w1 + w2 + w3 + w4 + w5 + w6, x = TRUE, family = binomial, data = dat2)
  }

  if (dim_vec[i] == 7) {
    mod_naiv1 <- lm(Y1 ~ w1 + w2 + w3 + w4 + w5 + w6 + w7, x = TRUE, data = dat1)
    mod_naiv2 <- glm(Y2 ~ w1 + w2 + w3 + w4 + w5 + w6 + w7, x = TRUE, family = binomial, data = dat2)
  }

  if (dim_vec[i] == 8) {
    mod_naiv1 <- lm(Y1 ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8, x = TRUE, data = dat1)
    mod_naiv2 <- glm(Y2 ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8, x = TRUE, family = binomial, data = dat2)
  }
  if (dim_vec[i] == 9) {
    mod_naiv1 <- lm(Y1 ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9, x = TRUE, data = dat1)
    mod_naiv2 <- glm(Y2 ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9, x = TRUE, family = binomial, data = dat2)
  }
  if (dim_vec[i] == 10) {
    mod_naiv1 <- lm(Y1 ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, x = TRUE, data = dat1)
    mod_naiv2 <- glm(Y2 ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10, x = TRUE, family = binomial, data = dat2)
  }

  est1 <- refitME(mod_naiv1, sigma.sq.u, B)
  est2 <- refitME(mod_naiv2, sigma.sq.u, B)

  eff_vec1 <- c(eff_vec1, mean(est1$eff.samp.size, na.rm = T)/B)
  eff_vec2 <- c(eff_vec2, mean(est2$eff.samp.size, na.rm = T)/B)
}

dats1 <- as.data.frame(cbind((dim_vec), (eff_vec1), (1/dim_vec)))
dats2 <- as.data.frame(cbind((dim_vec), (eff_vec2), (1/dim_vec)))

colnames(dats1) = colnames(dats2) = c("dim_vec", "ESS", "dim_inv")

p1 <- ggplot(dats1, aes(x = dim_vec, y = ESS)) + geom_point() + geom_line(size = 1) +
  geom_line(aes(x = dim_vec, y = dim_inv), color = "red", linetype = "dotted", colour ="ESS/B", size = 2) +
  labs(x = bquote("p"), y = bquote("ESS/B"), title = "Gaussian response") +
  theme(axis.title.y = element_text(size = 20, angle = 90), text = element_text(size = 20),
        axis.text = element_text(size = 20), axis.title = element_text(size = 20, face = "bold")) +
  scale_y_continuous(trans = "log", breaks = trans_breaks("log", function(x) round(exp(x), digits = 2))) +
  scale_x_continuous(trans = "log", breaks = trans_breaks("log", function(x) round(exp(x), digits = 0)))

p2 <- ggplot(dats2, aes(x = dim_vec, y = ESS)) + geom_point() + geom_line(size = 1) +
  geom_line(aes(x = dim_vec, y = dim_inv), color = "red", linetype = "dotted", size = 2) +
  labs(x = bquote("p"), y = bquote("ESS/B"), title = "Binary response") +
  theme(axis.title.y = element_text(size = 20, angle = 90), legend.title = element_blank(), text = element_text(size = 20),
        legend.position = "none", axis.text = element_text(size = 20), axis.title = element_text(size = 20, face = "bold")) +
  scale_y_continuous(trans = "log", breaks = trans_breaks("log", function(x) round(exp(x), digits = 2))) +
  scale_x_continuous(trans = "log", breaks = trans_breaks("log", function(x) round(exp(x), digits = 0)))

multiplot(p1, p2, layout = matrix(c(1, 2), ncol = 2, nrow = 1, byrow = T))

#--------------------------------------------------------------

# Figure 2: Effective sample size (EFF) against measurement error variance.

library(refitME)

setwd(system.file(package = 'refitME'))
source('MCEM_prog.R')

set.seed(2019)

B <- 50

n <- 10000

sigma.sq.u_vec <- c(seq(0.1, 1, 0.1))

sigma.sq.e <- diag(c(rep(1, 2))) # True predictor covariance (sigma).

alpha <- c(rep(0, 2))

beta_vec <- c(0.5, rep(1, 2))

eff_vec1 <- eff_vec2 <- c()

for(i in 1:length(sigma.sq.u_vec)) {
  sigma.sq.u <- diag(c(rep(sigma.sq.u_vec[i], 2)))

  x1 <- rnorm(n, alpha[1], sd = 1)
  x2 <- rnorm(n, alpha[2], sd = 1)

  X <- cbind(rep(1, n), x1, x2)

  w1 <- x1 + rnorm(n, 0, sd = sqrt(sigma.sq.u[1, 1]))
  w2 <- x2 + rnorm(n, 0, sd = sqrt(sigma.sq.u[1, 1]))

  W <- cbind(w1, w2)

  # Gaussian response.

  eta <- X%*%beta_vec
  mu_Y <- eta
  Y <- c(mu_Y) + rnorm(n, 0, 1)
  dat <- data.frame(cbind(Y, X, W))
  names(dat)[2] <- "(Intercept)"

  mod_naiv1 <- lm(Y ~ w1 + w2, x = TRUE, data = dat)

  est <- refitME(mod_naiv1, sigma.sq.u, B)
  eff_vec1 <- c(eff_vec1, mean(est$eff.samp.size, na.rm = T)/B)

  # Binomial response.

  mu_Y <- exp(eta)/(1 + exp(eta))
  Y <- rbinom(n, 1, prob = mu_Y)
  dat <- data.frame(cbind(Y, X, W))
  names(dat)[2] <- "(Intercept)"

  mod_naiv1 <- glm(Y ~ w1 + w2, x = TRUE, family = binomial, data = dat)

  est <- refitME(mod_naiv1, sigma.sq.u, B)
  eff_vec2 <- c(eff_vec2, mean(est$eff.samp.size, na.rm = T)/B)
}

dat1 <- as.data.frame(cbind(sigma.sq.u_vec, eff_vec1))
dat2 <- as.data.frame(cbind(sigma.sq.u_vec, eff_vec2))

colnames(dat1) = colnames(dat2) = c("sigma.sq.u_big", "ESS")

p1 <- ggplot(dat1, aes(x = sigma.sq.u_vec, y = ESS)) + geom_point() + geom_line(size = 1) +
  labs(x = bquote("ME variance"), y = bquote("ESS/B"), title = "Gaussian response") +
  xlab(expression(sigma[u]^2)) +
  theme(axis.title.y = element_text(size = 20, angle = 90), text = element_text(size = 20),
        axis.text = element_text(size = 20), axis.title = element_text(size = 20, face = "bold"))

p2 <- ggplot(dat2, aes(x = sigma.sq.u_vec, y = ESS)) + geom_point() + geom_line(size = 1) +
  labs(x = bquote("ME variance"), y = bquote("ESS/B"), title = "Binary response") +
  xlab(expression(sigma[u]^2)) +
  theme(axis.title.y = element_text(size = 20, angle = 90), legend.title = element_blank(), text = element_text(size = 20),
        legend.position = "none", axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, face = "bold"))

par(mfrow = c(1, 2), las = 1)

multiplot(p1, p2, layout = matrix(c(1, 2), ncol = 2, nrow = 1, byrow = T))

#-----------------------------------------------------------------------------

# Sections 5.1 and 5.2: Simulation study examining bias, efficiency and
# confidence interval coverage for model parameters, and prediction.

library(refitME)
library(mgcv)

setwd(system.file(package = 'refitME'))
source('MCEM_prog.R')

set.seed(2016)

N.sim <- 200

B <- 50

n <- 1000
n_train <- 800
n_test <- n - n_train

# Set true parameter values here.

sigma.sq.e <- c(1, 1) # True predictor covariance (sigma).
alpha <- c(0, 0)

beta <- c(0.5, 1, -0.3); par.int <- 3

par.int <- "spline" # Remove/toggle the hash (#) if you want to fit spline models.

scen_par <- 0.5  # Additional shift in "future" predictor values.

test.W <- F

sigma.sq.u_vec <- c(0.0001, 0.05, 0.1, 0.25, 0.5)

bias_mat <- c()
CP_mat <- c()
SE_mat <- c()
SD_mat <- c()
RMSE_mat <- c()
CV_mat1 <- c()
CV_mat2 <- c()

#family <- "binomial"
family <- "poisson"

start <- Sys.time()

for(i in 1:length(sigma.sq.u_vec)) {
  nocrash <- c()

  sigma.sq.u1 <- sigma.sq.u_vec[i]  # Measurement error variance.

  beta_mat0 <- c()
  beta_mat1 <- c()
  beta_mat2 <- c()
  beta_mat3 <- c()

  beta_mat.se0 <- c()
  beta_mat.se1 <- c()
  beta_mat.se2 <- c()
  beta_mat.se3 <- c()

  RMSE_CV0a <- c()
  RMSE_CV1a <- c()
  RMSE_CV2a <- c()
  RMSE_CV3a <- c()

  RMSE_CV0b <- c()
  RMSE_CV1b <- c()
  RMSE_CV2b <- c()
  RMSE_CV3b <- c()

  c.sim <- TRUE
  k <- 0

  print("___________________________________")
  print(sigma.sq.u1)

  while(c.sim) {
    x1 <- rnorm(n, alpha[1], sd = sqrt(sigma.sq.e[1]))

    if (par.int == 3) X <- cbind(rep(1, n), x1, x1^2)
    if (par.int == "spline") X <- NA

    w1 <- x1 + rnorm(n, 0, sd = sqrt(sigma.sq.u1))  # Add/contaminate the predictor with error.

    if (par.int == 3) {
      W <- cbind(rep(1, n), w1, w1^2)
      colnames(W) <- c("Intercept", "w1", "I(w1^2)")
    }

    if (par.int == "spline") W <- w1

    if (par.int != "spline") eta <- X%*%beta
    if (par.int == "spline" & family == "binomial") {
      eta <- cos(2*x1 + 0.25)
      mu_Y <- logit_fun(eta)
      Y <- rbinom(n, 1, prob = mu_Y)
    }
    if (par.int == "spline" & family == "poisson") {
      eta <- cos(2*x1 + 0.25)
      mu_Y <- exp(eta)
      Y <- rpois(n, lambda = mu_Y)
    }

    dat <- data.frame(cbind(Y, w1, x1))

    sigma.sq.u <- sigma.sq.u1
    sigma.sq.u_test <- sigma.sq.u1

    dat_train <- dat[1:n_train, ]
    dat_test <- data.frame(dat[(n_train + 1):n, ])
    eta_test <- eta[(n_train + 1):n]

    x1_test <- dat_test$x1

    if (par.int == 3) {
      X_test <- cbind(rep(1, nrow(dat_test)), x1_test, x1_test^2)
      W_test <- cbind(rep(1, nrow(dat_test)), dat_test$w1, dat_test$w1^2)
    }

    # Create new data (containing X, W and Y) as a future climate scenario i.e., add an extra bit to X_test.

    x1_test_scen <- x1_test + scen_par

    if (par.int == 3) X_test_scen <- cbind(rep(1, nrow(dat_test)), x1_test_scen, x1_test_scen^2)
    if (par.int == "spline") X_test_scen <- cbind(x1_test_scen)

    if (par.int != "spline") eta_test_scen <- as.matrix(X_test_scen)%*%beta
    if (par.int == "spline") eta_test_scen <- cos(x1_test_scen*2 + 0.25)

    w1_test_scen <- x1_test_scen + rnorm(n_test, 0, sd = sqrt(sigma.sq.u_test))

    if (par.int == 3) W_test_scen <- cbind(rep(1, n_test), w1_test_scen, w1_test_scen^2)

    dat_test_scen <- data.frame(cbind(w1_test_scen, x1_test_scen))
    colnames(dat_test_scen) <- colnames(dat_train)[-1]

    # Start fitting models here.

    if (par.int == 3) {
      mod_true1 <- glm(Y ~ x1 + I(x1^2), x = TRUE, family = binomial, data = dat_train)

      Y <- mod_true1$mod$Y
      x <- model.matrix(mod_true1)

      mod_naiv1 <- glm(Y ~ w1 + I(w1^2), x = TRUE, family = binomial, data = dat_train)

      mod_simex1 <- simex(mod_naiv1, SIMEXvariable = c("w1", "I(w1^2)"),
                          measurement.error = cbind(sqrt(sigma.sq.u), sqrt(sigma.sq.u))) # SIMEX.
    }

    if (par.int == "spline") {
      if (family == "binomial") {
        mod_true1 <- gam(Y ~ s(x1), family = binomial, data = dat_train, gamma = 1.4, select = T)
        mod_naiv1 <- gam(Y ~ s(w1), family = binomial, data = dat_train, gamma = 1.4, select = T)
      }
      if (family == "poisson") {
        mod_true1 <- gam(Y ~ s(x1), family = poisson, data = dat_train, gamma = 1.4, select = T)
        mod_naiv1 <- gam(Y ~ s(w1), family = poisson, data = dat_train, gamma = 1.4, select = T)
      }
    }

    est <- try(refitME(mod_naiv1, sigma.sq.u, B), silent = TRUE)  # MCEM.

    if (class(est)[1] == "try-error") next
    nocrash <- c(nocrash, 1)

    k <- k + 1
    print(k)

    if (par.int != "spline") {
      beta_mat0 <- rbind(beta_mat0, mod_true1$coef)
      beta_mat1 <- rbind(beta_mat1, mod_naiv1$coef)
      beta_mat2 <- rbind(beta_mat2, mod_simex1$coef)
      beta_mat3 <- rbind(beta_mat3, est$coef)

      beta_mat.se0 <- rbind(beta_mat.se0, sqrt(diag(vcov(mod_true1))))
      beta_mat.se1 <- rbind(beta_mat.se1, sqrt(diag(vcov(mod_naiv1))))
      beta_mat.se2 <- rbind(beta_mat.se2, sqrt(diag(mod_simex1$variance.jackknife)))
      beta_mat.se3 <- rbind(beta_mat.se3, est$se)

      # CV/MSE on test data.

      eta_pred_true <- as.matrix(X_test)%*%mod_true1$coef
      RMSE_CV0a <- c(RMSE_CV0a, sqrt(mean((eta_test - eta_pred_true)^2)))
      eta_pred_true <- as.matrix(X_test_scen)%*%mod_true1$coef
      RMSE_CV0b <- c(RMSE_CV0b, sqrt(mean((eta_test_scen - eta_pred_true)^2)))

      if (test.W == F) eta_pred_naiv <- as.matrix(X_test)%*%mod_naiv1$coef
      if (test.W == T) eta_pred_naiv <- as.matrix(W_test)%*%mod_naiv1$coef
      RMSE_CV1a <- c(RMSE_CV1a, sqrt(mean((eta_test - eta_pred_naiv)^2)))
      if (test.W == F) eta_pred_naiv <- as.matrix(X_test_scen)%*%mod_naiv1$coef
      if (test.W == T) eta_pred_naiv <- as.matrix(W_test_scen)%*%mod_naiv1$coef
      RMSE_CV1b <- c(RMSE_CV1b, sqrt(mean((eta_test_scen - eta_pred_naiv)^2)))

      if (test.W == F) eta_pred_simex <- as.matrix(X_test)%*%mod_simex1$coef
      if (test.W == T) eta_pred_simex <- as.matrix(W_test)%*%mod_simex1$coef
      RMSE_CV2a <- c(RMSE_CV2a, sqrt(mean((eta_test - eta_pred_simex)^2)))
      if (test.W == F) eta_pred_simex <- as.matrix(X_test_scen)%*%mod_simex1$coef
      if (test.W == T) eta_pred_simex <- as.matrix(W_test_scen)%*%mod_simex1$coef
      RMSE_CV2b <- c(RMSE_CV2b, sqrt(mean((eta_test_scen - eta_pred_simex)^2)))

      if (test.W == F) eta_pred_mcem <- as.matrix(X_test)%*%est$coef
      if (test.W == T) eta_pred_mcem <- as.matrix(W_test)%*%est$coef
      RMSE_CV3a <- c(RMSE_CV3a, sqrt(mean((eta_test - eta_pred_mcem)^2)))
      if (test.W == F) eta_pred_mcem <- as.matrix(X_test_scen)%*%est$coef
      if (test.W == T) eta_pred_mcem <- as.matrix(W_test_scen)%*%est$coef
      RMSE_CV3b <- c(RMSE_CV3b, sqrt(mean((eta_test_scen - eta_pred_mcem)^2)))
    }

    if (par.int == "spline") {
      if (test.W == F) dat_test$w1 <- dat_test$x1

      eta_pred_true <- predict(mod_true1, newdata = dat_test, type = "link")
      RMSE_CV0a <- c(RMSE_CV0a, sqrt(mean((eta_test - eta_pred_true)^2)))

      eta_pred_naiv <- predict(mod_naiv1, newdata = dat_test, type = "link")
      RMSE_CV1a <- c(RMSE_CV1a, sqrt(mean((eta_test - eta_pred_naiv)^2)))

      eta_pred_mcem <- predict(est, newdata = dat_test, type = "link")
      RMSE_CV3a <- c(RMSE_CV3a, sqrt(mean((eta_test - eta_pred_mcem)^2)))
    }

    if (length(nocrash) >= N.sim) c.sim <- FALSE
  }

  if (par.int != "spline") {
    SE_mat <- cbind(beta_mat.se0[, par.int], beta_mat.se1[, par.int], beta_mat.se2[, par.int], beta_mat.se3[, par.int])
    SD_mat <- c(sd(beta_mat0[, par.int], na.rm = T), sd(beta_mat1[, par.int], na.rm = T), sd(beta_mat2[, par.int], na.rm = T), sd(beta_mat3[, par.int], na.rm = T), sd(beta_mat3[, par.int], na.rm = T))

    bias_mat1 <- c(mean((beta_mat0[, par.int] - beta[par.int])/beta[par.int], na.rm = T), mean((beta_mat1[, par.int] - beta[par.int])/beta[par.int], na.rm = T), mean((beta_mat2[, par.int] - beta[par.int])/beta[par.int], na.rm = T), mean((beta_mat3[, par.int] - beta[par.int])/beta[par.int], na.rm = T))
    RMSE_mat1 <- c(sqrt(mean((beta_mat0[, par.int] - beta[par.int])^2, na.rm = T)), sqrt(mean((beta_mat1[, par.int] - beta[par.int])^2, na.rm = T)), sqrt(mean((beta_mat2[, par.int] - beta[par.int])^2, na.rm = T)), sqrt(mean((beta_mat3[, par.int] - beta[par.int])^2, na.rm = T)))
    CP_mat_SE3 <- which(is.nan(beta_mat.se3[, par.int]) | is.na(beta_mat.se3[, par.int]))
    if (length(CP_mat_SE3) == 0) CP_mat3 <- cov.pc(beta, beta_mat3, beta_mat.se3, N.sim - length(CP_mat_SE3))[par.int]
    if (length(CP_mat_SE3) > 0) CP_mat3 <- cov.pc(beta, beta_mat3[-CP_mat_SE3, ], beta_mat.se3[-CP_mat_SE3, ], N.sim - length(CP_mat_SE3))[par.int]
    CP_mat1 <- c(cov.pc(beta, beta_mat0, beta_mat.se0, N.sim)[par.int], cov.pc(beta, beta_mat1, beta_mat.se1, N.sim)[par.int], cov.pc(beta, beta_mat2, beta_mat.se2, N.sim)[par.int], CP_mat3)

    bias_mat <- rbind(bias_mat, bias_mat1)
    RMSE_mat <- rbind(RMSE_mat, RMSE_mat1)
    CP_mat <- rbind(CP_mat, CP_mat1)
  }

  if (par.int != "spline") {
    CV_mat1 <- rbind(CV_mat1, c(mean(RMSE_CV0a, na.rm = T), mean(RMSE_CV1a, na.rm = T), mean(RMSE_CV2a, na.rm = T), mean(RMSE_CV3a, na.rm = T)))
    CV_mat2 <- rbind(CV_mat2, c(mean(RMSE_CV0b, na.rm = T), mean(RMSE_CV1b, na.rm = T), mean(RMSE_CV2b, na.rm = T), mean(RMSE_CV3b, na.rm = T)))
  }

  if (par.int == "spline") CV_mat1 <- rbind(CV_mat1, c(mean(RMSE_CV0a, na.rm = T), mean(RMSE_CV1a, na.rm = T), mean(RMSE_CV3a, na.rm = T)))

  print("___________________________________")
}

end <- Sys.time()
end - start

if (par.int != "spline") est.names <- c("true GLM", "naive GLM", "SIMEX", "MCEM")
if (par.int == "spline") est.names <- c("true GAM", "naive GAM", "MCEM (GAM)")
meth.len <- length(est.names)

# Bias, RMSE and variance.

if (par.int != "spline") {
  colnames(bias_mat) = colnames(RMSE_mat) = colnames(CP_mat) = est.names

  bias_mat <- stack(as.data.frame(bias_mat))
  bias_mat <- cbind(bias_mat, rep(sigma.sq.u_vec, meth.len))
  colnames(bias_mat) <- c("bias", "model", "sigma.sq")

  RMSE_mat <- stack(as.data.frame(RMSE_mat))
  RMSE_mat <- cbind(RMSE_mat, rep(sigma.sq.u_vec, meth.len))
  colnames(RMSE_mat) <- c("RMSE", "model", "sigma.sq")

  CP_mat <- stack(as.data.frame(CP_mat))
  CP_mat <- cbind(CP_mat, rep(sigma.sq.u_vec, meth.len))
  colnames(CP_mat) <- c("coverage", "model", "sigma.sq")

  my.title1 <- bquote("(a) Relative bias for" ~ beta[X])
  my.title2 <- bquote("(b) RMSE for" ~ beta[X])
  my.title3 <- bquote("(c) 95% Coverage probability for" ~ beta[X])

  p1 <- ggplot(bias_mat, aes(x = sigma.sq, y = bias, colour = model, group = model)) +
    geom_line(size = 2, aes(linetype = model)) + ggtitle(my.title1) + xlab(expression(sigma[u]^2)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
    scale_linetype_manual(values = c(1, 2, 3, 4)) +
    scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
    geom_hline(aes(yintercept = 0), colour = "grey", linetype = "dashed") +
    theme(legend.position = "none") + theme(aspect.ratio = 1)

  p2 <- ggplot(RMSE_mat, aes(x = sigma.sq, y = RMSE, colour = model, group = model)) +
    geom_line(size = 2, aes(linetype = model)) + ggtitle(my.title2) + xlab(expression(sigma[u]^2)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
    scale_linetype_manual(values = c(1, 2, 3, 4)) +
    scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
    geom_hline(aes(yintercept = 0), colour = "grey", linetype = "dashed") +
    theme(legend.position = "none") + theme(aspect.ratio = 1)

  p3 <- ggplot(CP_mat, aes(x = sigma.sq, y = coverage, colour = model, group = model)) +
    geom_line(size = 2, aes(linetype = model)) + ggtitle(my.title3) + xlab(expression(sigma[u]^2)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
    scale_linetype_manual(values = c(1, 2, 3, 4)) +
    scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
    geom_hline(aes(yintercept = 0.95), colour = "grey", linetype = "dashed") +
    theme(legend.position = "left") +
    guides(fill = guide_legend(keywidth = 1, keyheight = 1),
           linetype = guide_legend(keywidth = 7, keyheight = 2), colour = guide_legend(keywidth = 3, keyheight = 1)) +
    theme(aspect.ratio = 1)

  multiplot(p1, p2, p3, layout = matrix(c(1, 2, 3, 3), ncol = 2, nrow = 2, byrow = T))
}

# Predictive performance on test data sets (RMSE).

colnames(CV_mat1) <- est.names
CV_mat1 <- stack(as.data.frame(CV_mat1))
CV_mat1 <- cbind(CV_mat1, rep(sigma.sq.u_vec, meth.len))
colnames(CV_mat1) <- c("RMSE", "model", "sigma.sq")

if (par.int != "spline") {
  colnames(CV_mat2) <- est.names
  CV_mat2 <- stack(as.data.frame(CV_mat2))
  CV_mat2 <- cbind(CV_mat2, rep(sigma.sq.u_vec, meth.len))
  colnames(CV_mat2) <- c("RMSE", "model", "sigma.sq")
}

if (test.W == T) my.title1 <- bquote("(i) RMSE on linear predictor using the true test data" ~ paste((W[t])))
if (test.W == F) my.title1 <- bquote("(i) RMSE on linear predictor using the true test data" ~ paste((X[t])))

if (par.int != "spline") {
  my.title2 <- bquote("(ii) RMSE on linear predictor using test data with a shift")

  p1_pred <- ggplot(CV_mat1, aes(x = sigma.sq, y = RMSE, colour = model, group = model)) +
    geom_line(size = 2, aes(linetype = model)) + ggtitle(my.title1) + xlab(expression(sigma[u]^2)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
    scale_linetype_manual(values = c(1, 2, 3, 4)) + scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
    theme(legend.position = "top", legend.title = element_blank()) +
    guides(fill = guide_legend(keywidth = 1, keyheight = 2), linetype = guide_legend(keywidth = 7, keyheight = 2), colour = guide_legend(keywidth = 3, keyheight = 1))

  p2_pred <- ggplot(CV_mat2, aes(x = sigma.sq, y = RMSE, colour = model, group = model)) +
    geom_line(size = 2, aes(linetype = model)) + ggtitle(my.title2) + xlab(expression(sigma[u]^2)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
    scale_linetype_manual(values = c(1, 2, 3, 4)) + scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
    theme(legend.position = "none", legend.title = element_blank())

  x11()

  multiplot(p1_pred, p2_pred, layout = matrix(c(1, 2), ncol = 1, nrow = 2, byrow = T))
}

if (par.int == "spline") {
  my.title1 <- bquote("RMSE on linear predictor using the true test data" ~ paste((X[t])))

  p1_pred <- ggplot(CV_mat1, aes(x = sigma.sq, y = RMSE, colour = model, group = model)) +
    geom_line(size = 2, aes(linetype = model)) + ggtitle(my.title1) + xlab(expression(sigma[u]^2)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
    scale_linetype_manual(values = c(1, 2, 4)) + scale_color_manual(values = c("Red", "Black", "Blueviolet")) +
    theme(legend.position = "top", legend.title = element_blank()) + guides(fill = guide_legend(keywidth = 1, keyheight = 2), linetype = guide_legend(keywidth = 7, keyheight = 2), colour = guide_legend(keywidth = 3, keyheight = 1))

  p1_pred
}

detach(package:mgcv)

#------------------------------------------------------------------------------------

# Section 5.3: Simulation study examining model robustness.

library(refitME)

setwd(system.file(package = 'refitME'))
source('MCEM_prog.R')

set.seed(2016)

N.sim <- 200

B <- 50

n <- 1000
n_train <- 800
n_test <- n - n_train

# Set true parameter values here.

sigma.sq.e <- c(1, 1) # True predictor covariance (sigma).

alpha <- c(0, 0)
beta <- c(0.5, 1, -0.3)
par.int <- 3

scen_par <- 0.5  # Additional shift in "future" predictor values.

sigma.sq.u_vec <- c(0.0001, 0.05, 0.1, 0.25, 0.5)

bias_mat <- c()
CP_mat <- c()
SE_mat <- c()
SD_mat <- c()
RMSE_mat <- c()
CV_mat1 <- c()
CV_mat2 <- c()

family <- "binomial"

start <- Sys.time()

for(i in 1:length(sigma.sq.u_vec)) {
  nocrash <- c()

  sigma.sq.u1 <- sigma.sq.u_vec[i]  # Measurement error variance.

  beta_mat0 <- c()
  beta_mat1 <- c()
  beta_mat2 <- c()
  beta_mat3 <- c()

  beta_mat.se0 <- c()
  beta_mat.se1 <- c()
  beta_mat.se2 <- c()
  beta_mat.se3 <- c()

  RMSE_CV0a <- c()
  RMSE_CV1a <- c()
  RMSE_CV2a <- c()
  RMSE_CV3a <- c()

  RMSE_CV0b <- c()
  RMSE_CV1b <- c()
  RMSE_CV2b <- c()
  RMSE_CV3b <- c()

  c.sim <- TRUE
  k <- 0

  print("___________________________________")
  print(sigma.sq.u1)

  while(c.sim) {

    # Alternate/change here to change the true distribution of X.

    #x1 <- (rchisq(n, 3) - 3)/sqrt(6)   # Case (i). True predictor is a Chi-square distribution.
    x1 <- rsnorm(n, mean = 0, sd = 1, xi = 3) # Case (ii). True predictor is a (skewed) normal distribution.

    X <- cbind(rep(1, n), x1, x1^2)

    w1 <- x1 + rnorm(n, 0, sd = sqrt(sigma.sq.u1))  # Add/contaminate the predictor with error.

    W <- cbind(rep(1, n), w1, w1^2)
    colnames(W) <- c("Intercept", "w1", "I(w1^2)")

    eta <- X%*%beta

    mu_Y <- logit_fun(eta)
    Y <- rbinom(n, 1, prob = mu_Y)

    dat <- data.frame(cbind(Y, w1, x1))

    sigma.sq.u <- sigma.sq.u1
    sigma.sq.u_test <- sigma.sq.u1

    dat_train <- dat[1:n_train, ]
    dat_test <- data.frame(dat[(n_train + 1):n, ])
    eta_test <- eta[(n_train + 1):n]
    mu_test <- logit_fun(eta_test)

    x1_test <- dat_test$x1

    X_test <- cbind(rep(1, nrow(dat_test)), x1_test, x1_test^2)
    W_test <- cbind(rep(1, nrow(dat_test)), dat_test$w1, dat_test$w1^2)

    # Create new data (containing X, W and Y) as a future climate scenario - i.e., slightly increase X_test.

    x1_test_scen <- x1_test + scen_par

    X_test_scen <- cbind(rep(1, nrow(dat_test)), x1_test_scen, x1_test_scen^2)

    eta_test_scen <- as.matrix(X_test_scen)%*%beta
    mu_test_scen <- logit_fun(eta_test_scen)

    w1_test_scen <- x1_test_scen + rnorm(n_test, 0, sd = sqrt(sigma.sq.u_test))

    W_test_scen <- cbind(rep(1, n_test), w1_test_scen, w1_test_scen^2)

    dat_test_scen <- data.frame(cbind(w1_test_scen, x1_test_scen))
    colnames(dat_test_scen) <- colnames(dat_train)[-1]

    # Start fitting models here.

    mod_true1 <- glm(Y ~ x1 + I(x1^2), x = TRUE, family = binomial, data = dat_train)
    mod_naiv1 <- glm(Y ~ w1 + I(w1^2), x = TRUE, family = binomial, data = dat_train)

    mod_simex1 <- simex(mod_naiv1, SIMEXvariable = c("w1", "I(w1^2)"),
                          measurement.error = cbind(sqrt(sigma.sq.u), sqrt(sigma.sq.u))) # SIMEX.

    est <- try(refitME(mod_naiv1, sigma.sq.u, B), silent = TRUE)  # MCEM GLM.

    if (class(est)[1] == "try-error") next
    nocrash <- c(nocrash, 1)

    k <- k + 1
    print(k)

    beta_mat0 <- rbind(beta_mat0, mod_true1$coef)
    beta_mat1 <- rbind(beta_mat1, mod_naiv1$coef)
    beta_mat2 <- rbind(beta_mat2, mod_simex1$coef)
    beta_mat3 <- rbind(beta_mat3, est$coef)

    beta_mat.se0 <- rbind(beta_mat.se0, sqrt(diag(vcov(mod_true1))))
    beta_mat.se1 <- rbind(beta_mat.se1, sqrt(diag(vcov(mod_naiv1))))
    beta_mat.se2 <- rbind(beta_mat.se2, sqrt(diag(mod_simex1$variance.jackknife)))
    beta_mat.se3 <- rbind(beta_mat.se3, est$se)

    # CV/MSE on test data.

    eta_pred_true <- as.matrix(X_test)%*%mod_true1$coef
    RMSE_CV0a <- c(RMSE_CV0a, sqrt(mean((eta_test - eta_pred_true)^2)))
    eta_pred_true <- as.matrix(X_test_scen)%*%mod_true1$coef
    RMSE_CV0b <- c(RMSE_CV0b, sqrt(mean((eta_test_scen - eta_pred_true)^2)))

    eta_pred_naiv <- as.matrix(W_test)%*%mod_naiv1$coef
    RMSE_CV1a <- c(RMSE_CV1a, sqrt(mean((eta_test - eta_pred_naiv)^2)))
    eta_pred_naiv <- as.matrix(W_test_scen)%*%mod_naiv1$coef
    RMSE_CV1b <- c(RMSE_CV1b, sqrt(mean((eta_test_scen - eta_pred_naiv)^2)))

    eta_pred_simex <- as.matrix(W_test)%*%mod_simex1$coef
    RMSE_CV2a <- c(RMSE_CV2a, sqrt(mean((eta_test - eta_pred_simex)^2)))
    eta_pred_simex <- as.matrix(W_test_scen)%*%mod_simex1$coef
    RMSE_CV2b <- c(RMSE_CV2b, sqrt(mean((eta_test_scen - eta_pred_simex)^2)))

    eta_pred_mcem <- as.matrix(W_test)%*%est$coef
    RMSE_CV3a <- c(RMSE_CV3a, sqrt(mean((eta_test - eta_pred_mcem)^2)))
    eta_pred_mcem <- as.matrix(W_test_scen)%*%est$coef
    RMSE_CV3b <- c(RMSE_CV3b, sqrt(mean((eta_test_scen - eta_pred_mcem)^2)))

    if (length(nocrash) >= N.sim) c.sim <- FALSE
  }

  bias_mat1 <- c(mean((beta_mat0[, par.int] - beta[par.int])/beta[par.int], na.rm = T), mean((beta_mat1[, par.int] - beta[par.int])/beta[par.int], na.rm = T), mean((beta_mat2[, par.int] - beta[par.int])/beta[par.int], na.rm = T), mean((beta_mat3[, par.int] - beta[par.int])/beta[par.int], na.rm = T))
  RMSE_mat1 <- c(sqrt(mean((beta_mat0[, par.int] - beta[par.int])^2, na.rm = T)), sqrt(mean((beta_mat1[, par.int] - beta[par.int])^2, na.rm = T)), sqrt(mean((beta_mat2[, par.int] - beta[par.int])^2, na.rm = T)), sqrt(mean((beta_mat3[, par.int] - beta[par.int])^2, na.rm = T)))
  CP_mat_SE3 <- which(is.nan(beta_mat.se3[, par.int]) | is.na(beta_mat.se3[, par.int]))

  if (length(CP_mat_SE3) == 0) CP_mat3 <- cov.pc(beta, beta_mat3, beta_mat.se3, N.sim - length(CP_mat_SE3))[par.int]
  if (length(CP_mat_SE3) > 0) CP_mat3 <- cov.pc(beta, beta_mat3[-CP_mat_SE3, ], beta_mat.se3[-CP_mat_SE3, ], N.sim - length(CP_mat_SE3))[par.int]
  CP_mat1 <- c(cov.pc(beta, beta_mat0, beta_mat.se0, N.sim)[par.int], cov.pc(beta, beta_mat1, beta_mat.se1, N.sim)[par.int], cov.pc(beta, beta_mat2, beta_mat.se2, N.sim)[par.int], CP_mat3)

  bias_mat <- rbind(bias_mat, bias_mat1)
  RMSE_mat <- rbind(RMSE_mat, RMSE_mat1)
  CP_mat <- rbind(CP_mat, CP_mat1)

  CV_mat1 <- rbind(CV_mat1, c(mean(RMSE_CV0a, na.rm = T), mean(RMSE_CV1a, na.rm = T), mean(RMSE_CV2a, na.rm = T), mean(RMSE_CV3a, na.rm = T)))
  CV_mat2 <- rbind(CV_mat2, c(mean(RMSE_CV0b, na.rm = T), mean(RMSE_CV1b, na.rm = T), mean(RMSE_CV2b, na.rm = T), mean(RMSE_CV3b, na.rm = T)))

  print("___________________________________")
}

end <- Sys.time()
end - start

est.names <- c("true GLM", "naive GLM", "SIMEX", "MCEM")
meth.len <- length(est.names)

# Bias, RMSE and variance.

colnames(bias_mat) = colnames(RMSE_mat) = colnames(CP_mat) = est.names

bias_mat <- stack(as.data.frame(bias_mat))
bias_mat <- cbind(bias_mat, rep(sigma.sq.u_vec, meth.len))
colnames(bias_mat) <- c("bias", "model", "sigma.sq")

RMSE_mat <- stack(as.data.frame(RMSE_mat))
RMSE_mat <- cbind(RMSE_mat, rep(sigma.sq.u_vec, meth.len))
colnames(RMSE_mat) <- c("RMSE", "model", "sigma.sq")

CP_mat <- stack(as.data.frame(CP_mat))
CP_mat <- cbind(CP_mat, rep(sigma.sq.u_vec, meth.len))
colnames(CP_mat) <- c("coverage", "model", "sigma.sq")

my.title1 <- bquote("(a) Relative bias for" ~ beta[X])
my.title2 <- bquote("(b) RMSE for" ~ beta[X])
my.title3 <- bquote("(c) 95% Coverage probability for" ~ beta[X])

p1 <- ggplot(bias_mat, aes(x = sigma.sq, y = bias, colour = model, group = model)) +
    geom_line(size = 2, aes(linetype = model)) + ggtitle(my.title1) + xlab(expression(sigma[u]^2)) +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
    scale_linetype_manual(values = c(1, 2, 3, 4)) + scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
    geom_hline(aes(yintercept = 0), colour = "grey", linetype = "dashed") + theme(legend.position = "non") + theme(aspect.ratio = 1)

p2 <- ggplot(RMSE_mat, aes(x = sigma.sq, y = RMSE, colour = model, group = model)) + geom_line(size = 2, aes(linetype = model)) +
    ggtitle(my.title2) + xlab(expression(sigma[u]^2)) + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
    scale_linetype_manual(values = c(1, 2, 3, 4)) + scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
    geom_hline(aes(yintercept = 0), colour = "grey", linetype = "dashed") + theme(legend.position = "none") + theme(aspect.ratio = 1)

p3 <- ggplot(CP_mat, aes(x = sigma.sq, y = coverage, colour = model, group = model)) + geom_line(size = 2, aes(linetype = model)) +
    ggtitle(my.title3) + xlab(expression(sigma[u]^2)) + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
    scale_linetype_manual(values = c(1, 2, 3, 4)) + scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
    geom_hline(aes(yintercept = 0.95), colour = "grey", linetype = "dashed") + theme(legend.position = "left") +
    guides(fill = guide_legend(keywidth = 1, keyheight = 1), linetype = guide_legend(keywidth = 7, keyheight = 2), colour = guide_legend(keywidth = 3, keyheight = 1)) +
    theme(aspect.ratio = 1)

multiplot(p1, p2, p3, layout = matrix(c(1, 2, 3, 3), ncol = 2, nrow = 2, byrow = T))

#-----------------------------------------------------------------------------------------

# Web Appendix A.3: Capture-recapture simulation study.

library(refitME)
library(VGAM)

setwd(system.file(package = 'refitME'))
source('MCEM_prog.R')

set.seed(2016)

N.sim <- 200

epsilon <- 0.00001
B <- 100

N <- 200
tau <- 7

# Set true parameter values here.

beta <- c(-1, 1);  par.int <- 2

# Un-hash/toggle below if you want to fit quadratic models.

beta <- c(0.5, 0, -0.3); par.int <- 3

sigma.sq.e <- 1

sigma.sq.u_vec <- c(0.0001, 0.05, 0.1, 0.25, 0.5)

bias_mat <- c()
CP_mat <- c()
SE_mat <- c()
SD_mat <- c()
RMSE_mat <- c()
bias_N_mat <- c()
CP_N_mat <- c()
SE_N_mat <- c()
SD_N_mat <- c()
RMSE_N_mat <- c()

start <- Sys.time()

for(i in 1:length(sigma.sq.u_vec)) {
  nocrash <- c()

  sigma.sq.u <- sigma.sq.u_vec[i]  # Measurement error variance.

  beta_mat0 = beta_mat1 = beta_mat2 = beta_mat3 <- c()
  beta_mat.se0 = beta_mat.se1 = beta_mat.se2 = beta_mat.se3 <- c()

  N_mat0 = N_mat1 = N_mat2 = N_mat3 <- c()
  N_mat.se0 = N_mat.se1 = N_mat.se2 = N_mat.se3 <- c()

  c.sim <- TRUE
  k <- 0

  print("___________________________________")
  print(sigma.sq.u)

  while(c.sim) {
    dats <- genCRdata(N, beta, tau, sigma.sq.u, par.int)

    cap.hist <- dats$hist
    X <- dats$X
    W <- dats$W
    y <- dats$y
    x1 <- X[, 2]
    w1 <- W[, 2]

    CR_dat0 <- data.frame(cbind(x1, y, tau - y))
    CR_dat1 <- data.frame(cbind(w1, y, tau - y))
    colnames(CR_dat0) <- c("x1", "cap", "noncap")
    colnames(CR_dat1) <- c("w1", "cap", "noncap")

    # Start fitting models here (model M_h).

    # True and naive models.

    if (par.int == 2) {
      mod_true1 <- vglm(cbind(cap, noncap) ~ x1, posbinomial(omit.constant = TRUE, parallel = TRUE ~ x1), data = CR_dat0, trace = F)
      mod_naiv1 <- vglm(cbind(cap, noncap) ~ w1, posbinomial(omit.constant = TRUE, parallel = TRUE ~ w1), data = CR_dat1, trace = F)
      mod_naiv_CS <- mod_naiv1
    }

    if (par.int == 3) {
      mod_true1 <- vglm(cbind(cap, noncap) ~ x1 + I(x1^2), posbinomial(omit.constant = TRUE, parallel = TRUE ~ x1 + I(x1^2)), data = CR_dat0, trace = F)
      mod_naiv1 <- vglm(cbind(cap, noncap) ~ w1 + I(w1^2), posbinomial(omit.constant = TRUE, parallel = TRUE ~ w1 + I(w1^2)), data = CR_dat1, trace = F)
      mod_naiv_CS <- vglm(cbind(cap, noncap) ~ w1, posbinomial(omit.constant = TRUE, parallel = TRUE ~ w1), data = CR_dat1, trace = F)
    }

    # Conditional score method.

    CS_beta.est <- nleqslv(coef(mod_naiv_CS), est.cs, y = y, w1 = w1, tau = tau, sigma.sq.u = sigma.sq.u, method = c("Newton"))$x
    CS_N.est <- N.CS.est(CS_beta.est, y, w1, tau, sigma.sq.u)

    var.ests1 <- var.CS(CS_beta.est, y, w1, tau, sigma.sq.u)
    if (par.int == 2) {
      CS_beta.est.se <- sqrt(var.ests1)[1:(par.int)]
      CS_N.est.se <- sqrt(var.ests1)[par.int + 1]
    }

    if (par.int == 3) {
      CS_beta.est.se <- sqrt(var.ests1)[1:(par.int- 1)]
      CS_N.est.se <- sqrt(var.ests1)[par.int]
    }

    # MCEM.

    est <- try(refitME(mod_naiv1, sigma.sq.u, B), silent = TRUE)
    if (class(est) == "try-error") next
    nocrash <- c(nocrash, 1)

    N.hat_MCEM <- est$N.est
    N.hat.se_MCEM <- est$N.est.se

    k <- k + 1
    print(k)

    beta_mat0 <- rbind(beta_mat0, coef(mod_true1))
    beta_mat1 <- rbind(beta_mat1, coef(mod_naiv1))
    if (par.int == 2) beta_mat2 <- rbind(beta_mat2, CS_beta.est)
    if (par.int == 3) beta_mat2 <- rbind(beta_mat2, c(CS_beta.est, 0))
    beta_mat3 <- rbind(beta_mat3, c(est$beta))

    beta_mat.se0 <- rbind(beta_mat.se0, sqrt(diag(vcov(mod_true1))))
    beta_mat.se1 <- rbind(beta_mat.se1, sqrt(diag(vcov(mod_naiv1))))
    if (par.int == 2) beta_mat.se2 <- rbind(beta_mat.se2, CS_beta.est.se)
    if (par.int == 3) beta_mat.se2 <- rbind(beta_mat.se2, c(CS_beta.est.se, 0))
    beta_mat.se3 <- rbind(beta_mat.se3, c(est$beta.se2))

    N_mat0 <- c(N_mat0, mod_true1@extra$N.hat)
    N_mat1 <- c(N_mat1, mod_naiv1@extra$N.hat)
    N_mat2 <- c(N_mat2, CS_N.est)
    N_mat3 <- c(N_mat3, N.hat_MCEM)

    N_mat.se0 <- rbind(N_mat.se0, mod_true1@extra$SE.N.hat)
    N_mat.se1 <- rbind(N_mat.se1, mod_naiv1@extra$SE.N.hat)
    N_mat.se2 <- rbind(N_mat.se2, CS_N.est.se)
    N_mat.se3 <- rbind(N_mat.se3, N.hat.se_MCEM)

    if (length(nocrash) >= N.sim) c.sim <- FALSE
  }

  bias_mat1 <- c(mean((beta_mat0[, par.int] - beta[par.int])/beta[par.int], na.rm = T), mean((beta_mat1[, par.int] - beta[par.int])/beta[par.int], na.rm = T), mean((beta_mat2[, par.int] - beta[par.int])/beta[par.int], na.rm = T), mean((beta_mat3[, par.int] - beta[par.int])/beta[par.int], na.rm = T))
  RMSE_mat1 <- c(sqrt(mean((beta_mat0[, par.int] - beta[par.int])^2, na.rm = T)), sqrt(mean((beta_mat1[, par.int] - beta[par.int])^2, na.rm = T)), sqrt(mean((beta_mat2[, par.int] - beta[par.int])^2, na.rm = T)), sqrt(mean((beta_mat3[, par.int] - beta[par.int])^2, na.rm = T)))
  CP_mat_SE3 <- which(is.nan(beta_mat.se3[, par.int]) | is.na(beta_mat.se3[, par.int]))
  if (length(CP_mat_SE3) == 0) CP_mat3 <- cov.pc(beta, beta_mat3, beta_mat.se3, N.sim - length(CP_mat_SE3))[par.int]
  if (length(CP_mat_SE3) > 0) CP_mat3 <- cov.pc(beta, beta_mat3[-CP_mat_SE3, ], beta_mat.se3[-CP_mat_SE3, ], N.sim - length(CP_mat_SE3))[par.int]
  CP_mat1 <- c(cov.pc(beta, beta_mat0, beta_mat.se0, N.sim)[par.int], cov.pc(beta, beta_mat1, beta_mat.se1, N.sim)[par.int], cov.pc(beta, beta_mat2, beta_mat.se2, N.sim)[par.int], CP_mat3)

  bias_mat2 <- c(mean((N_mat0 - N)/N, na.rm = T), mean((N_mat1 - N)/N, na.rm = T), mean((N_mat2 - N)/N, na.rm = T), mean((N_mat3 - N)/N, na.rm = T))
  RMSE_mat2 <- c(sqrt(mean((N_mat0 - N)^2, na.rm = T)), sqrt(mean((N_mat1 - N)^2, na.rm = T)), sqrt(mean((N_mat2 - N)^2, na.rm = T)), sqrt(mean((N_mat3 - N)^2, na.rm = T)))
  CP_mat_SE4 <- which(is.nan(N_mat.se3) | is.na(N_mat.se3))
  if (length(CP_mat_SE4) == 0) CP_mat4 <- cov.pc_N(N, N_mat3, N_mat.se3, N.sim - length(CP_mat_SE4))
  if (length(CP_mat_SE4) > 0) CP_mat4 <- cov.pc_N(N, N_mat3[-CP_mat_SE4, ], N_mat.se3[-CP_mat_SE4, ], N.sim - length(CP_mat_SE4))
  CP_mat2 <- c(cov.pc_N(N, N_mat0, N_mat.se0, N.sim), cov.pc_N(N, N_mat1, N_mat.se1, N.sim), cov.pc_N(N, N_mat2, N_mat.se2, N.sim), CP_mat4)

  bias_mat <- rbind(bias_mat, bias_mat1)
  RMSE_mat <- rbind(RMSE_mat, RMSE_mat1)
  CP_mat <- rbind(CP_mat, CP_mat1)

  bias_N_mat <- rbind(bias_N_mat, bias_mat2)
  RMSE_N_mat <- rbind(RMSE_N_mat, RMSE_mat2)
  CP_N_mat <- rbind(CP_N_mat, CP_mat2)

  print("___________________________________")
}

end <- Sys.time()
end - start

est.names <- c("true CL", "naive CL", "CS", "MCEM (CL)")
meth.len <- length(est.names)

# Bias, RMSE and variance.

colnames(bias_mat) = colnames(RMSE_mat) = colnames(CP_mat) = colnames(bias_N_mat) = colnames(RMSE_N_mat) = colnames(CP_N_mat) = est.names

bias_mat <- stack(as.data.frame(bias_mat))
bias_mat <- cbind(bias_mat, rep(sigma.sq.u_vec, meth.len))
colnames(bias_mat) <- c("bias", "model", "sigma.sq")

RMSE_mat <- stack(as.data.frame(RMSE_mat))
RMSE_mat <- cbind(RMSE_mat, rep(sigma.sq.u_vec, meth.len))
colnames(RMSE_mat) <- c("RMSE", "model", "sigma.sq")

CP_mat <- stack(as.data.frame(CP_mat))
CP_mat <- cbind(CP_mat, rep(sigma.sq.u_vec, meth.len))
colnames(CP_mat) <- c("coverage", "model", "sigma.sq")

bias_N_mat <- stack(as.data.frame(bias_N_mat))
bias_N_mat <- cbind(bias_N_mat, rep(sigma.sq.u_vec, meth.len))
colnames(bias_N_mat) <- c("bias", "model", "sigma.sq")

RMSE_N_mat <- stack(as.data.frame(RMSE_N_mat))
RMSE_N_mat <- cbind(RMSE_N_mat, rep(sigma.sq.u_vec, meth.len))
colnames(RMSE_N_mat) <- c("RMSE", "model", "sigma.sq")

CP_N_mat <- stack(as.data.frame(CP_N_mat))
CP_N_mat <- cbind(CP_N_mat, rep(sigma.sq.u_vec, meth.len))
colnames(CP_N_mat) <- c("coverage", "model", "sigma.sq")

my.title1 <- bquote("(a) Relative bias for" ~ beta[X])
my.title2 <- bquote("(b) RMSE for" ~ beta[X])
my.title3 <- bquote("(c) 95% Coverage probability for" ~ beta[X])
my.title4 <- bquote("(a) Relative bias for" ~ N)
my.title5 <- bquote("(b) RMSE for" ~ N)
my.title6 <- bquote("(c) 95% Coverage probability for" ~ N)

p1 <- ggplot(bias_mat, aes(x = sigma.sq, y = bias, colour = model, group = model)) +
  geom_line(size = 2, aes(linetype = model)) + ggtitle(my.title1) + xlab(expression(sigma[u]^2)) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
  scale_linetype_manual(values = c(1, 2, 3, 4)) + scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
  geom_hline(aes(yintercept = 0), colour = "grey", linetype = "dashed") + theme(legend.position = "none") + theme(aspect.ratio = 1)

p2 <- ggplot(RMSE_mat, aes(x = sigma.sq, y = RMSE, colour = model, group = model)) + geom_line(size = 2, aes(linetype = model)) +
  ggtitle(my.title2) + xlab(expression(sigma[u]^2)) + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
  scale_linetype_manual(values = c(1, 2, 3, 4)) + scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
  geom_hline(aes(yintercept = 0), colour = "grey", linetype = "dashed") + theme(legend.position = "none") + theme(aspect.ratio = 1)

p3 <- ggplot(CP_mat, aes(x = sigma.sq, y = coverage, colour = model, group = model)) + geom_line(size = 2, aes(linetype = model)) +
  ggtitle(my.title3) + xlab(expression(sigma[u]^2)) + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
  scale_linetype_manual(values = c(1, 2, 3, 4)) + scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
  geom_hline(aes(yintercept = 0.95), colour = "grey", linetype = "dashed") + theme(legend.position = "left") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1), linetype = guide_legend(keywidth = 7, keyheight = 2), colour = guide_legend(keywidth = 3, keyheight = 1)) +
  theme(aspect.ratio = 1)

x11()

p4 <- ggplot(bias_N_mat, aes(x = sigma.sq, y = bias, colour = model, group = model)) +
  geom_line(size = 2, aes(linetype = model)) + ggtitle(my.title4) + xlab(expression(sigma[u]^2)) +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
  scale_linetype_manual(values = c(1, 2, 3, 4)) + scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
  geom_hline(aes(yintercept = 0), colour = "grey", linetype = "dashed") + theme(legend.position = "none") + theme(aspect.ratio = 1)

p5 <- ggplot(RMSE_N_mat, aes(x = sigma.sq, y = RMSE, colour = model, group = model)) + geom_line(size = 2, aes(linetype = model)) +
  ggtitle(my.title5) + xlab(expression(sigma[u]^2)) + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
  scale_linetype_manual(values = c(1, 2, 3, 4)) + scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
  geom_hline(aes(yintercept = 0), colour = "grey", linetype = "dashed") + theme(legend.position = "none") + theme(aspect.ratio = 1)

p6 <- ggplot(CP_N_mat, aes(x = sigma.sq, y = coverage, colour = model, group = model)) + geom_line(size = 2, aes(linetype = model)) +
  ggtitle(my.title6) + xlab(expression(sigma[u]^2)) + theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16, face = "bold")) +
  scale_linetype_manual(values = c(1, 2, 3, 4)) + scale_color_manual(values = c("Red", "Black", "Forestgreen", "Blueviolet")) +
  geom_hline(aes(yintercept = 0.95), colour = "grey", linetype = "dashed") + theme(legend.position = "left") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1), linetype = guide_legend(keywidth = 7, keyheight = 2), colour = guide_legend(keywidth = 3, keyheight = 1)) +
  theme(aspect.ratio = 1)

multiplot(p4, p5, p6, layout = matrix(c(1, 2, 3, 3), ncol = 2, nrow = 2, byrow = T))

detach(package:VGAM)
