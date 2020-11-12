#-----------------------------------------------------------------------
# refitME_package.R
#
# Source file for the "refitME" R-package.
#-----------------------------------------------------------------------

# Required R-packages.

library(mvtnorm)
library(MASS)
library(mgcv)
library(Matrix)
library(VGAM)
library(expm)
suppressMessages(library(sandwich))

#' @title The Framingham heart study data set.
#' @description Data set consisting of records of male patients with coronary heart disease collected from the Framingham heart study. The \code{Framinghamdata} data consists of binary responses and four predictor variables collected on `n = 1615` patients.
#' @format A data set that contains: 5 columns with 1,615 observations. The columns are defined as follows:
#' \describe{
#' \item{\code{Y}}{Response indicator (binary variable) of first evidence of CHD status of patient.}
#' \item{\code{z1}}{Serum cholesterol level of patient.}
#' \item{\code{z2}}{Age of patient.}
#' \item{\code{z3}}{Smoking indicator - whether the patient smokes.}
#' \item{\code{w1}}{Systolic blood pressure (SBP) of patient - this is the error contaminated variable, calculated from mean scores. The measurement error is 0.00630, see pp. 112 of Carroll \emph{et al.} (2006).}
#' }
#' @source See Carroll \emph{et al.} (2006) for full details of the data and study. Also, see \url{https://github.com/JakubStats/refitME} for an RMarkdown tutorial of an example that uses the data.
#' @references Carroll, R. J., Ruppert, D., Stefanski, L. A., and Crainiceanu, C. M. (2006). \emph{Measurement Error in Nonlinear Models: A Modern Perspective.} 2nd Ed. London: Chapman \& Hall/CRC.
#' @examples # Load the data.
#'
#' data(Framinghamdata)
"Framinghamdata"

#' MCEMfit_glm
#'
#' Function for wrapping the MCEM algorithm on GLMs where covariates are subject to measurement error/error-in-variables.
#' @name MCEMfit_glm
#' @param mod : a glm object (this is the naive fitted model). Make sure the first input predictor variables are the selected error-contaminated variable (i.e., the \code{W}'s).
#' @param family : a specified family/distribution.
#' @param sigma.sq.u : measurement error variance. A scalar if there is only one error-contaminated variable, otherwise this must stored as a covariance matrix.
#' @param W : a matrix of error-contaminated covariates (if not specified, the default assumes all covariates in the naive fitted model are error-contaminated).
#' @param sigma.sq.e : variance of the true covariate (X).
#' @param B : the number of Monte Carlo replication values (default is set to 50).
#' @param epsilon : a set convergence threshold (default is set to 0.00001).
#' @param theta.est : an initial value for the dispersion parameter (this is required for fitting negative binomial models).
#' @param shape.est : an initial value for the shape parameter (this is required for fitting gamma models).
#' @return \code{MCEMfit_glm} returns the naive fitted model object where coefficient estimates, the covariance matrix, fitted values and residuals have been replaced with the final MCEM model fit. Standard errors and effective sample size have been additional included.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. \pkg{refitME}: Measurement Error Modelling using Monte Carlo Expectation Maximization in \proglang{R}.
#' @import mvtnorm MASS mgcv sandwich expm
#' @importFrom stats Gamma
#' @export
#' @seealso \code{\link{MCEMfit_gam}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown tutorial with examples.
MCEMfit_glm <- function(mod, family, sigma.sq.u, W = NULL, sigma.sq.e = 1, B = 50, epsilon = 0.00001, theta.est = 1, shape.est = 1) {

  options(warn = -1)

  mod_n <- mod

  reps <- 0
  cond <- TRUE

  mod.terms <- attr(mod$terms, "term.labels")

  d <- length(mod.terms)

  Y <- mod$model[, 1]

  bigY <- rep(Y, B)

  n <- length(Y)

  if (family == "gaussian") {
    if (is.null(stats::weights(mod))) p.wt <- rep(1, n)
    if (!is.null(stats::weights(mod))) p.wt <- stats::weights(mod)
  }
  if (family != "gaussian") p.wt <- stats::weights(mod)

  bigp.wt <- rep(p.wt, B)

  if(length(names(mod$model)) > (d + 1))  W1 <- as.matrix(mod$model[, -c(1, (d + 2))])
  if(length(names(mod$model)) == (d + 1))  W1 <- as.matrix(mod$model[, -1])

  if(is.null(W)) W <- W1

  muPred <- rep(stats::predict(mod, type = "response"), B)

  beta.est <- stats::coef(mod)

  if (is.matrix(sigma.sq.u) == F) {
    w1 <- mod$model[, 2]

    if (d == 1) p1 <- 1

    if (d > 1) {
      if (sum((w1)^2 != (W1[, 2])) == n) p1 <- 1
      if (sum((w1)^2 == (W1[, 2])) == n) p1 <- 2
    }

    sigma.sq.u1 <- sigma.sq.u
    sigma.sq.e1 <- sigma.sq.e

    U1_j <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1, B)))
    X1_j <- rep(w1, B) - U1_j

    X <- cbind(rep(1, B*n), X1_j)
    if (p1 == 2) X <- cbind(rep(1, B*n), X1_j, (X1_j)^2)

    if (!is.null(ncol(W1))) {
      if (d > 2 & p1 == 2) X <- cbind(X, do.call(rbind, replicate(B, W1[, -c(1:p1)], simplify = FALSE)))
      if (d > 2 & p1 == 1) X <- cbind(X, do.call(rbind, replicate(B, W1[, -c(1:p1)], simplify = FALSE)))
      if (d == 2 & p1 == 2) X <- X
    }

    X <- as.matrix(X)

    colnames(X) <- names(stats::coef(mod))

    colnames(X)[2] <- "x1"

    if (p1 == 2) colnames(X)[2:3] <- c("x1", "I(x1^2)")

    mu.e1 <- mean(X1_j)
  }

  if (is.matrix(sigma.sq.u) == T) {
    q1 <- dim(W)[2]

    if (d == q1) {
      p <- rep(1, q1)
      col.nameX <- paste0('x', 1:q1)

      XA <- W
    }
    if (d > q1) {
      p <- c()
      col.nameX <- c()
      XA <- c()

      kk1 <- 1

      for(kk in 1:q1) {
        if (sum((W[, kk])^2 != (W1[, kk1 + 1])) == n) {
          p1 <- 1
          kk1 <- kk
          col.nameX1 <- paste0('x', kk)
        }

        if (sum((W[, kk])^2 == (W1[, kk1 + 1])) == n) {
          p1 <- 2
          kk1 <- kk + 1
          col.nameX1 <- c(paste0('x', kk), paste0(paste0('I(x', 1), '^2)'))
        }

        p <- c(p, p1)
        col.nameX <- c(col.nameX, col.nameX1)
      }
    }

    if (d < q1) print("Error: Number of error-contaminated covariates exceeds total number of covariates!")

    p2 <- sum(p)
    sigma.sq.u1 <- sigma.sq.u
    sigma.sq.e1 <- diag(sigma.sq.e)

    w1 <- W[, 1]
    U1_j <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1[1, 1], B)))
    X1_j <- rep(w1, B) - U1_j

    X <- cbind(rep(1, B*n), X1_j)
    XA <- X[, 2]
    if (p[1] == 2) X <- cbind(rep(1, B*n), X1_j, (X1_j)^2)

    mu.e1 <- mean(X1_j)

    for(kk in 2:q1) {
      w1 <- W[, kk]
      U1_j <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1[kk, kk], B)))
      X1_j <- rep(w1, B) - U1_j

      XA <- cbind(XA, X1_j)

      Xa <- X1_j
      if (p[kk] == 2) Xa <- cbind(Xa, (X1_j)^2)
      X <- cbind(X, Xa)
      mu.e1 <- c(mu.e1, mean(X1_j))
    }

    if (d == p2) X <- X
    if (d > p2) X <- cbind(X, do.call(rbind, replicate(B, as.matrix(W1[, -c(1:p2)]), simplify = FALSE)))

    X <- as.matrix(X)

    colnames(X) <- names(stats::coef(mod))

    colnames(X)[2:(p2 + 1)] <- col.nameX
  }

  while(cond) {

    # MC and E-step.

    if (is.matrix(sigma.sq.u) == F) prX <- stats::dnorm(X1_j, mu.e1, sd = sqrt(sigma.sq.e1))
    if (is.matrix(sigma.sq.u) == T) prX <- mvtnorm::dmvnorm(XA, mu.e1, sigma = sqrt(sigma.sq.e1))

    if (family == "gaussian") prY <- stats::dnorm(bigY, muPred, 1)
    if (family == "binomial") prY <- stats::dbinom(bigY, 1, muPred)
    if (family == "poisson") prY <- stats::dpois(bigY, muPred)
    if (family == "Gamma") prY <- stats::dgamma(bigY, shape = shape.est, scale = muPred/shape.est)
    if (family == "negbin") prY <- stats::dnbinom(bigY, size = theta.est, mu = muPred)

    prY[prY == 0] <- 1.e-6

    bigW <- matrix(prY*prX, n, B)
    sumW <- rep(apply(bigW, 1, sum, na.rm = T), B)
    weights1 <- as.vector(bigW)/sumW

    weights2 <- weights1
    weights1 <- bigp.wt*weights1

    # M-step (updates).

    if (family == "gaussian") {
      mod <- stats::lm(bigY ~ X - 1, weights = weights1)
      sigma.sq.est <- (summary(mod)$sigma)^2
    }
    if (family == "binomial") mod <- stats::glm(bigY ~ X - 1, weights = weights1, family = "binomial")
    if (family == "poisson") mod <- stats::glm(bigY ~ X - 1, weights = weights1, family = "poisson")
    if (family == "Gamma") mod <- stats::glm(bigY ~ X - 1, weights = weights1, family = Gamma(link = "log"))
    if (family == "negbin") mod <- MASS::glm.nb(bigY ~ X - 1, weights = weights1, init.theta = theta.est)

    beta.update <- stats::coef(mod)

    if (family == "negbin") theta.update <- mod$theta

    muPred <- stats::predict(mod, type = "response")

    if (family == "Gamma") shape.update <- summary(mod)[14]$dispersion

    if (is.matrix(sigma.sq.u) == F) {
      sigma.sq.e1.update <- wt.var(X[, 2], w = weights2)
      mu.e1.update <- stats::weighted.mean(X[, 2], w = weights2)
    }

    if (is.matrix(sigma.sq.u) == T) {
      sigma.sq.e1.update <- c()
      mu.e1.update <- c()

      for(kk in 1:q1) {
        sigma.sq.e1.update1 <- wt.var(XA[, kk], w = weights2)
        sigma.sq.e1.update <- c(sigma.sq.e1.update, sigma.sq.e1.update1)

        mu.e1.update1 <- stats::weighted.mean(XA[, kk], w = weights2)
        mu.e1.update <- c(mu.e1.update, mu.e1.update1)
      }

      sigma.sq.e1.update <- diag(sigma.sq.e1.update)
    }

    # Convergence monitoring.

    beta.norm <- sum((beta.est - beta.update)^2)

    if (family == "negbin") theta.norm <- sum((theta.est - theta.update)^2)
    if (family == "Gamma") shape.norm <- sum((shape.est - shape.update)^2)

    if (is.matrix(sigma.sq.u) == F) diff.sig_e <- abs(sigma.sq.e1.update - sigma.sq.e1)
    if (is.matrix(sigma.sq.u) == T) diff.sig_e <- sum(abs(diag(sigma.sq.e1.update) - diag(sigma.sq.e1)))

    diff.mu_e <- sum((mu.e1.update - mu.e1)^2)

    reps <- reps + 1 # Keeps track of number of iterations.

    if (family == "binomial" | family == "poisson" | family == "gaussian") {
      if (diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon) {
        cond <- FALSE
        print("convergence :-)")
        print(reps)
        break
      }
    }

    if (family == "negbin") {
      if (diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon & theta.norm < epsilon) {
        cond <- FALSE
        print("convergence :-)")
        print(reps)
        break
      }
    }

    if (family == "Gamma") {
      if (diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon & shape.norm < epsilon) {
        cond <- FALSE
        print("convergence :-)")
        print(reps)
        break
      }
    }

    # Update parameters.

    beta.est <- beta.update
    if (family == "negbin") theta.est <- theta.update
    if (family == "Gamma") shape.est <- shape.update
    sigma.sq.e1 <- sigma.sq.e1.update
    mu.e1 <- mu.e1.update
  }

  mod_n$coefficients <- beta.est

  eta <- stats::predict(mod_n)

  if (family == "gaussian") {
    residuals <- Y - stats::fitted(mod_n)
    mod_n$residuals <- residuals
  }

  if (family != "gaussian") {
    mu <- mod_n$family$linkinv(eta)
    mu.eta <- mod_n$family$mu.eta(eta)
    residuals <- (Y - mu)/mod_n$family$mu.eta(eta)
    mod_n$residuals <- residuals
  }

  mod_n$linear.predictors <- stats::predict(mod_n)

  mod_n$fitted.values <- stats::predict(mod_n, type = "response")

  sumW <- apply(bigW, 1, sum, na.rm = T)
  weights1 <- bigW/sumW
  weights1[is.nan(weights1)] <- 0

  eff.samp.size <- 1/apply(weights1^2, 1, sum)
  eff.samp.size[is.infinite(eff.samp.size)] <- "NA"
  eff.samp.size <- as.numeric(eff.samp.size)

  mod_n$eff.samp.size <- eff.samp.size

  # Standard error calculations start here.

  estfun_mat <- sandwich::estfun(mod)
  if (family == "gaussian") estfun_mat <- sandwich::estfun(mod)*B

  K1 <- length(beta.update)

  S_1 <- matrix(0, nrow = K1, ncol = K1)
  SS_1 <- matrix(0, nrow = K1, ncol = K1)

  ind_mat <- matrix(1:(n*B), ncol = n, byrow = T)

  for(ii in 1:n) {
    index_vec <- ind_mat[, ii]
    S_1 <- S_1 + (apply(estfun_mat[index_vec, ], 2, sum))%*%t(apply(estfun_mat[index_vec, ], 2, sum))
  }

  if (family == "gaussian") {
    sand1 <- (sandwich::estfun(mod)*B/mod$weights)
    sand1[is.nan(sand1)] <- 1
    SS_1 <- t(sandwich::estfun(mod)*B)%*%sand1
    u.bar <- solve(stats::vcov(mod)*B)
    beta.est.se2 <- sqrt(diag(solve(u.bar - SS_1/B^2 + S_1/B^2)))
    AA <- expm::sqrtm((u.bar - SS_1/B^2 + S_1/B^2)*(stats::summary.lm(mod_n)$sigma)^2)
  }

  if (family != "gaussian") {
    sand1 <- (sandwich::estfun(mod)/mod$prior)
    sand1[is.nan(sand1)] <- 1
    SS_1 <- t(sandwich::estfun(mod))%*%sand1
    u.bar <- solve(stats::vcov(mod))
    beta.est.se2 <- sqrt(diag(solve(u.bar - SS_1 + S_1)))
    AA <- expm::sqrtm((u.bar - SS_1 + S_1)*stats::summary.glm(mod_n)$dispersion)
  }

  if (length(which(is.nan(beta.est.se2))) > 0) beta.est.se2 <- c(rep(NA, K1))

  mod_n$se <- beta.est.se2

  mod_n$qr <- qr(AA)
  dim_temp <- ncol(mod_n$qr$qr)*nrow(mod_n$qr$qr) - ncol(mod_n$qr$qr)*ncol(mod_n$qr$qr)
  mod_n$qr$qr <- rbind(qr(AA)$qr, matrix(rep(0, dim_temp), ncol = ncol(mod_n$qr$qr)))

  return(mod_n)
}

#' MCEMfit_gam
#'
#' Function for wrapping the MCEM algorithm on GAMs where covariates are subject to measurement error/error-in-variables.
#' @name MCEMfit_gam
#' @param mod : a gam object (this is the naive fitted model). Make sure the first input predictor variables are the selected error-contaminated variable (i.e., the \code{W}'s).
#' @param family : a specified family/distribution.
#' @param sigma.sq.u : measurement error variance. A scalar if there is only one error-contaminated variable, otherwise this must stored as a covariance matrix.
#' @param W : a matrix of error-contaminated covariates (if not specified, the default assumes all covariates in the naive fitted model are error-contaminated).
#' @param sigma.sq.e : variance of the true covariate (X).
#' @param B : the number of Monte Carlo replication values (default is set to 50).
#' @param epsilon : convergence threshold (default is set to 0.00001).
#' @param theta.est : an initial value for the dispersion parameter (this is required for fitting negative binomial models).
#' @param shape.est : an initial value for the shape parameter (this is required for fitting gamma models).
#' @return \code{MCEMfit_gam} returns the naive fitted model object where coefficient estimates, the covariance matrix, fitted values and residuals have been replaced with the final MCEM model fit. Standard errors and effective sample size have been additional included.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. \pkg{refitME}: Measurement Error Modelling using Monte Carlo Expectation Maximization in \proglang{R}.
#' @import mvtnorm MASS mgcv sandwich
#' @importFrom stats Gamma
#' @export
#' @seealso \code{\link{MCEMfit_glm}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown tutorial with examples.
MCEMfit_gam <- function(mod, family, sigma.sq.u, W = NULL, sigma.sq.e = 1, B = 50, epsilon = 0.00001, theta.est = 1, shape.est = 10) {

  options(warn = -1)

  reps <- 0
  cond <- TRUE

  mod.terms <- attr(mod$terms, "term.labels")

  d <- length(mod.terms)

  Y <- mod$model[, 1]

  bigY <- rep(Y, B)

  n <- length(Y)

  if (family == "gaussian") {
    if (is.null(stats::weights(mod))) p.wt <- rep(1, n)
    if (!is.null(stats::weights(mod))) p.wt <- stats::weights(mod)
  }
  if (family != "gaussian") p.wt <- stats::weights(mod)

  bigp.wt <- rep(p.wt, B)

  if(length(names(mod$model)) > (d + 1))  W1 <- as.matrix(mod$model[, -c(1, (d + 2))])
  if(length(names(mod$model)) == (d + 1))  W1 <- as.matrix(mod$model[, -1])

  W1 <- as.matrix(mod$model[, -1])

  if(is.null(W)) W <- W1

  muPred <- rep(stats::predict(mod, type = "response"), B)

  beta.est <- stats::coef(mod)

  form.name <- stats::formula(mod)

  if (is.matrix(sigma.sq.u) == F) {
    w1 <- W[, 1]
    sigma.sq.u1 <- sigma.sq.u
    sigma.sq.e1 <- sigma.sq.e

    U1_j <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1, B)))
    X1_j <- rep(w1, B) - U1_j

    mu.e1 <- mean(X1_j)

    X <- cbind(X1_j)
    x1 <- X1_j

    dat_new <- cbind(bigY, x1)

    colnames(dat_new)[2] <- "x1"

    form.name <- stats::as.formula(stats::update(stats::formula(mod), bigY ~ s(x1) +. - s(w1)))

    if (d > 1) {
      smooth.err <- which(W1[1, ]%in%w1[1])
      non.err.names <- mod.terms[-smooth.err]
      W2 <- as.matrix(W1[, -smooth.err])

      for(jj in 1:(d - 1)) {
        dat_new <- cbind(dat_new, rep(W2[, jj], B))
        colnames(dat_new)[dim(dat_new)[2]] <- non.err.names[jj]
      }
    }
  }

  if (is.matrix(sigma.sq.u) == T) {
    q1 <- dim(W)[2]
    if (q1 == d) {
      col.nameX <- paste0('x', 1:q1)
      X <- W
    }

    if (d > q1) {
      col.nameX <- c()
      for(kk in 1:q1) {
        col.nameX1 <- paste0('x', kk)
        col.nameX <- c(col.nameX, col.nameX1)
      }
    }

    if (d < q1) print("Error: Number of error-contaminated covariates exceeds total number of covariates!")

    sigma.sq.u1 <- sigma.sq.u
    sigma.sq.e1 <- diag(sigma.sq.e)
    w1 <- W[, 1]
    U1_j <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1[1, 1], B)))
    X1_j <- rep(w1, B) - U1_j
    X <- cbind(X1_j)
    mu.e1 <- mean(X1_j)

    for(kk in 2:q1) {
      w1 <- W[, kk]
      U1_j <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1[kk, kk], B)))
      X1_j <- rep(w1, B) - U1_j

      X <- cbind(X, X1_j)
      mu.e1 <- c(mu.e1, mean(X1_j))
    }

    dat_new <- cbind(bigY, X)

    colnames(dat_new)[-1] <- col.nameX

    if (d > q1) {
      smooth.err <- which(W1[1, ]%in%W[1, ])

      non.err.names <- mod.terms[-smooth.err]

      W2 <- as.matrix(W1[, -smooth.err])

      for(jj in 1:(d - q1)) {
        dat_new <- cbind(dat_new, rep(W2[, jj], B))
        colnames(dat_new)[dim(dat_new)[2]] <- non.err.names[jj]
      }
    }

    if (q1 == 2) form.name <- stats::as.formula(stats::update(stats::formula(mod), bigY ~ s(x1) + s(x2) +. - s(w1) - s(w2)))
    if (q1 == 3) form.name <- stats::as.formula(stats::update(stats::formula(mod), bigY ~ s(x1) + s(x2) + s(x3) +. - s(w1) - s(w2) - s(w3)))
    if (q1 == 4) form.name <- stats::as.formula(stats::update(stats::formula(mod), bigY ~ s(x1) + s(x2) + s(x3) + s(x4) +. - s(w1) - s(w2) - s(w3) - s(w4)))
    if (q1 == 5) form.name <- stats::as.formula(stats::update(stats::formula(mod), bigY ~ s(x1) + s(x2) + s(x3) + s(x4) + s(x5) +. - s(w1) - s(w2) - s(w3) - s(w4) - s(w5)))
    if (q1 > 5) print("Sorry, refitME cannot handle more than 5 error-contaminated varaibles with GAM! But we're working on this...")
  }

  dat_new <- as.data.frame(dat_new)

  while(cond) {

    # MC and E-step.

    if (is.matrix(sigma.sq.u) == F) prX <- stats::dnorm(X1_j, mu.e1, sd = sqrt(sigma.sq.e1))
    if (is.matrix(sigma.sq.u) == T) prX <- mvtnorm::dmvnorm(X, mu.e1, sigma = sqrt(sigma.sq.e1))

    if (family == "gaussian") prY <- stats::dnorm(bigY, muPred, 1)
    if (family == "binomial") prY <- stats::dbinom(bigY, 1, muPred)
    if (family == "poisson") prY <- stats::dpois(bigY, muPred)
    if (family == "Gamma") prY <- stats::dgamma(bigY, rate = shape.est/muPred, shape = shape.est)
    if (family == "negbin") prY <- stats::dnbinom(bigY, size = theta.est, mu = muPred)

    prY[prY == 0] <- 1.e-6

    bigW <- matrix(prY*prX, n, B)
    sumW <- rep(apply(bigW, 1, sum, na.rm = T), B)
    weights1 <- as.vector(bigW)/sumW

    weights2 <- weights1
    weights1 <- bigp.wt*weights1

    dat_new$weights1 <- weights1

    # M-step (updates).

    if (family == "gaussian") {
      mod <- mgcv::gam(formula = form.name, weights = weights1, data = dat_new)
      sigma.sq.est <- (summary(mod)$dispersion)
    }
    if (family == "binomial") mod <- mgcv::gam(formula = form.name, weights = weights1, family = "binomial", gamma = 1.4, data = dat_new)
    if (family == "poisson") mod <- mgcv::gam(formula = form.name, weights = weights1, family = "poisson", data = dat_new)
    if (family == "Gamma") mod <- mgcv::gam(formula = form.name, weights = weights1, family = Gamma(link = "log"), data = dat_new)
    if (family == "negbin") mod <- mgcv::gam(form.name, family = nb(), weights = weights1, data = dat_new)

    beta.update <- stats::coef(mod)

    if (family == "negbin") theta.update <- mod$family$getTheta(TRUE)
    if (family == "Gamma") shape.update <- MASS::gamma.shape(mod)[1]$alpha

    muPred <- mgcv::predict.gam(mod, type = "response")

    if (is.matrix(sigma.sq.u) == F) {
      sigma.sq.e1.update <- wt.var(X[, 1], w = weights2)
      mu.e1.update <- stats::weighted.mean(X[, 1], w = weights2)
    }

    if (is.matrix(sigma.sq.u) == T) {
      sigma.sq.e1.update <- c()
      mu.e1.update <- c()

      for(kk in 1:q1) {
        sigma.sq.e1.update1 <- wt.var(X[, kk], w = weights2)
        sigma.sq.e1.update <- c(sigma.sq.e1.update, sigma.sq.e1.update1)

        mu.e1.update1 <- stats::weighted.mean(X[, kk], w = weights2)
        mu.e1.update <- c(mu.e1.update, mu.e1.update1)
      }
      sigma.sq.e1.update <- diag(sigma.sq.e1.update)
    }

    # Convergence monitoring.

    beta.norm <- sum((beta.est - beta.update)^2)

    if (family == "negbin") theta.norm <- sum((theta.est - theta.update)^2)
    if (family == "Gamma") shape.norm <- sum((shape.est - shape.update)^2)

    if (is.matrix(sigma.sq.u) == F) diff.sig_e <- abs(sigma.sq.e1.update - sigma.sq.e1)
    if (is.matrix(sigma.sq.u) == T) diff.sig_e <- sum(abs(diag(sigma.sq.e1.update) - diag(sigma.sq.e1)))

    diff.mu_e <- sum((mu.e1.update - mu.e1)^2)

    reps <- reps + 1 # Keeps track of number of iterations.

    if (family == "binomial" | family == "poisson" | family == "gaussian") {
      if ((diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon) | reps > 50) {
        cond <- FALSE
        print("convergence :-)")
        print(reps)
        break
      }
    }

    if (family == "negbin") {
      if ((diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon & theta.norm < epsilon) | reps > 50) {
        cond <- FALSE
        print("convergence :-)")
        print(reps)
        break
      }
    }

    if (family == "Gamma") {
      if ((diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon & shape.norm < epsilon) | reps > 50) {
        cond <- FALSE
        print("convergence :-)")
        print(reps)
        break
      }
    }

    # Update parameters.

    beta.est <- beta.update
    if (family == "negbin") theta.est <- theta.update
    if (family == "Gamma") shape.est <- shape.update
    sigma.sq.e1 <- sigma.sq.e1.update
    mu.e1 <- mu.e1.update
  }

  sumW <- apply(bigW, 1, sum, na.rm = T)
  weights1 <- bigW/sumW
  weights1[is.nan(weights1)] <- 0

  eff.samp.size <- 1/apply(weights1^2, 1, sum)
  eff.samp.size[is.infinite(eff.samp.size)] <- "NA"
  eff.samp.size <- as.numeric(eff.samp.size)

  mod$eff.samp.size <- eff.samp.size

  # Standard error calculations start here.

  beta.est <- stats::coef(mod)

  estfun_mat <- sandwich::estfun(mod)

  if (family == "gaussian") estfun_mat <- sandwich::estfun(mod)*B

  X <- stats::model.matrix(mod)
  K1 <- ncol(X)

  ind_mat <- matrix(1:(n*B), ncol = n, byrow = T)

  S_1 <- matrix(0, nrow = K1, ncol = K1)
  SS_1 <- matrix(0, nrow = K1, ncol = K1)

  for(ii in 1:n) {
    index_vec <- ind_mat[, ii]
    S_1 <- S_1 + (apply(estfun_mat[index_vec, ], 2, sum))%*%t(apply(estfun_mat[index_vec, ], 2, sum))
  }

  if (family == "gaussian") {
    u.bar <- solve(stats::vcov(mod)*B)
    beta.est.se2 <- sqrt(diag(solve(u.bar - SS_1/B^2 + S_1/B^2)))
    mod$Vp <- solve(u.bar - SS_1/B^2 + S_1/B^2)
  }

  if (family == "binomial" | family == "poisson" | family == "Gamma" | family == "negbin") {
    sand1 <- (sandwich::estfun(mod)/mod$prior)
    sand1[is.nan(sand1)] <- 1
    SS_1 <- t(sandwich::estfun(mod))%*%sand1
    u.bar <- solve(stats::vcov(mod))
    beta.est.se2 <- sqrt(diag(solve(u.bar - SS_1 + S_1)))
    mod$Vp <- solve(u.bar - SS_1 + S_1)
  }

  if (length(which(is.nan(beta.est.se2))) > 0) beta.est.se2 <- c(rep(NA, K1))

  mod$se <- beta.est.se2

  return(mod)
}

#' MCEMfit_CR
#'
#' Function for fitting VGAM capture-recapture (CR) model using the MCEM algorithm where covariates have measurement error.
#' @name MCEMfit_CR
#' @param mod : a vglm/vgam object (this is the naive fitted CR model). Make sure the first input predictor variables are the selected error-contaminated variable (i.e., the \code{W}'s).
#' @param sigma.sq.u : measurement error variance. A scalar if there is only one error-contaminated variable, otherwise this must stored as a covariance matrix.
#' @param sigma.sq.e : variance of the true covariate (X).
#' @param B : the number of Monte Carlo replication values (default is set to 50).
#' @param epsilon : a set convergence threshold (default is set to 0.00001).
#' @return \code{MCEMfit_CR} returns model coefficient and population size estimates with standard errors and the effective sample size.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. \pkg{refitME}: Measurement Error Modelling using Monte Carlo Expectation Maximization in \proglang{R}.
#' @import MASS VGAM
#' @importFrom VGAM s
#' @export
#' @seealso \code{\link{MCEMfit_glm}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown tutorial with examples.
MCEMfit_CR <- function(mod, sigma.sq.u, sigma.sq.e = 1, B = 100, epsilon = 0.00001) {

  options(warn = -1)

  reps <- 0
  cond <- TRUE

  tau <- mod@extra$tau[1]

  mod.terms <- attr(mod@terms$terms, "term.labels")

  d <- length(mod.terms)

  n <- mod@misc$n

  beta.est <- coef(mod)

  w1 <- mod@x[, 2]

  W1 <- model.matrix(mod)[, -1]

  if (d == 1 & mod@misc$function.name == "vglm") p1 <- 2

  if (d == 1 & mod@misc$function.name == "vgam") p1 <- "spline"

  if (d > 1 & mod@misc$function.name == "vglm") {
    if (sum((w1)^2 != (w1[, 2])) == n) p1 <- 1
    if (sum((w1)^2 == (W1[, 2])) == n) p1 <- 2
  }

  sigma.sq.u1 <- sigma.sq.u
  sigma.sq.e1 <- sigma.sq.e

  U1_j <- stats::rnorm(n*B, 0, sd = sqrt(sigma.sq.u1))
  X1_j <- rep(w1, B) - U1_j

  mu.e1 <- mean(X1_j)

  if (p1 == 2) {
    X <- cbind(rep(1, B*n), X1_j)
    muPred <- rep(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))), B)
  }

  if (p1 == 3) {
    X <- cbind(rep(1, B*n), X1_j, (X1_j)^2)
    muPred <- rep(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))), B)
  }

  if (p1 == "spline") {
    X <- cbind(rep(1, B*n), X1_j)
    muPred <- rep(1/(1 + exp(-VGAM::predict.vgam(mod, type = "link"))), B)
  }

  if (p1 == 2) colnames(X) <- c(names(coef(mod))[1], "x1")
  if (p1 == 3) colnames(X) <- c(names(coef(mod))[1], "x1", "x2")

  bigY <- rep(mod@y*tau, B)

  CR_dat <- data.frame(cbind(X1_j, bigY, tau - bigY))
  colnames(CR_dat) <- c("x1", "cap", "noncap")

  while(cond) {

    # MC and E-step.

    prX <- stats::dnorm(X1_j, mu.e1, sd = sqrt(sigma.sq.e1))
    prY <- VGAM::dposbinom(bigY, tau, muPred)

    # M-step (updates).

    bigW <- matrix(prY*prX, n, B)

    sumW <- rep(apply(bigW, 1, sum), B)

    weights1 <- as.vector(bigW)/sumW

    if (p1 == "spline") mod <- VGAM::vgam(cbind(cap, noncap) ~ s(x1, df = 2), VGAM::posbinomial(omit.constant = TRUE, parallel = TRUE ~ s(x1, df = 2)), data = CR_dat, trace = F, weights = as.vector(bigW)/sumW)
    if (p1 == 2) mod <- VGAM::vglm(cbind(cap, noncap) ~ x1, VGAM::posbinomial(omit.constant = TRUE, parallel = TRUE ~ x1), data = CR_dat, trace = F, weights = weights1)
    if (p1 == 3) mod <- VGAM::vglm(cbind(cap, noncap) ~ x1 + I(x1^2), VGAM::posbinomial(omit.constant = TRUE, parallel = TRUE ~ x1 + I(x1^2)), data = CR_dat, trace = F, weights = as.vector(bigW)/sumW)

    beta.update <- VGAM::coef(mod)
    if (p1 == 2 | p1 == 3) muPred <- 1/(1 + exp(-VGAM::predictvglm(mod, type = "link")))
    if (p1 == "spline") muPred <- 1/(1 + exp(-VGAM::predict.vgam(mod, type = "link")))

    sigma.sq.e1.update <- wt.var(X[, 2], w = as.vector(bigW)/sumW)
    mu.e1.update <- stats::weighted.mean(X[, 2], w = ((as.vector(bigW)/sumW)))

    # Convergence monitoring.

    beta.norm <- sum((beta.est - beta.update)^2)
    diff.sig_e <- abs(sigma.sq.e1.update - sigma.sq.e1)
    diff.mu_e <- sum((mu.e1.update - mu.e1)^2)

    reps <- reps + 1   # Keeps track of number of iterations.

    if ((diff.sig_e < epsilon & beta.norm < epsilon) | reps > 50) {
      cond <- FALSE
      print("convergence :-)")
      print(reps)
      break
    }

    # Update parameters.

    beta.est <- beta.update
    sigma.sq.e1 <- sigma.sq.e1.update
    mu.e1 <- mu.e1.update
  }

  beta.est.se1 <- sqrt(diag(VGAM::vcov(mod))) # Naive SE estimator.

  sumW <- apply(bigW, 1, sum, na.rm = T)
  weights1 <- bigW/sumW
  weights1[is.nan(weights1)] <- 0

  eff.samp.size <- 1/apply(weights1^2, 1, sum)
  eff.samp.size[is.infinite(eff.samp.size)] <- "NA"
  eff.samp.size <- as.numeric(eff.samp.size)

  if (p1 == 2 | p1 == 3) {
    pr.est <- 1/(1 + exp(-VGAM::predictvglm(mod, type = "link")))
    pi.est <- 1 - (1 - pr.est)^tau
  }

  if (p1 == "spline") {
    pr.est <- 1/(1 + exp(-VGAM::predict.vgam(mod, type = "link")))
    pi.est <- 1 - (1 - pr.est)^tau
  }

  w.est <- as.numeric(weights1)
  N.est <- sum(w.est/pi.est)

  # Standard error calculations start here.

  K1 <- length(beta.update)

  S_1 <- matrix(0, nrow = K1, ncol = K1)
  SS_1 <- matrix(0, nrow = K1, ncol = K1)
  II_1 <- matrix(0, nrow = K1, ncol = K1)

  pi.est_1 <- c()
  pr.est_1 <- c()

  iii1 <- 1

  ind_mat <- matrix(1:(n*B), ncol = n, byrow = T)

  for(iii in 1:n) {
    index_vec <- ind_mat[, iii]
    x <- X[index_vec, ]

    if (p1 == 2 | p1 == 3) {
      SS_1 <- SS_1 + t(x*(((bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))))[index_vec]/pi.est[index_vec]))*VGAM::weights(mod, type = "working")[index_vec]))%*%(x*((bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))))[index_vec]/pi.est[index_vec])))
      S_1 <- S_1 + (apply(x*(bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))))[index_vec]/pi.est[index_vec])*VGAM::weights(mod, type = "working")[index_vec], 2, sum))%*%t(apply(x*(bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))))[index_vec]/pi.est[index_vec])*VGAM::weights(mod, type = "working")[index_vec], 2, sum))
    }

    if (p1 == "spline") {
      SS_1 <- SS_1 + t(x*(((bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predict.vgam(mod, type = "link"))))[index_vec]/pi.est[index_vec]))*VGAM::weights(mod, type = "working")[index_vec]))%*%(x*((bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predict.vgam(mod, type = "link"))))[index_vec]/pi.est[index_vec])))
      S_1 <- S_1 + (apply(x*(bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predict.vgam(mod, type = "link"))))[index_vec]/pi.est[index_vec])*VGAM::weights(mod, type = "working")[index_vec], 2, sum))%*%t(apply(x*(bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predict.vgam(mod, type = "link"))))[index_vec]/pi.est[index_vec])*VGAM::weights(mod, type = "working")[index_vec], 2, sum))
    }
  }

  u.bar <- solve(VGAM::vcov(mod))
  beta.est.se <- sqrt(diag(solve(u.bar - SS_1 + S_1)))

  if (length(which(is.nan(beta.est.se))) > 0) beta.est.se2 <- c(rep(NA, K1))

  idG <- solve(u.bar - SS_1 + S_1)
  X.s <- VGAM::model.matrix(mod)

  dpi <- w.est*tau*pr.est*(1 - pi.est)/pi.est/pi.est
  dNhat <- apply(X.s*c(dpi), 2, sum)
  N.est.se <- sqrt(sum(w.est*(1 - pi.est)/pi.est/pi.est) + t(dNhat)%*%idG%*%dNhat)

  values <- list(beta = beta.est, beta.se = beta.est.se, N.est = N.est, N.est.se = N.est.se, mod = mod, eff.samp.size = eff.samp.size)

  return(values)
}

#' refitME
#'
#' Function that extracts the fitted (naive) model object and wraps the MCEM algorithm to correct for measurement error/error-in-variables (currently available for \code{lm()}, \code{glm()} and \code{gam()}, excludes \code{lme()}, \code{nlme()} and \code{polr()} models).
#' @name refitME
#' @param mod : a glm/gam object (this is the naive fitted model). Make sure the first input predictor variables are the selected error-contaminated variable (i.e., the \code{W}'s).
#' @param sigma.sq.u : measurement error variance. A scalar if there is only one error-contaminated variable, otherwise must be a covariance matrix
#' @param W : a matrix of error-contaminated covariates (if not specified, the default assumes all covariates in the naive fitted model are error-contaminated).
#' @param B : the number of Monte Carlo replication values (default is set 50).
#' @param epsilon : convergence threshold (default is set to 0.00001).
#' @return \code{refitME} returns the naive fitted model object where coefficient estimates, the covariance matrix, fitted values and residuals have been replaced with the final MCEM model fit. Standard errors and effective sample size (which diagnose how closely the proposal distribution matches the posterior, see equation (2) of Stoklosa, Hwang and Warton) have been additional included.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. \pkg{refitME}: Measurement Error Modelling using Monte Carlo Expectation Maximization in \proglang{R}.
#' @import mvtnorm MASS mgcv sandwich VGAM expm
#' @importFrom stats Gamma
#' @export
#' @seealso \code{\link{MCEMfit_glm}} and \code{\link{MCEMfit_gam}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown tutorial with examples.
refitME <- function(mod, sigma.sq.u, W = NULL, B = 50, epsilon = 0.00001) {
  if (is.matrix(sigma.sq.u) == F) {
    print("One specified error-contaminated covariate.")

    if (attr(mod, "class")[1] == "lm" | attr(mod, "class")[1] == "glm" | attr(mod, "class")[1] == "negbin" | attr(mod, "class")[1] == "gam") sigma.sq.e <- stats::var(mod$model[, 2]) - sigma.sq.u
    if (attr(mod, "class")[1] == "vglm" | attr(mod, "class")[1] == "vgam") sigma.sq.e <- stats::var(mod@x[, 2]) - sigma.sq.u
  }

  if (is.matrix(sigma.sq.u) == T) {
    print("Multiple specified error-contaminated covariates.")
    q1 <- dim(W)[2]
    sigma.sq.e <- c()

    for(kk in 1:q1) {
      sigma.sq.e1 <- stats::var(W[, kk]) - sigma.sq.u[kk, kk]
      sigma.sq.e <- c(sigma.sq.e, sigma.sq.e1)
    }
  }

  ob.type <- attr(mod, "class")[1]

  if ((ob.type == "lm" | ob.type == "glm" | ob.type == "negbin")) {
    if (ob.type == "lm") family <- "gaussian"
    if (ob.type == "glm") family <- mod$family$family
    if (ob.type == "negbin") family <- "negbin"

    return(MCEMfit_glm(mod, family, sigma.sq.u, W, sigma.sq.e, B, epsilon))
  }

  if (ob.type == "gam") {
    family <- mod$family$family
    if (strsplit(family, NULL)[[1]][1] == "N") family <- "negbin"

    return(MCEMfit_gam(mod, family, sigma.sq.u, W, sigma.sq.e, B, epsilon))
  }

  if ((ob.type == "vglm" | ob.type == "vgam")) return(MCEMfit_CR(mod, sigma.sq.u, sigma.sq.e, B, epsilon))

}

#' wt.var
#'
#' Function that calculates a weighted variance.
#' @name wt.var
#' @param x : a vector of numerical data.
#' @param w : a vector of equal length to \code{x} representing the weights.
#' @return \code{wt.var} returns a single value from analysis requested.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
#' # Define simple data
#' x = 1:25 # Set of numbers.
#' wt = runif(25) # Some arbitrary weights.
#'
#' # Display variances (unweighted and then weighted).
#' var(x)
#' wt.var(x, wt)
#' @export
#' @source See \url{https://rdrr.io/cran/SDMTools/src/R/wt.mean.R}
wt.var <- function(x, w) {
  s <- which(is.finite(x + w))
  wt <- w[s]
  x <- x[s] # Remove NA info.
  xbar <- stats::weighted.mean(x, w = wt) # Get the weighted mean.

  return(sum(wt*(x - xbar)^2)*(sum(wt)/(sum(wt)^2 - sum(wt^2)))) # Return the variance.
}
