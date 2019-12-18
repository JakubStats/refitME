#-----------------------------------------------------------------------
# refitME_package.R
#
# Source file for the "refitME" R-package.
#-----------------------------------------------------------------------

# Required R-packages.

library(mvtnorm)
library(MASS)
library(SDMTools)
library(mgcv)
library(Matrix)
suppressMessages(library(sandwich))

#' @title The Framingham heart study data set.
#' @description Data set consisting of records of male patients with coronary heart disease collected from the Framingham heart study. The \code{Framinghamdata} data consists of binary responses and four predictor variables collected on `n=1615` patients.
#' @format A data set that contains: 5 columns with 1615 observations. The columns are defined as follows:
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
#' Function for wrapping the MCEM algorithm on GLMs where covariates are subject to measurement error/error-in-varaibles.
#' @name MCEMfit_glm
#' @param mod : a glm object (this is the naive fitted model). Make sure the first input predictor variables are the selected error-contaminated varaible (i.e., the W's).
#' @param family : a specified family/distribution.
#' @param sigma.sq.u : measurement error variance. A scaler if there is only one error-contaminated variable, otherwise this must stored as a covaraince matrix.
#' @param W a matrix of error-contaminated covariates.
#' @param sigma.sq.e : variance of the true covariate (X).
#' @param B : the number of Monte Carlo replication values (default is set to 50).
#' @param epsilon : a set convergence threshold (default is set to 0.00001).
#' @param theta.est : an initial value for the dispersion parameter (this is required for fitting negative binomial models).
#' @param shape.est : an initial value for the shape parameter (this is required for fitting gamma models).
#' @return \code{MCEMfit_glm} returns model coef estimates with standard errors and the effective sample size.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. (2019). A general algorithm for error-in-variables modelling using Monte Carlo expectation maximization.
#' @import mvtnorm MASS SDMTools mgcv sandwich
#' @importFrom stats Gamma
#' @export
#' @seealso \code{\link{MCEMfit_gam}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown tutorial with examples.
MCEMfit_glm <- function(mod, family, sigma.sq.u, W, sigma.sq.e = 1, B = 50, epsilon = 0.00001, theta.est = 1, shape.est = 1) {

  options(warn = -1)

  reps <- 0
  cond <- TRUE

  mod.terms <- attr(mod$terms, "term.labels")

  d <- length(mod.terms)

  Y <- mod$model[, 1]
  bigY <- rep(Y, B)

  n <- length(Y)

  W1 <- mod$model[, -1]

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

    if (q1 == d) {
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

    # M-step (updates).

    bigW <- matrix(prY*prX, n, B)
    sumW <- rep(apply(bigW, 1, sum, na.rm = T), B)

    weights1 <- as.vector(bigW)/sumW
    weights1[is.nan(weights1)] <- 0

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
    if (family == "Gamma") {
      shape.update <- summary(mod)[14]$dispersion
    }

    if (is.matrix(sigma.sq.u) == F) {
      sigma.sq.e1.update <- SDMTools::wt.var(X[, 2], w = weights1)
      mu.e1.update <- stats::weighted.mean(X[, 2], w = weights1)
    }

    if (is.matrix(sigma.sq.u) == T) {
      sigma.sq.e1.update <- c()
      mu.e1.update <- c()

      for(kk in 1:q1) {
        sigma.sq.e1.update1 <- SDMTools::wt.var(XA[, kk], w = weights1)
        sigma.sq.e1.update <- c(sigma.sq.e1.update, sigma.sq.e1.update1)

        mu.e1.update1 <- stats::weighted.mean(XA[, kk], w = weights1)
        mu.e1.update <- c(mu.e1.update1, mu.e1.update)
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

  sumW <- apply(bigW, 1, sum, na.rm = T)
  weights1 <- bigW/sumW
  weights1[is.nan(weights1)] <- 0

  eff.samp.size <- 1/apply(weights1^2, 1, sum)
  eff.samp.size[is.infinite(eff.samp.size)] <- "NA"
  eff.samp.size <- as.numeric(eff.samp.size)

  # Standard error calculations start here.

  beta.est.se1 <- sqrt(diag(stats::vcov(mod)))  # Naive SE estimator.

  if (family == "gaussian") beta.est.se1 <- sqrt(diag(stats::vcov(mod)*B))

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
  }

  if (family!="gaussian") {
    sand1 <- (sandwich::estfun(mod)/mod$prior)
    sand1[is.nan(sand1)] <- 1
    SS_1 <- t(sandwich::estfun(mod))%*%sand1
    u.bar <- solve(stats::vcov(mod))
    beta.est.se2 <- sqrt(diag(solve(u.bar - SS_1 + S_1)))
  }

  if (length(which(is.nan(beta.est.se2))) > 0) beta.est.se2 <- c(rep(NA, K1))

  values <- list(beta = beta.est, beta.se1 = beta.est.se1, beta.se2 = beta.est.se2, mod = mod, eff.samp.size = eff.samp.size)

  return(values)
}

#' MCEMfit_gam
#'
#' Function for wrapping the MCEM algorithm on GAMs where covariates are subject to measurement error/error-in-varaibles.
#' @name MCEMfit_gam
#' @param mod : a gam object (this is the naive fitted model). Make sure the first input predictor variables are the selected error-contaminated varaible (i.e., the W's).
#' @param family : a specified family/distribution.
#' @param sigma.sq.u : measurement error variance. A scaler if there is only one error-contaminated variable, otherwise this must stored as a covaraince matrix.
#' @param W a matrix of error-contaminated covariates.
#' @param sigma.sq.e : variance of the true covariate (X).
#' @param B : the number of Monte Carlo replication values (default is set to 50).
#' @param epsilon : convergence threshold (default is set to 0.00001).
#' @return \code{refitME} returns model coef estimates with standard errors and the effective sample size.
#' @param theta.est : an initial value for the dispersion parameter (this is required for fitting negative binomial models).
#' @param shape.est : an initial value for the shape parameter (this is required for fitting gamma models).
#' @return \code{MCEMfit_glm} returns model coef estimates with standard errors.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. (2019). A general algorithm for error-in-variables modelling using Monte Carlo expectation maximization.
#' @import mvtnorm MASS SDMTools mgcv sandwich
#' @importFrom stats Gamma
#' @export
#' @seealso \code{\link{MCEMfit_glm}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown tutorial with examples.
MCEMfit_gam <- function(mod, family, sigma.sq.u, W, sigma.sq.e = 1, B = 50, epsilon = 0.00001, theta.est = 1, shape.est = 10) {

  options(warn = -1)

  reps <- 0
  cond <- TRUE

  mod.terms <- attr(mod$terms, "term.labels")

  d <- length(mod.terms)

  Y <- mod$model[, 1]
  bigY <- rep(Y, B)

  n <- length(Y)

  W1 <- as.matrix(mod$model[, -1])

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

    # M-step (updates).

    bigW <- matrix(prY*prX, n, B)
    sumW <- rep(apply(bigW, 1, sum), B)

    weights1 <- as.vector(bigW)/sumW
    weights1[is.nan(weights1)] <- 0

    dat_new$weights1 <- weights1

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
    muPred <- stats::predict(mod, type = "response")

    if (is.matrix(sigma.sq.u) == F) {
      sigma.sq.e1.update <- SDMTools::wt.var(X[, 1], w = weights1)
      mu.e1.update <- stats::weighted.mean(X[, 1], w = weights1)
    }

    if (is.matrix(sigma.sq.u) == T) {
      sigma.sq.e1.update <- c()
      mu.e1.update <- c()

      for(kk in 1:q1) {
        sigma.sq.e1.update1 <- SDMTools::wt.var(X[, kk], w = weights1)
        sigma.sq.e1.update <- c(sigma.sq.e1.update, sigma.sq.e1.update1)

        mu.e1.update1 <- stats::weighted.mean(X[, kk], w = weights1)
        mu.e1.update <- c(mu.e1.update1, mu.e1.update)
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

  # Standard error calculations start here.

  beta.est <- stats::coef(mod)
  beta.est.se1 <- sqrt(diag(stats::vcov(mod)))   # Naive SE estimator.

  if (family == "gaussian") beta.est.se1 <- sqrt(diag(stats::vcov(mod)*B))

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
  }

  if (family == "binomial" | family == "poisson" | family == "Gamma" | family == "negbin") {
    sand1 <- (sandwich::estfun(mod)/mod$prior)
    sand1[is.nan(sand1)] <- 1
    SS_1 <- t(sandwich::estfun(mod))%*%sand1
    u.bar <- solve(stats::vcov(mod))
    beta.est.se2 <- sqrt(diag(solve(u.bar - SS_1 + S_1)))
  }

  if (length(which(is.nan(beta.est.se2))) > 0) beta.est.se2 <- c(rep(NA, K1))

  values <- list(beta = beta.est, beta.se1 = beta.est.se1, beta.se2 = beta.est.se2, mod = mod, eff.samp.size = eff.samp.size)

  return(values)
}

#' refitME
#'
#' Function that extracts the fitted (naive) model object and wraps the MCEM algorithm to correct for measurement error/error-in-varaibles (currently available for lm, glm and gam, excludes quassi models).
#' @name refitME
#' @param mod : a glm/gam object (this is the naive fitted model). Make sure the first input predictor variables are the selected error-contaminated varaible (i.e., the W's).
#' @param sigma.sq.u : measurement error variance. A scaler if there is only one error-contaminated variable, otherwise must be a covaraince matrix
#' @param W a matrix of error-contaminated covariates.
#' @param B : the number of Monte Carlo replication values (default is set 50).
#' @param epsilon : convergence threshold (default is set to 0.00001).
#' @return \code{refitME} returns model coef estimates with standard errors and the effective sample size.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. (2019). A general algorithm for error-in-variables modelling using Monte Carlo expectation maximization.
#' @import mvtnorm MASS SDMTools mgcv sandwich
#' @importFrom stats Gamma
#' @export
#' @seealso \code{\link{MCEMfit_glm}} and \code{\link{MCEMfit_gam}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown tutorial with examples.
refitME <- function(mod, sigma.sq.u, W, B = 50, epsilon = 0.00001) {

  if (is.matrix(sigma.sq.u) == F) {
    print("One specified error-contaminated covariate.")
    sigma.sq.e <- stats::var(mod$model[, 2]) - sigma.sq.u
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

  if (ob.type == "lm" | ob.type == "glm" | ob.type == "negbin") {
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
}
