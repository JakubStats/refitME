#-----------------------------------------------------------------------
# refitME_package.R
#
# Source file for the "refitME" R-package.
#-----------------------------------------------------------------------

#' @title The Framingham heart study data set
#' @description Data set consisting of records of male patients with coronary heart disease collected from the Framingham heart study. The \code{Framinghamdata} data consists of binary responses and four predictor variables collected on `n = 1615` patients.
#' @format A data set that contains: 5 columns with 1,615 observations. The columns are defined as follows:
#' \describe{
#' \item{\code{Y}}{Response indicator (binary variable) of first evidence of CHD status of patient.}
#' \item{\code{z1}}{Serum cholesterol level of patient.}
#' \item{\code{z2}}{Age of patient.}
#' \item{\code{z3}}{Smoking indicator - whether the patient smokes.}
#' \item{\code{w1}}{Systolic blood pressure (SBP) of patient - this is the error contaminated variable, calculated from mean scores. The measurement error is 0.00630, see pp. 112 of Carroll \emph{et al.} (2006).}
#' }
#' @source See Carroll \emph{et al.} (2006) for full details on the data and study. Also, see \url{https://github.com/JakubStats/refitME} for an RMarkdown vignette of an example that uses the data.
#' @references Carroll, R. J., Ruppert, D., Stefanski, L. A., and Crainiceanu, C. M. (2006). \emph{Measurement Error in Nonlinear Models: A Modern Perspective.} 2nd Ed. London: Chapman & Hall/CRC.
#' @examples # Load the data.
#'
#' data(Framinghamdata)
"Framinghamdata"

#' @title The Corymbia eximia presence-only data set
#' @description Data set consisting of presence-only records for the plant species \emph{Corymbia eximia}, site coordinates 5 covariates for each site.
#' @format A data set that contains: 8 columns with 86,316 observations (or sites). The columns are defined as follows:
#' \describe{
#' \item{\code{X}}{Longitude coordinate.}
#' \item{\code{Y}}{Latitude coordinate.}
#' \item{\code{FC}}{Recorded number of fire counts for each site.}
#' \item{\code{MNT}}{Recorded minimum temperatures for each site.}
#' \item{\code{MXT}}{Recorded maximum temperature for each site.}
#' \item{\code{Rain}}{Recorded rainfall for each site.}
#' \item{\code{D.Main}}{Recorded distance from nearest major road.}
#' \item{\code{Y.obs}}{Presences for the plant species \emph{Corymbia eximia} for each site.}
#' }
#' @source See Renner and Warton (2013) for full details on the data and study.
#' @references Renner, I. W. and Warton, D. I. (2013). Equivalence of MAXENT and Poisson point process models for species distribution modeling in ecology. \emph{Biometrics}, \strong{69}, 274–281.
#' @examples # Load the data.
#'
#' data(Corymbiaeximiadata)
"Corymbiaeximiadata"

#' @title The yellow-bellied Prinia \emph{Prinia flaviventris} capture-recapture data
#' @description Data set consisting of capture-recapture histories 164 uniquely captured birds across 17 weekly capture occasions. Bird wing lengths were also measured in the study.
#' @format A data set that contains: 3 columns with 164 observations. The columns are defined as follows:
#' \describe{
#' \item{\code{w1}}{Bird wing lengths.}
#' \item{\code{cap}}{Number of times the individual was captured.}
#' \item{\code{noncap}}{Number of times the individual was not captured.}
#' }
#' @source See Hwang, Huang and Wang (2007) for full details on the data and study.
#' @references Hwang, W. H., Huang, S. Y. H., and Wang, C. (2007). Effects of measurement error and conditional score estimation in capture--recapture models. \emph{Statistica Sinica}, \strong{17}, 301-316.
#' @examples # Load the data.
#'
#' data(Priniadata)
"Priniadata"

#' A wrapper function for correcting measurement error in predictor variables via the MCEM algorithm
#'
#' Function that extracts the fitted (naive) model object and wraps the MCEM algorithm to correct for measurement error/error-in-variables in predictors.
#' @name refitME
#' @param mod : any (S3 class) fitted object that responds to the generic functions \code{family}, \code{model.frame}, \code{update} and \code{predict}, and accepts weighted observations via \code{weights}. The \code{mod} argument specifies the naive fitted model. Make sure the first \eqn{p} input predictor variables in the naive model are the selected error-contaminated predictors variables. Also, the \code{mod} argument allows \code{vlgm/vgam} (S4 class) model objects when using the \code{posbinomial} family -- this is a specific function developed for fitting closed population capture--recapture models, see \code{\link{MCEMfit_CR}}.
#' @param sigma.sq.u : measurement error (ME) variance. A scalar if there is only one error-contaminated predictor variable, otherwise this must be stored as a vector (of known ME variances) or a matrix if the ME covariance matrix is known.
#' @param B : the number of Monte Carlo replication values (default is set 50).
#' @param epsilon : convergence threshold (default is set to 0.00001).
#' @param silent : if \code{TRUE}, the convergence message (which tells the user if the model has converged and reports the number of iterations required) is suppressed (default is set to \code{FALSE}).
#' @param ... : further arguments passed through to the function that was used to fit \code{mod}, that will be used in refitting. These need only be specified if making changes to the arguments as compared to the original call that produced \code{mod}.
#' @return \code{refitME} returns the naive fitted model object where coefficient estimates, the covariance matrix, fitted values, the log-likelihood, and residuals have been replaced with the final MCEM model fit. Standard errors are included and returned, if \code{mod} is a class of object accepted by the \pkg{sandwich} package (such as \code{glm}, \code{gam}, \code{survreg} and many more). The effective sample size (which diagnose how closely the proposal distribution matches the posterior, see equation (2) of Stoklosa, Hwang and Warton) have also been included as outputs.
#' @author Jakub Stoklosa, Wen-Han Hwang and David I. Warton.
#' @references Carroll, R. J., Ruppert, D., Stefanski, L. A., and Crainiceanu, C. M. (2006). \emph{Measurement Error in Nonlinear Models: A Modern Perspective.} 2nd Ed. London: Chapman & Hall/CRC.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. \pkg{refitME}: Measurement Error Modelling using Monte Carlo Expectation Maximization in \proglang{R}.
#' @export
#' @seealso \code{\link{MCEMfit_glm}}, \code{\link{MCEMfit_gam}} and \code{\link{MCEMfit_gen}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown vignette with examples.
#' @examples # A GLM example I - binary response data.
#'
#' library(refitME)
#'
#' data(Framinghamdata)
#'
#' glm_naiv <- glm(Y ~ w1 + z1 + z2 + z3, x = TRUE, family = binomial, data = Framinghamdata)
#'
#' # The error-contaminated predictor variable in this example is systolic blood pressure (w1).
#'
#' sigma.sq.u <- 0.01259/2 # ME variance, as obtained from Carroll et al. (2006) monograph.
#'
#' B <- 50  # The number of Monte Carlo replication values.
#'
#' glm_MCEM <- refitME(glm_naiv, sigma.sq.u, B)
#'
refitME <- function(mod, sigma.sq.u, B = 50, epsilon = 0.00001, silent = FALSE, ...) {
  if (!isS4(mod)) {
    if (stats::formula(mod$model)[-2] == ~1) stop("The fitted naive model is an intercept-only model. Please specify/include the error-contaminated predictor/covariate in your model fit.", call. = TRUE)

    if (stats::formula(mod$model)[-2] != ~1) {
      if (is.matrix(sigma.sq.u) == FALSE) {

        n <- dim(stats::model.frame(mod))[1]

        if (length(sigma.sq.u) == 1) {
          W1 <- as.data.frame(stats::model.frame(mod)[, -1])[1]

          if (dim(as.matrix(W1))[2] != 1) {
            if (class(stats::model.frame(mod)[, -1][, 1])[1] == "poly") W1 <- W1[, 1][, 1]
          }

          message("One specified error-contaminated predictor/covariate.")
        }

        if (length(sigma.sq.u) > 1) {
          message("Multiple specified error-contaminated predictors/covariates.")

          sigma.sq.u <- diag(sigma.sq.u)

          q1 <- dim(sigma.sq.u)[2]
        }
      }

      if (is.matrix(sigma.sq.u) == TRUE) {
        message("Multiple specified error-contaminated predictors/covariates.")

        q1 <- dim(sigma.sq.u)[2]
      }
    }

    ob.type <- attr(mod, "class")[1]

    if (ob.type == "lm" | ob.type == "glm" | ob.type == "negbin") {
      if (ob.type == "lm") family <- "gaussian"
      if (ob.type == "glm") family <- mod$family$family
      if (ob.type == "negbin") family <- "negbin"

      return(MCEMfit_glm(mod, family, sigma.sq.u, B, epsilon, silent, ...))
    }

    if (ob.type == "gam") {
      family <- mod$family$family
      if (strsplit(family, NULL)[[1]][1] == "N") family <- "negbin"

      return(MCEMfit_gam(mod, family, sigma.sq.u, B, epsilon, silent, ...))
    }

    if (ob.type != "lm" | ob.type != "glm" | ob.type != "negbin" | ob.type != "gam") {
      family <- mod$family$family
      if (strsplit(family, NULL)[[1]][1] == "N") family <- "negbin"

      return(MCEMfit_gen(mod, family, sigma.sq.u, B, epsilon, silent, ...))
    }
  }

  if (isS4(mod)) return(MCEMfit_CR(mod, sigma.sq.u, B, epsilon, silent))
}

#' Function for wrapping the MCEM algorithm on \code{lm} or \code{glm} objects
#'
#' Function for wrapping the MCEM algorithm on GLMs where predictors are subject to measurement error/error-in-variables.
#' @name MCEMfit_glm
#' @param mod : a \code{lm/glm} object (this is the naive fitted model). Make sure the first \eqn{p} input predictor variables entered in the naive model are the specified error-contaminated variables. These \eqn{p} predictors also need the measurement error variance to be specified in \code{sigma.sq.u}, see below.
#' @param family : a specified family/distribution.
#' @param sigma.sq.u : measurement error (ME) variance. A scalar if there is only one error-contaminated predictor variable, otherwise this must be stored as a vector (of ME variances) or a matrix if the ME covariance matrix is known.
#' @param B : the number of Monte Carlo replication values (default is set to 50).
#' @param epsilon : a set convergence threshold (default is set to 0.00001).
#' @param silent : if \code{TRUE}, the convergence message (which tells the user if the model has converged and reports the number of iterations required) is suppressed (default is set to \code{FALSE}).
#' @param ... : further arguments passed to \code{lm} or \code{glm}.
#' @return \code{MCEMfit_glm} returns the naive fitted model object where coefficient estimates, the covariance matrix, fitted values, the log-likelihood, and residuals have been replaced with the final MCEM model fit. Standard errors and the effective sample size (which diagnose how closely the proposal distribution matches the posterior, see equation (2) of Stoklosa, Hwang and Warton) have also been included as outputs.
#' @author Jakub Stoklosa, Wen-Han Hwang and David I. Warton.
#' @references Carroll, R. J., Ruppert, D., Stefanski, L. A., and Crainiceanu, C. M. (2006). \emph{Measurement Error in Nonlinear Models: A Modern Perspective.} 2nd Ed. London: Chapman & Hall/CRC.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. \pkg{refitME}: Measurement Error Modelling using Monte Carlo Expectation Maximization in \proglang{R}.
#' @importFrom stats Gamma
#' @importFrom mvtnorm dmvnorm
#' @importFrom MASS glm.nb
#' @importFrom sandwich estfun
#' @importFrom expm sqrtm
#' @export
#' @seealso \code{\link{MCEMfit_gam}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown vignette with examples.
#' @examples # A GLM example I - binary response data.
#'
#' library(refitME)
#'
#' data(Framinghamdata)
#'
#' glm_naiv <- glm(Y ~ w1 + z1 + z2 + z3, x = TRUE, family = binomial, data = Framinghamdata)
#'
#' # The error-contaminated predictor in this example is systolic blood pressure (w1).
#'
#' sigma.sq.u <- 0.006295  # ME variance, as obtained from Carroll et al. (2006) monograph.
#'
#' B <- 50 # The number of Monte Carlo replication values.
#'
#' glm_MCEM <- refitME(glm_naiv, sigma.sq.u, B)
#'
MCEMfit_glm <- function(mod, family, sigma.sq.u, B = 50, epsilon = 0.00001, silent = FALSE, ...) {

  mod_n <- mod

  if (family == "gaussian") {
    if (length(class(mod)) == 2) mod_n <- stats::lm(mod$formula)
    if (length(class(mod)) == 1) mod_n <- mod
    mod$data <- eval(stats::getCall(mod)$data, environment(stats::formula(mod)))
  }

  if (family == "negbin") {
    mod$data <- eval(stats::getCall(mod)$data, environment(stats::formula(mod)))
    theta.est <- mod$theta
  }

  if (family == "Gamma") shape.est <- summary(mod)[14]$dispersion

  reps <- 0
  cond <- TRUE

  Y <- stats::model.frame(mod)[, 1]
  bigY <- rep(Y, B)
  n <- length(Y)

  form.name <- stats::formula(mod)

  if (family == "gaussian") {
    if (is.null(stats::weights(mod))) p.wt <- rep(1, n)
    if (!is.null(stats::weights(mod))) p.wt <- stats::weights(mod)
  }
  if (family != "gaussian") p.wt <- stats::weights(mod)

  bigp.wt <- rep(p.wt, B)

  W <- as.data.frame(stats::model.frame(mod)[, -1])
  names(W) <- names(stats::model.frame(mod))[-1]

  d <- ncol(W)

  muPred <- rep(stats::predict(mod, type = "response"), B)
  beta.est <- stats::coef(mod)

  if (is.matrix(sigma.sq.u) == F) {
    w <- W[, 1]

    names.w <- names(W)[1]

    if (dim(as.matrix(W))[2] != 1) {
      if (class(stats::model.frame(mod)[, -1])[1] == "poly" | class(stats::model.frame(mod)[, -1][, 1])[1] == "poly") {
        names.w <- names(mod$data)[2]
        names(mod$data)[2] <- names.w
      }
    }

    if (class(w)[1] == "poly") w <- w[, 1]

    sigma.sq.u1 <- sigma.sq.u
    sigma.sq.e1 <- stats::var(w) - sigma.sq.u1

    U <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1, B)))
    X <- rep(w, B) - U

    mu.e1 <- mean(X)

    new.dat.temp <- as.data.frame(mod$data[, -c(1, which(names(mod$data) %in% names.w))])

    new.dat.temp0 <- which(!(names(mod$data) %in% names.w))[-1]

    if (length(new.dat.temp0) != 0) {
      new.dat <- as.data.frame(cbind(bigY, matrix(NA, nrow = B*n, ncol = (ncol(mod$data) - 1))))
      new.dat[, which(names.w == names(mod$data))] <- X

      for (j in 1:length(new.dat.temp0)) {
        W2 <- new.dat.temp[, j]

        if (class(W2)[1] == "poly") W2 <- W2[, 1]

        new.dat[, new.dat.temp0[j]] <- rep(W2, B)
        names.w3 <- names(mod$data)[new.dat.temp0[j]]
        names(new.dat)[new.dat.temp0[j]] <- names.w3
      }
    }

    if (length(new.dat.temp0) == 0) new.dat <- as.data.frame(cbind(bigY, X))

    names(new.dat)[c(1, which(names(mod$data) %in% names.w))] <- names(mod$data)[c(1, which(names(mod$data) %in% names.w))]
  }

  if (is.matrix(sigma.sq.u) == T) {
    q1 <- dim(sigma.sq.u)[2]

    names.w <- names(W)[1]

    w <- W[, 1]

    W1 <- W

    poly.trig <- FALSE

    if (class(stats::model.frame(mod)[, -1][, 1])[1] == "poly") {
      names.w <- names(mod$data)[2]
      W1 <- W[, 1]

      d <- ncol(W)
      poly.trig <- TRUE
    }

    if (class(w)[1] == "poly") {
      w <- w[, 1]
      poly.trig <- TRUE
    }

    d1 <- 0

    if (ncol(as.matrix(W1)) > 1) {
      if ((sum(W1[, 2] == w^2) == n)) {    # Check if quadratic.
        if (class(stats::model.frame(mod)[, -1][, 1])[1] != "poly") d <- d - 1
        d1 <- 1
      }
    }

    if (ncol(as.matrix(W1)) > 2) {    # Check if cubic.
      if (sum(W1[, 3] == w^3) == n) {
        if (class(stats::model.frame(mod)[, -1][, 1])[1] != "poly") d <- d - 1
        d1 <- 2
      }
    }
    if (ncol(as.matrix(W1)) > 3) {     # Check if quartic.
      if (sum(W1[, 4] == w^4) == n) {
        if (class(stats::model.frame(mod)[, -1][, 1])[1] != "poly") d <- d - 1
        d1 <- 3
      }
    }

    if (d < q1) stop("Number of error-contaminated covariates exceeds total number of covariates!")

    sigma.sq.u1 <- sigma.sq.u
    sigma.sq.e2 <- stats::var(w) - sigma.sq.u[1, 1]

    U <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1[1, 1], B)))
    X <- cbind(rep(w, B) - U)

    mu.e1 <- mean(X)

    for (kk in 2:q1) {

      if (poly.trig  == TRUE) {
        W2 <- stats::model.frame(mod)[, -c(1:kk)]

        if (ncol(as.matrix(W2)) == 1) {
          W1 <- w <- W[, ncol(W)]
          names.w1 <- names(mod$data)[length(names(mod$data))]

          poly.trig0 <- FALSE
        }

        if (ncol(as.matrix(W2)) > 1) {
          if (class(W2)[1] == "poly" | class(W2[, 1])[1] == "poly") {
              W1 <- W[, kk]
              w <- W1[, 1]

              poly.trig0 <- TRUE
          }

          if (class(W2)[1] != "poly" & class(W2[, 1])[1] != "poly") {
            W1 <- as.matrix(W[, -c(1:(kk + d1 - 1))])
            w <- as.matrix(W[, kk + d1])

            poly.trig0 <- FALSE
          }

          names.w1 <- names(mod$data)[kk + 1]
        }

        names.w <- c(names.w, names.w1)
      }

      if (poly.trig  == FALSE) {
        W2 <- stats::model.frame(mod)[, -c(1:(kk + d1))]

        if (ncol(as.matrix(W2)) == 1) {
          W1 <- w1 <- W[, ncol(W)]
          names.w1 <- colnames(W)[ncol(W)]

          poly.trig0 <- FALSE
        }

        if (ncol(as.matrix(W2)) > 1) {
          if (class(W2)[1] == "poly" | class(W2[, 1])[1] == "poly") {
            names.w1 <- names(mod$data)[kk + 1]
            W1 <- W[, kk + d1]
            w <- W1[, 1]

            poly.trig0 <- TRUE
          }

          if (class(W2)[1] != "poly" & class(W2[, 1])[1] != "poly") {
            W1 <- as.matrix(W2)
            w <- W2[, 1]
            names.w1 <- colnames(W1)[1]

            poly.trig0 <- FALSE
          }
        }

        names.w <- c(names.w, names.w1)
      }

      sigma.sq.e <- stats::var(w) - sigma.sq.u[kk, kk]
      sigma.sq.e2 <- c(sigma.sq.e2, sigma.sq.e)

      if (ncol(as.matrix(W1)) > 1) {
        if (sum(W1[, 2] == w^2) == n) d1 <- d1 + 1  # Check if quadratic.
      }

      if (ncol(as.matrix(W1)) > 2) {   # Check if cubic.
        if (sum(W1[, 3] == w^3) == n) d1 <- d1 + 2
      }

      if (ncol(as.matrix(W1)) > 3) {    # Check if quartic.
        if (sum(W1[, 4] == w^4) == n) d1 <- d1 + 3
      }

      poly.trig <- poly.trig0

      U <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1[kk, kk], B)))
      X1 <- rep(w, B) - U

      X <- cbind(X, X1)
      mu.e1 <- c(mu.e1, mean(X1))
    }

    sigma.sq.e1 <- diag(sigma.sq.e2)

    new.dat.temp <- as.matrix(mod$data[, -c(1, which(names(mod$data) %in% names.w))])
    new.dat.temp0 <- which(!(names(mod$data) %in% names.w))[-1]

    if (length(new.dat.temp0) != 0) {
      new.dat <- as.data.frame(cbind(bigY, matrix(NA, nrow = B*n, ncol = (ncol(mod$data) - 1))))
      new.dat[, which(names(mod$data) %in% names.w)] <- X

      for(j in 1:length(new.dat.temp0)) {
        W2a <- new.dat.temp[, j]

        names.w3 <- names(mod$data)[new.dat.temp0]

        new.dat[, new.dat.temp0[j]] <- rep(W2a, B)
        names(new.dat)[new.dat.temp0[j]] <- names.w3[j]
      }
    }

    if (length(new.dat.temp0) == 0) new.dat <- as.data.frame(cbind(bigY, X))

    names(new.dat)[c(1, which(names(mod$data) %in% names.w))] <- names(mod$data)[c(1, which(names(mod$data) %in% names.w))]
  }

  # MC and E-step.

  while(cond) {

    # MC and E-step.

    if (is.matrix(sigma.sq.u) == F) prX <- stats::dnorm(X, mu.e1, sd = sqrt(sigma.sq.e1))
    if (is.matrix(sigma.sq.u) == T) prX <- mvtnorm::dmvnorm(X, mu.e1, sigma = sqrt(sigma.sq.e1))

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

    new.dat$weights1 <- weights1

    # M-step (updates).

    if (family == "gaussian") {
      mod <- stats::lm(formula = form.name, data = new.dat, weights = weights1, ...)
      sigma.sq.est <- (summary(mod)$sigma)^2
    }
    if (family == "binomial") mod <- suppressWarnings(stats::glm(formula = form.name, family = "binomial", data = new.dat, weights = weights1, ...))
    if (family == "poisson") mod <- suppressWarnings(stats::glm(formula = form.name, family = "poisson", data = new.dat, weights = weights1))
    if (family == "Gamma") mod <- suppressWarnings(stats::glm(formula = form.name, family = Gamma(link = "log"), data = new.dat, weights = weights1, ...))
    if (family == "negbin") mod <- suppressWarnings(MASS::glm.nb(formula = form.name, init.theta = theta.est, data = new.dat, weights = weights1, ...))

    beta.update <- stats::coef(mod)

    if (family == "negbin") theta.update <- mod$theta
    if (family == "Gamma") shape.update <- summary(mod)[14]$dispersion

    muPred <- stats::predict(mod, type = "response")

    if (is.matrix(sigma.sq.u) == F) {
      sigma.sq.e1.update <- wt.var(X, w = weights2)
      mu.e1.update <- stats::weighted.mean(X, w = weights2)
    }

    if (is.matrix(sigma.sq.u) == T) {
      sigma.sq.e1.update <- c()
      mu.e1.update <- c()

      for (kk in 1:q1) {
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
      if (diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon) {
        cond <- FALSE
        if (silent == FALSE) {
          message("convergence :-)")
          print(reps)
        }
        if (silent == TRUE) reps
        break
      }
    }

    if (family == "negbin") {
      if (diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon & theta.norm < epsilon) {
        cond <- FALSE
        if (silent == FALSE) {
          message("convergence :-)")
          print(reps)
        }
        if (silent == TRUE) reps
        break
      }
    }

    if (family == "Gamma") {
      if (diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon & shape.norm < epsilon) {
        cond <- FALSE
        if (silent == FALSE) {
          message("convergence :-)")
          print(reps)
        }
        if (silent == TRUE) reps
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

  mod_n$sigma.sq.u <- sigma.sq.u

  mod_n$B <- B

  names(beta.est) <- names(stats::coef(mod_n))
  mod_n$coefficients <- beta.est

  if (class(W[, 1])[1] == "poly") {
    mod_n$linear.predictors <- eta <- stats::predict(mod_n, newdata = eval(stats::getCall(mod_n)$data, environment(stats::formula(mod_n))))
    mod_n$fitted.values <- stats::predict(mod_n, type = "response", newdata = eval(stats::getCall(mod_n)$data, environment(stats::formula(mod_n))))
  }
  if (class(W[, 1])[1] != "poly") {
    mod_n$linear.predictors <- eta <- stats::predict(mod_n, newdata = mod_n$model)
    mod_n$fitted.values <- stats::predict(mod_n, type = "response", newdata = mod_n$model)
  }

  if (family == "gaussian") residuals <- Y - mod_n$linear.predictors

  if (family != "gaussian") {
    mu <- mod_n$family$linkinv(eta)
    mu.eta <- mod_n$family$mu.eta(eta)
    residuals <- (Y - mu)/mod_n$family$mu.eta(eta)
  }

  mod_n$residuals <- residuals

  qq <- mod_n$rank

  sumW <- apply(bigW, 1, sum, na.rm = T)
  weights1 <- as.vector(bigW)/sumW
  weights1[is.nan(weights1)] <- 0

  entropy <- sum(weights1*log(weights1), na.rm = TRUE)
  mod_n$entropy <- entropy/B

  if (family %in% c("gaussian", "Gamma")) qq <- qq + 1

  if (family != "gaussian") {
    mod_n$deviance <- mod$deviance + 2*mod_n$entropy
    mod_n$aic <- mod_n$deviance + 2*qq
    logLik.value <- mod_n$deviance/(-2)
  }

  if (family == "gaussian") {
    logLik.value <- 0.5*(sum(log(p.wt)) - n*(log(2*pi) + 1 - log(n) + log(sum(p.wt*residuals^2)))) - mod_n$entropy
    mod_n$aic <- -2*logLik.value + 2*qq
  }

  class(logLik.value) <- "logLik"
  mod_n$logLik <- logLik.value

  eff.samp.size <- 1/apply((bigW/sumW)^2, 1, sum)
  eff.samp.size <- as.numeric(eff.samp.size)
  mod_n$eff.samp.size <- eff.samp.size

  # Standard error calculations start here.

  estfun_mat <- sandwich::estfun(mod)
  if (family == "gaussian") estfun_mat <- sandwich::estfun(mod)*B

  K1 <- length(beta.update)

  S_1 <- matrix(0, nrow = K1, ncol = K1)
  SS_1 <- matrix(0, nrow = K1, ncol = K1)

  ind_mat <- matrix(1:(n*B), ncol = n, byrow = T)

  for (ii in 1:n) {
    index_vec <- ind_mat[, ii]
    S_1 <- S_1 + (apply(estfun_mat[index_vec, ], 2, sum))%*%t(apply(estfun_mat[index_vec, ], 2, sum))
  }

  if (family == "gaussian") {
    sand1 <- (sandwich::estfun(mod)*B/mod$weights)
    sand1[is.nan(sand1)] <- 1
    SS_1 <- t(sandwich::estfun(mod)*B)%*%sand1
    u.bar <- solve(stats::vcov(mod)*B)
    beta.est.se2 <- sqrt.na(diag(solve(u.bar - SS_1/B^2 + S_1/B^2)))
    AA <- expm::sqrtm((u.bar - SS_1/B^2 + S_1/B^2)*(stats::summary.lm(mod_n)$sigma)^2)
  }

  if (family != "gaussian") {
    sand1 <- (sandwich::estfun(mod)/mod$prior)
    sand1[is.nan(sand1)] <- 1
    SS_1 <- t(sandwich::estfun(mod))%*%sand1
    u.bar <- solve(stats::vcov(mod))
    beta.est.se2 <- sqrt.na(diag(solve(u.bar - SS_1 + S_1)))
    AA <- expm::sqrtm((u.bar - SS_1 + S_1)*stats::summary.glm(mod_n)$dispersion)
  }

  if (length(which(is.nan(beta.est.se2))) > 0) beta.est.se2 <- c(rep(NA, K1))

  mod_n$se <- beta.est.se2

  mod_n$qr <- qr(AA)
  dim_temp <- n - mod_n$rank
  mod_n$qr$qr <- rbind(qr(AA)$qr, matrix(rep(0, dim_temp), nrow = dim_temp, ncol = ncol(mod_n$qr$qr)))
  effects_names <- names(mod_n$effects)
  mod_n$effects <- mod$effects[1:n]
  names(mod_n$effects) <- effects_names

  class(mod_n) <- c("refitME", class(mod_n))

  mod_n$call <- match.call()

  return(mod_n)
}

#' Function for wrapping the MCEM algorithm on \code{gam} objects
#'
#' Function for wrapping the MCEM algorithm on GAMs where predictors are subject to measurement error/error-in-variables.
#' @name MCEMfit_gam
#' @param mod : a \code{gam} object (this is the naive fitted model). Make sure the first \eqn{p} input predictor variables entered in the naive model are the specified error-contaminated variables. These \eqn{p} predictors also need the measurement error variance to be specified in \code{sigma.sq.u}, see below.
#' @param family : a specified family/distribution.
#' @param sigma.sq.u : measurement error (ME) variance. A scalar if there is only one error-contaminated predictor variable, otherwise this must be stored as a vector (of ME variances) or a matrix if the ME covariance matrix is known.
#' @param B : the number of Monte Carlo replication values (default is set to 50).
#' @param epsilon : convergence threshold (default is set to 0.00001).
#' @param silent : if \code{TRUE}, the convergence message (which tells the user if the model has converged and reports the number of iterations required) is suppressed (default is set to \code{FALSE}).
#' @param ... : further arguments passed to \code{gam}.
#' @return \code{MCEMfit_gam} returns the original naive fitted model object but coefficient estimates and the covariance matrix have been replaced with the final MCEM model fit. Standard errors and the effective sample size (which diagnose how closely the proposal distribution matches the posterior, see equation (2) of Stoklosa, Hwang and Warton) have also been included as outputs.
#' @author Jakub Stoklosa, Wen-Han Hwang and David I. Warton.
#' @references Ganguli, B, Staudenmayer, J., and Wand, M. P. (2005). Additive models with predictors subject to measurement error. \emph{Australian & New Zealand Journal of Statistics}, \strong{47}, 193–202.
#' @references Wand, M. P. (2018). \pkg{SemiPar}: Semiparametic Regression. \proglang{R} package version 1.0-4.2., URL \url{https://CRAN.R-project.org/package=SemiPar}.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. \pkg{refitME}: Measurement Error Modelling using Monte Carlo Expectation Maximization in \proglang{R}.
#' @importFrom stats Gamma
#' @importFrom mgcv nb
#' @importFrom mvtnorm dmvnorm
#' @importFrom sandwich estfun
#' @export
#' @seealso \code{\link{MCEMfit_glm}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown vignette with examples.
#' @examples # A GAM example using the air pollution data set from the SemiPar package.
#'
#' library(refitME)
#' library(SemiPar)
#' library(mgcv)
#' library(dplyr)
#'
#' data(milan.mort)
#'
#' dat.air <- sample_n(milan.mort, 100) # Takes a random sample of size 100.
#'
#' Y <- dat.air[, 6]  # Mortality counts.
#'
#' n <- length(Y)
#'
#' z1 <- (dat.air[, 1])
#' z2 <- (dat.air[, 4])
#' z3 <- (dat.air[, 5])
#' w1 <- log(dat.air[, 9])  # The error-contaminated predictor (total suspended particles).
#'
#' dat <- data.frame(cbind(Y, w1, z1, z2, z3))
#'
#' gam_naiv <- gam(Y ~ s(w1), family = "poisson", data = dat)
#'
#' sigma.sq.u <- 0.0915 # Measurement error variance.
#'
#' B <- 10  # Consider increasing this if you want a more accurate answer.
#'
#' gam_MCEM <- refitME(gam_naiv, sigma.sq.u, B)
#'
#' plot(gam_MCEM, select = 1)
#'
#' detach(package:mgcv)
#'
MCEMfit_gam <- function(mod, family, sigma.sq.u, B = 50, epsilon = 0.00001, silent = FALSE, ...) {

  mod_n <- mod

  reps <- 0
  cond <- TRUE

  mod$data <- eval(stats::getCall(mod)$data, environment(stats::formula(mod)))

  Y <- stats::model.frame(mod)[, 1]
  bigY <- rep(Y, B)
  n <- length(Y)

  if (family == "gaussian") {
    if (is.null(stats::weights(mod))) p.wt <- rep(1, n)
    if (!is.null(stats::weights(mod))) p.wt <- stats::weights(mod)
  }
  if (family != "gaussian") p.wt <- stats::weights(mod)

  if (family == "negbin") theta.est <- mod$family$getTheta(TRUE)
  if (family == "Gamma") shape.est <- MASS::gamma.shape(mod)[1]$alpha

  bigp.wt <- rep(p.wt, B)

  W <- as.data.frame(stats::model.frame(mod)[, -1])
  names(W) <- names(stats::model.frame(mod))[-1]

  d <- ncol(W)

  muPred <- rep(stats::predict(mod, type = "response"), B)
  beta.est <- stats::coef(mod)

  form.name <- stats::formula(mod)

  if (is.matrix(sigma.sq.u) == F) {
    w <- W[, 1]

    names.w <- names(W)[1]

    sigma.sq.u1 <- sigma.sq.u
    sigma.sq.e1 <- stats::var(w) - sigma.sq.u1

    U <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1, B)))
    X <- rep(w, B) - U

    mu.e1 <- mean(X)

    new.dat.temp <- as.data.frame(mod$data[, -c(1, which(names(mod$data) %in% names.w))])

    new.dat.temp0 <- which(!(names(mod$data) %in% names.w))[-1]

    if (length(new.dat.temp0) != 0) {
      new.dat <- as.data.frame(cbind(bigY, matrix(NA, nrow = B*n, ncol = (ncol(mod$data) - 1))))
      new.dat[, which(names.w == names(mod$data))] <- X

      for (j in 1:length(new.dat.temp0)) {
        W2 <- new.dat.temp[, j]

        if (class(W2)[1] == "poly") W2 <- W2[, 1]

        new.dat[, new.dat.temp0[j]] <- rep(W2, B)
        names.w3 <- names(mod$data)[new.dat.temp0[j]]
        names(new.dat)[new.dat.temp0[j]] <- names.w3
      }
    }

    if (length(new.dat.temp0) == 0) new.dat <- as.data.frame(cbind(bigY, X))

    names(new.dat)[c(1, which(names(mod$data) %in% names.w))] <- names(mod$data)[c(1, which(names(mod$data) %in% names.w))]
  }

  if (is.matrix(sigma.sq.u) == T) {
    q1 <- dim(sigma.sq.u)[2]

    if (d < q1) stop("Number of error-contaminated covariates exceeds total number of covariates!")

    w <- W[, 1]

    sigma.sq.u1 <- sigma.sq.u
    sigma.sq.e2 <- stats::var(w) - sigma.sq.u1[1, 1]

    U <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1[1, 1], B)))
    X <- cbind(rep(w, B) - U)

    mu.e1 <- mean(X)

    for (kk in 2:q1) {
      w <- W[, kk]

      sigma.sq.e <- stats::var(w) - sigma.sq.u[kk, kk]
      sigma.sq.e2 <- c(sigma.sq.e2, sigma.sq.e)

      U <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1[kk, kk], B)))
      X1 <- rep(w, B) - U

      X <- cbind(X, X1)
      mu.e1 <- c(mu.e1, mean(X1))
    }

    sigma.sq.e1 <- diag(sigma.sq.e2)

    names.w <- names(W)[1:q1]

    new.dat.temp <- as.matrix(mod$data[, -c(1, which(names(mod$data) %in% names.w))])
    new.dat.temp0 <- which(!(names(mod$data) %in% names.w))[-1]

    if (length(new.dat.temp0) != 0) {
      new.dat <- as.data.frame(cbind(bigY, matrix(NA, nrow = B*n, ncol = (ncol(mod$data) - 1))))
      new.dat[, which(names(mod$data) %in% names.w)] <- X

      for(j in 1:length(new.dat.temp0)) {
        W2a <- new.dat.temp[, j]

        names.w3 <- names(mod$data)[new.dat.temp0]

        new.dat[, new.dat.temp0[j]] <- rep(W2a, B)
        names(new.dat)[new.dat.temp0[j]] <- names.w3[j]
      }
    }

    if (length(new.dat.temp0) == 0) new.dat <- as.data.frame(cbind(bigY, X))

    names(new.dat)[c(1, which(names(mod$data) %in% names.w))] <- names(mod$data)[c(1, which(names(mod$data) %in% names.w))]
  }

  while(cond) {

    # MC and E-step.

    if (is.matrix(sigma.sq.u) == F) prX <- stats::dnorm(X, mu.e1, sd = sqrt(sigma.sq.e1))
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

    new.dat$weights1 <- weights1

    # M-step (updates).

    if (family == "gaussian") {
      mod <- mgcv::gam(formula = form.name, data = new.dat, weights = weights1, ...)
      sigma.sq.est <- (summary(mod)$dispersion)
    }
    if (family == "binomial") mod <- suppressWarnings(mgcv::gam(formula = form.name, family = "binomial", gamma = 1.4, data = new.dat, weights = weights1, ...))
    if (family == "poisson") mod <- suppressWarnings(mgcv::gam(formula = form.name, family = "poisson", data = new.dat, weights = weights1, ...))
    if (family == "Gamma") mod <- suppressWarnings(mgcv::gam(formula = form.name, family = Gamma(link = "log"), data = new.dat, weights = weights1, ...))
    if (family == "negbin") mod <- suppressWarnings(mgcv::gam(form.name, family = mgcv::nb(), data = new.dat, weights = weights1, ...))

    beta.update <- stats::coef(mod)

    if (family == "negbin") theta.update <- mod$family$getTheta(TRUE)
    if (family == "Gamma") shape.update <- MASS::gamma.shape(mod)[1]$alpha

    muPred <- mgcv::predict.gam(mod, type = "response")

    if (is.matrix(sigma.sq.u) == F) {
      sigma.sq.e1.update <- wt.var(X, w = weights2)
      mu.e1.update <- stats::weighted.mean(X, w = weights2)
    }

    if (is.matrix(sigma.sq.u) == T) {
      sigma.sq.e1.update <- c()
      mu.e1.update <- c()

      for (kk in 1:q1) {
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
        if (silent == FALSE) {
          message("convergence :-)")
          print(reps)
        }
        if (silent == TRUE) reps
        break
      }
    }

    if (family == "negbin") {
      if ((diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon & theta.norm < epsilon) | reps > 50) {
        cond <- FALSE
        if (silent == FALSE) {
          print("convergence :-)")
          print(reps)
        }
        if (silent == TRUE) reps
        break
      }
    }

    if (family == "Gamma") {
      if ((diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon & shape.norm < epsilon) | reps > 50) {
        if (silent == FALSE) {
          print("convergence :-)")
          print(reps)
        }
        if (silent == TRUE) reps
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

  mod$sigma.sq.u <- sigma.sq.u

  mod$B <- B

  beta.est <- stats::coef(mod)
  names(beta.est) <- names(stats::coef(mod_n))
  mod$coefficients <- beta.est

  #mod_n$linear.predictors <- eta <- stats::predict(mod_n, newdata = eval(stats::getCall(mod_n)$data, environment(stats::formula(mod_n))))
  #mod_n$fitted.values <- stats::predict(mod_n, type = "response", newdata = eval(stats::getCall(mod_n)$data, environment(stats::formula(mod_n)))

  sumW <- apply(bigW, 1, sum, na.rm = T)
  weights1 <- as.vector(bigW)/sumW
  weights1[is.nan(weights1)] <- 0

  entropy <- sum(weights1*log(weights1), na.rm = TRUE)
  mod$entropy <- entropy/B

  sc.p <- as.numeric(mod$scale.estimated)
  qq <- sum(mod$edf) + sc.p

  if (family %in% c("gaussian", "Gamma")) qq <- qq + mod$family$n.theta

  logLik.value1 <- qq - mod$aic/2
  logLik.value <- logLik.value1 - mod$entropy
  mod$aic <- -2*logLik.value + 2*qq

  class(logLik.value) <- "logLik"
  mod$logLik <- logLik.value

  eff.samp.size <- 1/apply((bigW/sumW)^2, 1, sum)
  eff.samp.size <- as.numeric(eff.samp.size)
  mod$eff.samp.size <- eff.samp.size

  # Standard error calculations start here.

  estfun_mat <- sandwich::estfun(mod)

  if (family == "gaussian") estfun_mat <- sandwich::estfun(mod)*B

  X <- stats::model.matrix(mod)
  K1 <- ncol(X)

  ind_mat <- matrix(1:(n*B), ncol = n, byrow = T)

  S_1 <- matrix(0, nrow = K1, ncol = K1)
  SS_1 <- matrix(0, nrow = K1, ncol = K1)

  for (ii in 1:n) {
    index_vec <- ind_mat[, ii]
    S_1 <- S_1 + (apply(estfun_mat[index_vec, ], 2, sum))%*%t(apply(estfun_mat[index_vec, ], 2, sum))
  }

  if (family == "gaussian") {
    u.bar <- solve(stats::vcov(mod)*B)
    mod$Vp <- solve(u.bar - SS_1/B^2 + S_1/B^2)
  }

  if (family == "binomial" | family == "poisson" | family == "Gamma" | family == "negbin") {
    sand1 <- (sandwich::estfun(mod)/mod$prior)
    sand1[is.nan(sand1)] <- 1
    SS_1 <- t(sandwich::estfun(mod))%*%sand1
    u.bar <- solve(stats::vcov(mod))
    mod$Vp <- solve(u.bar - SS_1 + S_1)
  }

  beta.est.se2 <- sqrt.na(diag(mod$Vp))

  if (length(which(is.nan(beta.est.se2))) > 0) beta.est.se2 <- c(rep(NA, K1))

  mod$se <- beta.est.se2

  return(mod)
}

#' Function for fitting any likelihood-based model using the MCEM algorithm
#'
#' Function for wrapping the MCEM algorithm on any likelihood-based model where predictors are subject to measurement error/error-in-variables.
#' @name MCEMfit_gen
#' @param mod : a model object (this is the naive fitted model). Make sure the first \eqn{p} input predictor variables entered in the naive model are the specified error-contaminated variables. These \eqn{p} predictors also need the measurement error variance to be specified in \code{sigma.sq.u}, see below.
#' @param family : a specified family/distribution.
#' @param sigma.sq.u : measurement error (ME) variance. A scalar if there is only one error-contaminated predictor variable, otherwise this must be stored as a vector (of ME variances) or a matrix if the ME covariance matrix is known.
#' @param B : the number of Monte Carlo replication values (default is set to 50).
#' @param epsilon : a set convergence threshold (default is set to 0.00001).
#' @param silent : if \code{TRUE}, the convergence message (which tells the user if the model has converged and reports the number of iterations required) is suppressed (default is set to \code{FALSE}).
#' @param theta.est : an initial value for the dispersion parameter (this is required for fitting negative binomial models).
#' @param shape.est : an initial value for the shape parameter (this is required for fitting gamma models).
#' @param ... : further arguments passed through to the function that was used to fit \code{mod}, that will be used in refitting. These need only be specified if making changes to the arguments as compared to the original call that produced \code{mod}.
#' @return \code{MCEMfit_gen} returns the original naive fitted model object but coefficient estimates and residuals have been replaced with the final MCEM model fit. Standard errors are included and returned, if \code{mod} is a class of object accepted by the \pkg{sandwich} package (such as \code{glm}, \code{gam}, \code{survreg} and many more).
#' @author Jakub Stoklosa, Wen-Han Hwang and David I. Warton.
#' @references Carroll, R. J., Ruppert, D., Stefanski, L. A., and Crainiceanu, C. M. (2006). \emph{Measurement Error in Nonlinear Models: A Modern Perspective.} 2nd Ed. London: Chapman & Hall/CRC.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. \pkg{refitME}: Measurement Error Modelling using Monte Carlo Expectation Maximization in \proglang{R}.
#' @importFrom stats Gamma
#' @importFrom mvtnorm dmvnorm
#' @importFrom MASS glm.nb
#' @importFrom sandwich estfun
#' @importFrom mgcv nb
#' @export
#' @seealso \code{\link{MCEMfit_glm}} and \code{\link{MCEMfit_gam}}
#'
MCEMfit_gen <- function(mod, family, sigma.sq.u, B = 50, epsilon = 0.00001, silent = FALSE, theta.est = 1, shape.est = 1, ...) {

  mod_n <- mod

  if (family == "gaussian" | family == "negbin") mod$data <- eval(stats::getCall(mod)$data, environment(stats::formula(mod)))

  if (family == "negbin") {
    mod$data <- eval(stats::getCall(mod)$data, environment(stats::formula(mod)))
    theta.est <- mod$theta
  }

  if (family == "Gamma") shape.est <- summary(mod)[14]$dispersion

  if (is.null(mod$data)) mod$data <- eval(stats::getCall(mod)$data, environment(stats::formula(mod)))

  reps <- 0
  cond <- TRUE

  Y <- stats::model.frame(mod)[, 1]
  bigY <- rep(Y, B)
  n <- length(Y)

  if (family == "gaussian") {
    if (is.null(stats::weights(mod))) p.wt <- rep(1, n)
    if (!is.null(stats::weights(mod))) p.wt <- stats::weights(mod)
  }
  if (family != "gaussian") p.wt <- stats::weights(mod)

  bigp.wt <- rep(p.wt, B)

  W <- as.data.frame(stats::model.frame(mod)[, -1])
  names(W) <- names(stats::model.frame(mod))[-1]

  d <- ncol(W)

  muPred <- rep(stats::predict(mod, type = "response"), B)
  beta.est <- stats::coef(mod)

  if (is.matrix(sigma.sq.u) == F) {
    w <- W[, 1]

    names.w <- names(W)[1]

    if (dim(as.matrix(W))[2] != 1) {
      if (class(stats::model.frame(mod)[, -1])[1] == "poly" | class(stats::model.frame(mod)[, -1][, 1])[1] == "poly") {
        names.w <- names(mod$data)[2]
        names(mod$data)[2] <- names.w
      }
    }

    if (class(w)[1] == "poly") w <- w[, 1]

    sigma.sq.u1 <- sigma.sq.u
    sigma.sq.e1 <- stats::var(w) - sigma.sq.u1

    U <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1, B)))
    X <- rep(w, B) - U

    mu.e1 <- mean(X)

    new.dat.temp <- as.data.frame(mod$data[, -c(1, which(names(mod$data) %in% names.w))])

    new.dat.temp0 <- which(!(names(mod$data) %in% names.w))[-1]

    if (length(new.dat.temp0) != 0) {
      new.dat <- as.data.frame(cbind(bigY, matrix(NA, nrow = B*n, ncol = (ncol(mod$data) - 1))))
      new.dat[, which(names.w == names(mod$data))] <- X

      for (j in 1:length(new.dat.temp0)) {
        W2 <- new.dat.temp[, j]

        if (class(W2)[1] == "poly") W2 <- W2[, 1]

        new.dat[, new.dat.temp0[j]] <- rep(W2, B)
        names.w3 <- names(mod$data)[new.dat.temp0[j]]
        names(new.dat)[new.dat.temp0[j]] <- names.w3
      }
    }

    if (length(new.dat.temp0) == 0) new.dat <- as.data.frame(cbind(bigY, X))

    names(new.dat)[c(1, which(names(mod$data) %in% names.w))] <- names(mod$data)[c(1, which(names(mod$data) %in% names.w))]
  }

  if (is.matrix(sigma.sq.u) == T) {
    q1 <- dim(sigma.sq.u)[2]

    names.w <- names(W)[1]

    w <- W[, 1]

    W1 <- W

    poly.trig <- FALSE

    if (class(stats::model.frame(mod)[, -1][, 1])[1] == "poly") {
      names.w <- names(mod$data)[2]
      W1 <- W[, 1]

      d <- ncol(W)
      poly.trig <- TRUE
    }

    if (class(w)[1] == "poly") {
      w <- w[, 1]
      poly.trig <- TRUE
    }

    d1 <- 0

    if (ncol(as.matrix(W1)) > 1) {
      if ((sum(W1[, 2] == w^2) == n)) {    # Check if quadratic.
        if (class(stats::model.frame(mod)[, -1][, 1])[1] != "poly") d <- d - 1
        d1 <- 1
      }
    }

    if (ncol(as.matrix(W1)) > 2) {    # Check if cubic.
      if (sum(W1[, 3] == w^3) == n) {
        if (class(stats::model.frame(mod)[, -1][, 1])[1] != "poly") d <- d - 1
        d1 <- 2
      }
    }
    if (ncol(as.matrix(W1)) > 3) {     # Check if quartic.
      if (sum(W1[, 4] == w^4) == n) {
        if (class(stats::model.frame(mod)[, -1][, 1])[1] != "poly") d <- d - 1
        d1 <- 3
      }
    }

    if (d < q1) stop("Number of error-contaminated covariates exceeds total number of covariates!")

    sigma.sq.u1 <- sigma.sq.u
    sigma.sq.e2 <- stats::var(w) - sigma.sq.u[1, 1]

    U <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1[1, 1], B)))
    X <- cbind(rep(w, B) - U)

    mu.e1 <- mean(X)

    for (kk in 2:q1) {

      if (poly.trig  == TRUE) {
        W2 <- stats::model.frame(mod)[, -c(1:kk)]

        if (ncol(as.matrix(W2)) == 1) {
          W1 <- w <- W[, ncol(W)]
          names.w1 <- names(mod$data)[length(names(mod$data))]

          poly.trig0 <- FALSE
        }

        if (ncol(as.matrix(W2)) > 1) {
          if (class(W2)[1] == "poly" | class(W2[, 1])[1] == "poly") {
            W1 <- W[, kk]
            w <- W1[, 1]

            poly.trig0 <- TRUE
          }

          if (class(W2)[1] != "poly" & class(W2[, 1])[1] != "poly") {
            W1 <- as.matrix(W[, -c(1:(kk + d1 - 1))])
            w <- as.matrix(W[, kk + d1])

            poly.trig0 <- FALSE
          }

          names.w1 <- names(mod$data)[kk + 1]
        }

        names.w <- c(names.w, names.w1)
      }

      if (poly.trig  == FALSE) {
        W2 <- stats::model.frame(mod)[, -c(1:(kk + d1))]

        if (ncol(as.matrix(W2)) == 1) {
          W1 <- w1 <- W[, ncol(W)]
          names.w1 <- colnames(W)[ncol(W)]

          poly.trig0 <- FALSE
        }

        if (ncol(as.matrix(W2)) > 1) {
          if (class(W2)[1] == "poly" | class(W2[, 1])[1] == "poly") {
            names.w1 <- names(mod$data)[kk + 1]
            W1 <- W[, kk + d1]
            w <- W1[, 1]

            poly.trig0 <- TRUE
          }

          if (class(W2)[1] != "poly" & class(W2[, 1])[1] != "poly") {
            W1 <- as.matrix(W2)
            w <- W2[, 1]
            names.w1 <- colnames(W1)[1]

            poly.trig0 <- FALSE
          }
        }

        names.w <- c(names.w, names.w1)
      }

      sigma.sq.e <- stats::var(w) - sigma.sq.u[kk, kk]
      sigma.sq.e2 <- c(sigma.sq.e2, sigma.sq.e)

      if (ncol(as.matrix(W1)) > 1) {
        if (sum(W1[, 2] == w^2) == n) d1 <- d1 + 1  # Check if quadratic.
      }

      if (ncol(as.matrix(W1)) > 2) {   # Check if cubic.
        if (sum(W1[, 3] == w^3) == n) d1 <- d1 + 2
      }

      if (ncol(as.matrix(W1)) > 3) {    # Check if quartic.
        if (sum(W1[, 4] == w^4) == n) d1 <- d1 + 3
      }

      poly.trig <- poly.trig0

      U <- stats::rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1[kk, kk], B)))
      X1 <- rep(w, B) - U

      X <- cbind(X, X1)
      mu.e1 <- c(mu.e1, mean(X1))
    }

    sigma.sq.e1 <- diag(sigma.sq.e2)

    new.dat.temp <- as.matrix(mod$data[, -c(1, which(names(mod$data) %in% names.w))])
    new.dat.temp0 <- which(!(names(mod$data) %in% names.w))[-1]

    if (length(new.dat.temp0) != 0) {
      new.dat <- as.data.frame(cbind(bigY, matrix(NA, nrow = B*n, ncol = (ncol(mod$data) - 1))))
      new.dat[, which(names(mod$data) %in% names.w)] <- X

      for(j in 1:length(new.dat.temp0)) {
        W2a <- new.dat.temp[, j]

        names.w3 <- names(mod$data)[new.dat.temp0]

        new.dat[, new.dat.temp0[j]] <- rep(W2a, B)
        names(new.dat)[new.dat.temp0[j]] <- names.w3[j]
      }
    }

    if (length(new.dat.temp0) == 0) new.dat <- as.data.frame(cbind(bigY, X))

    names(new.dat)[c(1, which(names(mod$data) %in% names.w))] <- names(mod$data)[c(1, which(names(mod$data) %in% names.w))]
  }

  while(cond) {

    # MC and E-step.

    if (is.matrix(sigma.sq.u) == F) prX <- stats::dnorm(X, mu.e1, sd = sqrt(sigma.sq.e1))
    if (is.matrix(sigma.sq.u) == T) prX <- mvtnorm::dmvnorm(X, mu.e1, sigma = sqrt(sigma.sq.e1))

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

    new.dat$weights1 <- weights1

    # M-step (updates).

    mod <- suppressWarnings(stats::update(mod, data = new.dat, weights = weights1, evaluate = TRUE, ...))

    beta.update <- stats::coef(mod)

    if (family == "negbin") theta.update <- mod$theta
    if (family == "Gamma") shape.update <- summary(mod)[14]$dispersion

    muPred <- stats::predict(mod, type = "response")

    if (is.matrix(sigma.sq.u) == F) {
      sigma.sq.e1.update <- wt.var(X, w = weights2)
      mu.e1.update <- stats::weighted.mean(X, w = weights2)
    }

    if (is.matrix(sigma.sq.u) == T) {
      sigma.sq.e1.update <- c()
      mu.e1.update <- c()

      for (kk in 1:q1) {
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
      if (diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon) {
        cond <- FALSE
        if (silent == FALSE) {
          message("convergence :-)")
          print(reps)
        }
        if (silent == TRUE) reps
        break
      }
    }

    if (family == "negbin") {
      if (diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon & theta.norm < epsilon) {
        cond <- FALSE
        if (silent == FALSE) {
          message("convergence :-)")
          print(reps)
        }
        if (silent == TRUE) reps
        break
      }
    }

    if (family == "Gamma") {
      if (diff.mu_e < epsilon & diff.sig_e < epsilon & beta.norm < epsilon & shape.norm < epsilon) {
        cond <- FALSE
        if (silent == FALSE) {
          message("convergence :-)")
          print(reps)
        }
        if (silent == TRUE) reps
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

  mod_n$sigma.sq.u <- sigma.sq.u

  mod_n$B <- B

  names(beta.est) <- names(stats::coef(mod_n))
  mod_n$coefficients <- beta.est

  if (class(W[, 1])[1] == "poly") {
    mod_n$linear.predictors <- eta <- stats::predict(mod_n, newdata = eval(stats::getCall(mod_n)$data, environment(stats::formula(mod_n))))
    mod_n$fitted.values <- stats::predict(mod_n, type = "response", newdata = eval(stats::getCall(mod_n)$data, environment(stats::formula(mod_n))))
  }

  if (class(W[, 1])[1] != "poly") {
    mod_n$linear.predictors <- eta <- stats::predict(mod_n, newdata = mod_n$model)
    mod_n$fitted.values <- stats::predict(mod_n, type = "response", newdata = mod_n$model)
  }
  sumW <- apply(bigW, 1, sum, na.rm = T)
  weights1 <- as.vector(bigW)/sumW

  entropy <- sum(weights1*log(weights1), na.rm = TRUE)
  mod_n$entropy <- entropy/B

  # Standard error calculations start here.

  flag.se <- TRUE

  estfun_mat <- try(sandwich::estfun(mod), silent = TRUE)
  if (class(estfun_mat)[1] == "try-error") {
    flag.se <- FALSE
    message("Note: The fitted model object is not supported by the sandwich R-package, so we cannot compute standard errors. Only parameter estimates will be provided in the output.")
  }

  if (flag.se == TRUE) {
    if (family == "gaussian") estfun_mat <- sandwich::estfun(mod)*B

    K1 <- length(beta.update)

    S_1 <- matrix(0, nrow = K1, ncol = K1)
    SS_1 <- matrix(0, nrow = K1, ncol = K1)

    ind_mat <- matrix(1:(n*B), ncol = n, byrow = T)

    for (ii in 1:n) {
      index_vec <- ind_mat[, ii]
      S_1 <- S_1 + (apply(estfun_mat[index_vec, ], 2, sum))%*%t(apply(estfun_mat[index_vec, ], 2, sum))
    }

    if (family == "gaussian") {
      sand1 <- (sandwich::estfun(mod)*B/mod$weights)
      sand1[is.nan(sand1)] <- 1
      SS_1 <- t(sandwich::estfun(mod)*B)%*%sand1
      u.bar <- solve(stats::vcov(mod)*B)
      beta.est.se2 <- sqrt.na(diag(solve(u.bar - SS_1/B^2 + S_1/B^2)))
    }

    if (family != "gaussian") {
      sand1 <- (sandwich::estfun(mod)/mod$prior)
      sand1[is.nan(sand1)] <- 1
      SS_1 <- t(sandwich::estfun(mod))%*%sand1
      u.bar <- solve(stats::vcov(mod))
      beta.est.se2 <- sqrt.na(diag(solve(u.bar - SS_1 + S_1)))
    }

    if (length(which(is.nan(beta.est.se2))) > 0) beta.est.se2 <- c(rep(NA, K1))

    mod_n$se <- beta.est.se2
  }

  class(mod_n) <- c("refitME", class(mod_n))

  mod_n$call <- match.call()

  return(mod_n)
}

#' Function for fitting \code{VGAM} capture-recapture (CR) model using the MCEM algorithm
#'
#' Function for fitting \code{VGAM} capture-recapture (CR) model using the MCEM algorithm where covariates have measurement error.
#' @name MCEMfit_CR
#' @section Warning:
#' This function is still under development. Currently the function can only fit the CR model used in the manuscript. IT DOES NOT SUPPORT ALL \code{VGAM} families.
#' @param mod : a \code{vglm/vgam} object (this is the naive CR model). Make sure the first \eqn{p} input predictor variables in the naive model are the selected error-contaminated variables.
#' @param sigma.sq.u : measurement error (ME) variance. A scalar if there is only one error-contaminated predictor variable, otherwise this must be stored as a vector (of ME variances) or a matrix if the ME covariance matrix is known.
#' @param B : the number of Monte Carlo replication values (default is set to 50).
#' @param epsilon : a set convergence threshold (default is set to 0.00001).
#' @param silent : if \code{TRUE}, the convergence message (which tells the user if the model has converged and reports the number of iterations required) is suppressed (default is set to \code{FALSE}).
#' @return \code{MCEMfit_CR} returns model coefficient and population size estimates with standard errors and the effective sample size.
#' @author Jakub Stoklosa, Wen-Han Hwang and David I. Warton.
#' @references Stoklosa, J., Hwang, W-H., and Warton, D.I. \pkg{refitME}: Measurement Error Modelling using Monte Carlo Expectation Maximization in \proglang{R}.
#' @importFrom VGAM s posbinomial
#' @importFrom VGAMdata dposbinom
#' @export
#' @seealso \code{\link{MCEMfit_glm}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown vignette with examples.
#' @examples # A VGAM example using the Prinia flaviventris capture-recapture data.
#'
#' library(refitME)
#' library(VGAM)
#'
#' data(Priniadata)
#'
#' tau <- 17   # No. of capture occasions.
#' w1 <- Priniadata$w1 # Bird wing length predictor.
#'
#' CR_naiv <- vglm(cbind(cap, noncap) ~ w1,
#'    VGAM::posbinomial(omit.constant = TRUE, parallel = TRUE ~ w1),
#'    data = Priniadata, trace = FALSE)
#'
#' sigma.sq.u <- 0.37 # ME variance.
#'
#' CR_MCEM <- refitME(CR_naiv, sigma.sq.u)
#'
#' detach(package:VGAM)
#'
MCEMfit_CR <- function(mod, sigma.sq.u, B = 50, epsilon = 0.00001, silent = FALSE) {

  reps <- 0
  cond <- TRUE

  tau <- mod@extra$tau[1]

  mod.terms <- attr(mod@terms$terms, "term.labels")

  d <- length(mod.terms)

  n <- mod@misc$n

  beta.est <- stats::coef(mod)

  w1 <- mod@x[, 2]

  W1 <- stats::model.matrix(mod)[, -1]

  if (d == 1 & mod@misc$function.name == "vglm") p1 <- 1

  if (d == 1 & mod@misc$function.name == "vgam") p1 <- "spline"

  if (d > 1 & mod@misc$function.name == "vglm") {
    if (sum((w1)^2 != (W1[, 2])) == n) p1 <- 1
    if (sum((w1)^2 == (W1[, 2])) == n) p1 <- 2
    if (d > 2) {
      if ((sum((w1)^2 == (W1[, 2])) == n) & sum((w1)^3 == (W1[, 3])) == n) p1 <- 3
      if ((sum((w1)^2 == (W1[, 2])) == n) & sum((w1)^3 != (W1[, 3])) == n) p1 <- 2
    }
  }

  sigma.sq.u1 <- sigma.sq.u
  sigma.sq.e1 <- stats::var(mod@x[, 2]) - sigma.sq.u

  U1_j <- stats::rnorm(n*B, 0, sd = sqrt(sigma.sq.u1))
  X1_j <- rep(w1, B) - U1_j

  mu.e1 <- mean(X1_j)

  if (p1 == 1) {
    X <- cbind(rep(1, B*n), X1_j)
    muPred <- rep(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))), B)
  }

  if (p1 == 2) {
    X <- cbind(rep(1, B*n), X1_j, (X1_j)^2)
    muPred <- rep(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))), B)
  }

  if (p1 == 3) {
    X <- cbind(rep(1, B*n), X1_j, (X1_j)^2, (X1_j)^3)
    muPred <- rep(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))), B)
  }

  if (p1 == "spline") {
    X <- cbind(rep(1, B*n), X1_j)
    muPred <- rep(1/(1 + exp(-VGAM::predict.vgam(mod, type = "link"))), B)
  }

  if (p1 == 1) colnames(X) <- c(names(stats::coef(mod))[1], "x1")
  if (p1 == 2) colnames(X) <- c(names(stats::coef(mod))[1], "x1", "x2")
  if (p1 == 3) colnames(X) <- c(names(stats::coef(mod))[1], "x1", "x2", "x3")

  bigY <- rep(mod@y*tau, B)

  CR_dat <- data.frame(cbind(X1_j, bigY, tau - bigY))
  colnames(CR_dat) <- c("x1", "cap", "noncap")

  while(cond) {

    # MC and E-step.

    prX <- stats::dnorm(X1_j, mu.e1, sd = sqrt(sigma.sq.e1))
    prY <- VGAMdata::dposbinom(bigY, tau, muPred)

    # M-step (updates).

    bigW <- matrix(prY*prX, n, B)

    sumW <- rep(apply(bigW, 1, sum), B)

    weights1 <- as.vector(bigW)/sumW

    if (p1 == 1) mod <- suppressWarnings(VGAM::vglm(cbind(cap, noncap) ~ x1, VGAM::posbinomial(omit.constant = TRUE, parallel = TRUE ~ x1), data = CR_dat, trace = F, weights = weights1))
    if (p1 == 2) mod <- suppressWarnings(VGAM::vglm(cbind(cap, noncap) ~ x1 + I(x1^2), VGAM::posbinomial(omit.constant = TRUE, parallel = TRUE ~ x1 + I(x1^2)), data = CR_dat, trace = F, weights = weights1))
    if (p1 == 3) mod <- suppressWarnings(VGAM::vglm(cbind(cap, noncap) ~ x1 + I(x1^2) + I(x1^3), VGAM::posbinomial(omit.constant = TRUE, parallel = TRUE ~ x1 + I(x1^2)), data = CR_dat, trace = F, weights = weights1))
    if (p1 == "spline") mod <- suppressWarnings(VGAM::vgam(cbind(cap, noncap) ~ s(x1, df = 2), VGAM::posbinomial(omit.constant = TRUE, parallel = TRUE ~ s(x1, df = 2)), data = CR_dat, trace = F, weights = weights1))

    beta.update <- stats::coef(mod)
    if (p1 == 1 | p1 == 2 | p1 == 3) muPred <- 1/(1 + exp(-VGAM::predictvglm(mod, type = "link")))
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
      if (silent == FALSE) {
        print("convergence :-)")
        print(reps)
      }
      if (silent == TRUE) reps
      break
    }

    # Update parameters.

    beta.est <- beta.update
    sigma.sq.e1 <- sigma.sq.e1.update
    mu.e1 <- mu.e1.update
  }

  beta.est.se1 <- sqrt.na(diag(VGAM::vcov(mod))) # Naive SE estimator.

  sumW <- apply(bigW, 1, sum, na.rm = T)
  weights1 <- as.vector(bigW)/sumW
  weights1[is.nan(weights1)] <- 0
  entropy <- sum(weights1*log(weights1), na.rm = TRUE)/B

  aic.value <- -2*(stats::logLik(mod) - entropy) + (VGAM::AIC(mod) + 2*stats::logLik(mod))

  eff.samp.size <- 1/apply((bigW/sumW)^2, 1, sum)
  eff.samp.size <- as.numeric(eff.samp.size)

  if (p1 == 1 | p1 == 2 | p1 == 3) pr.est <- 1/(1 + exp(-VGAM::predictvglm(mod, type = "link")))

  if (p1 == "spline") pr.est <- 1/(1 + exp(-VGAM::predict.vgam(mod, type = "link")))

  pi.est <- 1 - (1 - pr.est)^tau

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

  for (iii in 1:n) {
    index_vec <- ind_mat[, iii]
    x <- X[index_vec, ]

    if (p1 == 1 | p1 == 2 | p1 == 3) {
      SS_1 <- SS_1 + t(x*(((bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))))[index_vec]/pi.est[index_vec]))*VGAM::weights(mod, type = "working")[index_vec]))%*%(x*((bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))))[index_vec]/pi.est[index_vec])))
      S_1 <- S_1 + (apply(x*(bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))))[index_vec]/pi.est[index_vec])*VGAM::weights(mod, type = "working")[index_vec], 2, sum))%*%t(apply(x*(bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predictvglm(mod, type = "link"))))[index_vec]/pi.est[index_vec])*VGAM::weights(mod, type = "working")[index_vec], 2, sum))
    }

    if (p1 == "spline") {
      SS_1 <- SS_1 + t(x*(((bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predict.vgam(mod, type = "link"))))[index_vec]/pi.est[index_vec]))*VGAM::weights(mod, type = "working")[index_vec]))%*%(x*((bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predict.vgam(mod, type = "link"))))[index_vec]/pi.est[index_vec])))
      S_1 <- S_1 + (apply(x*(bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predict.vgam(mod, type = "link"))))[index_vec]/pi.est[index_vec])*VGAM::weights(mod, type = "working")[index_vec], 2, sum))%*%t(apply(x*(bigY[index_vec] - tau*(1/(1 + exp(-VGAM::predict.vgam(mod, type = "link"))))[index_vec]/pi.est[index_vec])*VGAM::weights(mod, type = "working")[index_vec], 2, sum))
    }
  }

  u.bar <- solve(VGAM::vcov(mod))
  beta.est.se <- sqrt.na(diag(solve(u.bar - SS_1 + S_1)))

  if (length(which(is.nan(beta.est.se))) > 0) beta.est.se2 <- c(rep(NA, K1))

  idG <- solve(u.bar - SS_1 + S_1)
  X.s <- stats::model.matrix(mod)

  dpi <- w.est*tau*pr.est*(1 - pi.est)/pi.est/pi.est
  dNhat <- apply(X.s*c(dpi), 2, sum)
  N.est.se <- sqrt(sum(w.est*(1 - pi.est)/pi.est/pi.est) + t(dNhat)%*%idG%*%dNhat)

  values <- list(beta = beta.est, beta.se = beta.est.se, N.est = N.est, N.est.se = N.est.se, mod = mod, aic.value = aic.value, eff.samp.size = eff.samp.size)

  return(values)
}

#' Function that calculates a weighted variance
#'
#' This function that calculates a weighted variance for a given vector.
#' @name wt.var
#' @param x : a vector of numerical data.
#' @param w : a vector of equal length to \code{x} representing the weights.
#' @return \code{wt.var} returns a single value from analysis requested.
#' @examples # Define simple data
#' x = 1:25 # Set of numbers.
#' wt = runif(25) # Some arbitrary weights.
#'
#' # Display variances (unweighted and then weighted).
#' var(x)
#' wt.var(x, wt)
#' @export
#' @source The developer of this function is Jeremy VanDerWal. See \url{https://rdrr.io/cran/SDMTools/src/R/wt.mean.R}
#'
wt.var <- function(x, w) {
  s <- which(is.finite(x + w))
  wt <- w[s]
  x <- x[s] # Remove NA info.
  xbar <- stats::weighted.mean(x, w = wt) # Get the weighted mean.

  return(sum(wt*(x - xbar)^2)*(sum(wt)/(sum(wt)^2 - sum(wt^2)))) # Return the variance.
}

#' Function that replaces NA with zero for a matrix
#'
#' This function replaces NA with zero for a matrix.
#' @name sqrt.na
#' @param x : a matrix
#' @return \code{sqrt.na} returns a matrix.
#' @author Jakub Stoklosa
#' @export
#'
sqrt.na <- function(x) {
  id <- (x > 0) & (is.na(x) == F)
  sx <- rep(0, length(x))
  if (length(dim(x)) == 2) sx <- matrix(sx, dim(x))
  sx[id] <- sqrt(x[id])
  sx
}

#' An ANOVA function for fitted \code{refitME} objects
#'
#' An ANOVA function for fitted \code{refitME} objects.
#' @name anova.refitME
#' @param object : fitted model objects of class \code{refitME}.
#' @param ... : further arguments passed through to \code{lm} or \code{glm}.
#' @param dispersion : the dispersion parameter for the fitting family. By default it is obtained from the object(s).
#' @param test : a character string, (partially) matching one of "\code{Chisq}", "\code{LRT}", "\code{Rao}", "\code{F}" or "\code{Cp}". See \code{\link{stat.anova}}.
#' @return \code{anova.refitME} produces output identical to \code{anova.lm}, \code{anova.glm} or \code{anova.gam}.
#' @author Jakub Stoklosa, Wen-Han Hwang and David I. Warton.
#' @importFrom stats stat.anova pf
#' @export
#' @seealso \code{\link{anova}}
#'
anova.refitME <- function(object, ..., dispersion = NULL, test = NULL) {
  anova.lmlist <- function (object, ..., scale = 0, test = "F") {
    objects <- list(object, ...)
    responses <- as.character(lapply(objects, function(x) deparse(x$terms[[2L]])))
    sameresp <- responses == responses[1L]
    if (!all(sameresp)) {
      objects <- objects[sameresp]
      warning(gettextf("models with response %s removed because response differs from model 1",
                       sQuote(deparse(responses[!sameresp]))), domain = NA)
    }
    ns <- sapply(objects, function(x) length(x$residuals))
    if (any(ns != ns[1L]))
      stop("models were not all fitted to the same size of dataset")
    nmodels <- length(objects)
    if (nmodels == 1)
      return(anova.lm(object))
    resdf <- as.numeric(lapply(objects, stats::df.residual))
    resdev <- as.numeric(lapply(objects, stats::deviance))
    table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA,
                                                              -diff(resdev)))
    variables <- lapply(objects, function(x) paste(deparse(stats::formula(x)),
                                                   collapse = "\n"))
    dimnames(table) <- list(1L:nmodels, c("Res.Df", "RSS",
                                          "Df", "Sum of Sq"))
    title <- "Analysis of Variance Table\n"
    topnote <- paste0("Model ", format(1L:nmodels), ": ",
                      variables, collapse = "\n")
    if (!is.null(test)) {
      bigmodel <- order(resdf)[1L]
      scale <- if (scale > 0)
        scale
      else resdev[bigmodel]/resdf[bigmodel]
      table <- stats::stat.anova(table = table, test = test, scale = scale,
                          df.scale = resdf[bigmodel], n = length(objects[[bigmodel]]$residuals))
    }
    structure(table, heading = c(title, topnote), class = c("anova",
                                                            "data.frame"))
  }

  qr.lm <- function (x, ...) {
    if (is.null(r <- x$qr))
      stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
    r
  }

  anova.lm <- function (object, ...) {
    if (length(list(object, ...)) > 1L)
      return(anova.lmlist(object, ...))
    if (!inherits(object, "lm"))
      warning("calling anova.lm(<fake-lm-object>) ...")
    w <- object$weights
    ssr <- sum(if (is.null(w)) object$residuals^2 else w * object$residuals^2)
    mss <- sum(if (is.null(w)) object$fitted.values^2 else w *
                 object$fitted.values^2)
    if (ssr < 1e-10 * mss)
      warning("ANOVA F-tests on an essentially perfect fit are unreliable")
    dfr <- stats::df.residual(object)
    p <- object$rank
    if (p > 0L) {
      p1 <- 1L:p
      comp <- object$effects[p1]
      asgn <- object$assign[qr.lm(object)$pivot][p1]
      nmeffects <- c("(Intercept)", attr(object$terms,
                                         "term.labels"))
      tlabels <- nmeffects[1 + unique(asgn)]
      ss <- c(vapply(split(comp^2, asgn), sum, 1), ssr)
      df <- c(lengths(split(asgn, asgn)), dfr)
    }
    else {
      ss <- ssr
      df <- dfr
      tlabels <- character()
    }
    ms <- ss/df
    f <- ms/(ssr/dfr)
    P <- pf(f, df, dfr, lower.tail = FALSE)
    table <- data.frame(df, ss, ms, f, P)
    table[length(P), 4:5] <- NA
    dimnames(table) <- list(c(tlabels, "Residuals"), c("Df",
                                                       "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
    if (attr(object$terms, "intercept"))
      table <- table[-1, ]
    structure(table, heading = c("Analysis of Variance Table\n",
                                 paste("Response:", deparse(stats::formula(object)[[2L]]))),
              class = c("anova", "data.frame"))
  }

  classObject <- class(object)[2]
  anovaFn <- switch(classObject, "lm" = anova.lm, "glm" = anova_MCEMfit_glm, ...)

  return(anovaFn(object))
}

#' Extract log-Likelihoods for \code{refitME} model objects
#'
#' Extract log-Likelihoods for \code{refitME} model objects. This function subtracts the entropy term from the observed likelihood.
#' @name logLik.refitME
#' @param object : fitted model objects of class \code{refitME}.
#' @param ... : further arguments passed through to \code{lm} or \code{glm}.
#' @return \code{logLik.refitME} produces identical output to \code{logLik} but for \code{refitME} model objects.
#' @author Jakub Stoklosa, Wen-Han Hwang and David I. Warton.
#' @importFrom stats family
#' @export
#' @seealso \code{\link{logLik}}
#'
logLik.refitME <- function(object, ...) {
  logLik.glm <- function (object, ...) {
    if (!missing(...))
      warning("extra arguments discarded")
    fam <- stats::family(object)$family
    p <- object$rank
    if (fam %in% c("gaussian", "Gamma", "inverse.gaussian"))
      p <- p + 1

    val <- p - object$aic/2
    attr(val, "nobs") <- sum(!is.na(object$residuals))
    attr(val, "df") <- p
    class(val) <- "logLik"
    val
  }

  classObject <- class(object)[2]
  logLikFn <- switch(classObject, "lm" = logLik_MCEMfit_lm, "glm" = logLik.glm, ...)

  return(logLikFn(object))
}

#' An ANOVA function for fitted \code{MCEMfit_glm} objects
#'
#' An ANOVA function for fitted \code{MCEMfit_glm} objects.
#' @name anova_MCEMfit_glm
#' @param object : fitted model objects of class \code{MCEMfit_glm}.
#' @param ... : further arguments passed through to \code{glm}.
#' @param dispersion : the dispersion parameter for the fitting family. By default it is obtained from the object(s).
#' @param test : a character string, (partially) matching one of "\code{Chisq}", "\code{LRT}", "\code{Rao}", "\code{F}" or "\code{Cp}".
#' @return \code{anova_MCEMfit_glm} produces output identical to \code{anova.glm}.
#' @author Jakub Stoklosa, Wen-Han Hwang and David I. Warton.
#' @export
#' @seealso \code{\link{anova.glm}}
#'
anova_MCEMfit_glm <- function(object, ..., dispersion = NULL, test = NULL) {
  B <- object$B
  family <- object$family$family
  sigma.sq.u <- object$sigma.sq.u
  dotargs <- list(...)

  named <- if (is.null(names(dotargs))) rep_len(FALSE, length(dotargs))
  else (names(dotargs) != "")

  if (any(named)) warning("the following arguments to 'anova.refitME.glm' are invalid and dropped: ", paste(deparse(dotargs[named]), collapse = ", "))

  dotargs <- dotargs[!named]
  is.glm <- vapply(dotargs, function(x) inherits(x, "glm"), NA)
  dotargs <- dotargs[is.glm]

  doscore <- !is.null(test) && test == "Rao"
  varlist <- attr(object$terms, "variables")

  x <- if (n <- match("x", names(object), 0L)) object[[n]]
  else stats::model.matrix(object)

  varseq <- attr(x, "assign")
  nvars <- max(0, varseq)
  resdev <- resdf <- NULL

  if (doscore) {
    score <- numeric(nvars)
    y <- object$y
    fit <- stats::glm(y ~ 1, family = family)
    r <- fit$residuals
    w <- fit$weights
    icpt <- attr(object$terms, "intercept")
  }

  if (nvars > 1 || doscore) {
    method <- "MCEMfit_glm"
    y <- object$y

    if (is.null(y)) {
      mu.eta <- object$family$mu.eta
      eta <- object$linear.predictors
      y <- object$fitted.values + object$residuals*mu.eta(eta)
    }

    for (i in seq_len(max(nvars - 1L, 0))) {
      x1 <- x[, varseq <= i, drop = FALSE]
      dat.n <- as.data.frame(cbind(y, x1[, -1]))
      colnames(dat.n) <- c("y", colnames(x1)[-1])
      mod <- stats::glm(y ~ ., family = family, data = dat.n)
      fit <- eval(call(if (is.function(method)) "method" else method, mod = mod, family = family, sigma.sq.u = sigma.sq.u, B = B, silent = TRUE))

      if (doscore) {
        x1 <- x[, varseq <= i, drop = FALSE]
        dat.n <- as.data.frame(cbind(y, x1[, -1]))
        colnames(dat.n) <- c("y", colnames(x1)[-1])
        mod <- stats::glm(y ~ ., family = family, data = dat.n)
        zz <- eval(call(if (is.function(method)) "method" else method, mod = mod, family = family, sigma.sq.u = sigma.sq.u, B = B, silent = TRUE, y = r, weights = w, intercept = icpt))
        score[i] <- zz$null.deviance - zz$deviance
        r <- fit$residuals
        w <- fit$weights
      }

      resdev <- c(resdev, fit$deviance)
      resdf <- c(resdf, fit$df.residual)
    }

    if (doscore) {
      zz <- eval(call(if (is.function(method)) "method" else method, x = x, y = r, weights = w, intercept = icpt))
      score[nvars] <- zz$null.deviance - zz$deviance
    }
  }

  resdf <- c(object$df.null, resdf, object$df.residual)
  resdev <- c(object$null.deviance, resdev, object$deviance)

  table <- data.frame(c(NA, -diff(resdf)), c(NA, pmax(0, -diff(resdev))), resdf, resdev)

  tl <- attr(object$terms, "term.labels")

  if (length(tl) == 0L) table <- table[1, , drop = FALSE]

  dimnames(table) <- list(c("NULL", tl), c("Df", "Deviance", "Resid. Df", "Resid. Dev"))

  if (doscore) table <- cbind(table, Rao = c(NA, score))

  title <- paste0("Analysis of Deviance Table", "\n\nModel: ", family, ", link: ", object$family$link,
                  "\n\nResponse: ", as.character(varlist[-1L])[1L], "\n\nTerms added sequentially (first to last)\n\n")
  df.dispersion <- Inf

  if (is.null(dispersion)) {
    dispersion <- summary(object, dispersion = dispersion)$dispersion
    df.dispersion <- if (dispersion == 1) Inf
    else object$df.residual
  }

  if (!is.null(test)) {
    if (test == "F" && df.dispersion == Inf) {
      fam <- object$family$family
      if (fam == "binomial" || fam == "poisson") warning(gettextf("using F test with a '%s' family is inappropriate", fam), domain = NA)
      else warning("using F test with a fixed dispersion is inappropriate")
    }

    table <- stats::stat.anova(table = table, test = test, scale = dispersion, df.scale = df.dispersion, n = NROW(x))
  }

  structure(table, heading = title, class = c("anova", "data.frame"))
}

#' Extract log-Likelihoods for \code{MCEMfit_lm} model objects
#'
#' Extract log-Likelihoods for \code{MCEMfit_lm} model objects. This function subtracts the entropy term from the observed likelihood.
#' @name logLik_MCEMfit_lm
#' @param object : fitted model objects of class \code{MCEMfit_lm}.
#' @param ... : further arguments passed through to \code{lm}.
#' @param REML : an optional logical value. If \code{TRUE} the restricted log-likelihood is returned, else, if \code{FALSE}, the log-likelihood is returned. Defaults to \code{FALSE}.
#' @return \code{logLik_MCEMfit_lm} produces output identical to \code{logLik}.
#' @importFrom stats logLik
#' @author Jakub Stoklosa, Wen-Han Hwang and David I. Warton.
#' @export
#' @seealso \code{\link{logLik}}
#'
logLik_MCEMfit_lm <- function (object, REML = FALSE, ...) {
  if (inherits(object, "mlm")) stop("'logLik_MCEMfit_lm' does not support multiple responses")
  res <- object$residuals
  p <- object$rank
  N <- length(res)
  if (is.null(w <- object$weights)) w <- rep.int(1, N)
  else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      N <- length(res)
      w <- w[!excl]
    }
  }

  N0 <- N
  if (REML) N <- N - p
  val <- 0.5*(sum(log(w)) - N*(log(2*pi) + 1 - log(N) + log(sum(w*res^2)))) - object$entropy
  if (REML) val <- val - sum(log(abs(diag(object$qr$qr)[1L:p])))
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  attr(val, "df") <- p + 1
  class(val) <- "logLik"
  val
}
