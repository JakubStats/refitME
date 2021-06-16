#-------------------------------------------------------------------------------------------------------
# Some functions used in the MCEM algorithm measurement error paper. These function are mainly required
# for the capture-recapture modelling example (Section 4.4) and the simulation studies (Section 5) of
# the manuscript.
#-------------------------------------------------------------------------------------------------------

# Load R-packages.

library(splines)
library(simex)
library(rootSolve)
library(numDeriv)
library(ggplot2)
library(nleqslv)
library(fGarch)
library(gamlss)
library(scales)

#------------------------------------------------------------------------------------------------------

# A logit function.

logit_fun <- function(xb) 1/(1 + exp(-xb))

# Some conditional score (CS) functions required for the example 4.

# Since we need to solve (non-linear) estimating equations we require differentiation solvers
# (e.g., Newton's method). We need to set up the functions and set some tolerance parameter like
# the error tols. and some termination conditions.

# Numerical partial differential solver. Used in newton.1.

df.1 <- function(f, par0, y, w1, tau, sigma.sq.u, delta = 0.0005) {
  k <- length(par0)

  d1 <- diag(rep(delta, k))

  df.1 <- matrix(0, k, k)

  for(i in 1:k) {
    df.1[, i] <- (f(par0 + d1[i, ], y, w1, tau, sigma.sq.u) - f(par0 - d1[i, ], y, w1, tau, sigma.sq.u))/(2*delta)
  }

  return(df.1)
}

# Estimating functions for the conditional score method.

est.cs <- function(betain, y, w1, tau, sigma.sq.u) {
  D1 <- length(y)
  YK <- 1:tau

  Delta <- w1 + y*betain[2]*sigma.sq.u

  eta <- betain[1] + betain[2]*Delta

  expt1 <- exp(outer(eta, YK, "*") - 0.5*betain[2]^2*sigma.sq.u*matrix(YK^2, D1, tau, byrow = T))
  expt2 <- matrix(choose(tau, YK), D1, tau, byrow = T)*expt1

  normc <- apply(expt2, 1, sum)
  PY.CD <- expt2/normc
  EY.CD <- apply(matrix(YK, D1, tau, byrow = T)*PY.CD, 1, sum, na.rm = T)

  score.0 <- cbind(rep(1, D1), Delta)*(c(y - EY.CD))

  score.1 <- c()

  for(kk in 1:D1) {
    scoreAA <- score.0[kk, ]

    if (is.matrix(scoreAA)) score.1a <- apply(scoreAA, 2, sum)
    if (!is.matrix(scoreAA)) score.1a <- scoreAA

    score.1 <- rbind(score.1, score.1a)
  }

  score <- apply(score.1, 2, sum)

  return(score)
}

# Nhat calculation for the conditional score method.

N.CS.est <- function(beta.est, y, w1, tau, sigma.sq.u) {
  D1 <- length(y)
  YK <- 1:tau

  Delta <- w1 + y*beta.est[2]*sigma.sq.u

  eta <- beta.est[1] + beta.est[2]*Delta

  expt1 <- exp(outer(eta, YK, "*") - 0.5*beta.est[2]^2*sigma.sq.u*matrix(YK^2, D1, tau, byrow = T))
  expt2 <- matrix(choose(tau, YK), D1, tau, byrow = T)*expt1

  normc <- apply(expt2, 1, sum)

  pbar1 <- 1 - (1/(1 + normc))

  N.hat <- sum(1/pbar1)

  return(N.hat)
}

est.cs2 <- function(betain, y, w1, tau, sigma.sq.u) {
  D1 <- length(y)
  YK <- 1:tau

  Delta <- w1 + y*betain[2]*sigma.sq.u

  eta <- betain[1] + betain[2]*Delta

  expt1 <- exp(outer(eta, YK, "*") - 0.5*betain[2]^2*sigma.sq.u*matrix(YK^2, D1, tau, byrow = T))
  expt2 <- matrix(choose(tau, YK), D1, tau, byrow = T)*expt1

  normc <- apply(expt2, 1, sum)
  PY.CD <- expt2/normc
  EY.CD <- apply(matrix(YK, D1, tau, byrow = T)*PY.CD, 1, sum, na.rm = T)

  score.0 <- cbind(rep(1, D1), Delta)*(c(y - EY.CD))

  score <- apply(score.0, 2, sum)

  return(score)
}

est.cs3 <- function(betain, y, w1, tau, sigma.sq.u) {
  D1 <- length(y)
  YK <- 1:tau

  Delta <- w1 + y*betain[2]*sigma.sq.u

  eta <- betain[1] + betain[2]*Delta

  expt1 <- exp(outer(eta, YK, "*") - 0.5*betain[2]^2*sigma.sq.u*matrix(YK^2, D1, tau, byrow = T))
  expt2 <- matrix(choose(tau, YK), D1, tau, byrow = T)*expt1

  normc <- apply(expt2, 1, sum)
  PY.CD <- expt2/normc
  EY.CD <- apply(matrix(YK, D1, tau, byrow = T)*PY.CD, 1, sum, na.rm = T)

  score.0 <- cbind(rep(1, D1), Delta)*(c(y - EY.CD))

  return(score.0)
}

# CS variance calculation function.

var.CS <- function(beta.est, y, w1, tau, sigma.sq.u) {
  D1 <- length(y)
  YK <- 1:tau

  Delta <- w1 + y*beta.est[2]*sigma.sq.u

  eta <- beta.est[1] + beta.est[2]*Delta

  expt1 <- exp(outer(eta, YK, "*") - 0.5*beta.est[2]^2*sigma.sq.u*matrix(YK^2, D1, tau, byrow = T))
  expt2 <- matrix(choose(tau, YK), D1, tau, byrow = T)*expt1

  normc <- apply(expt2, 1, sum)

  cprob <- c()

  for(kk in 1:D1) {
    cprob1 <- sum(1/normc[kk])
    cprob <- c(cprob, cprob1)
  }

  pbar1 <- 1/(1 + (1/y)*(cprob))

  ScA <- ScA <- t(est.cs3(beta.est, y, w1, tau, sigma.sq.u))%*%(est.cs3(beta.est, y, w1, tau, sigma.sq.u))
  vv0 <- nleqslv(beta.est, est.cs2, y = y, w1 = w1, tau = tau, sigma.sq.u = sigma.sq.u, jacobian = T, method = c("Newton"))$jac

  vv <- (solve(vv0))%*%ScA%*%t(solve(vv0))

  aa <- df.1(N.CS.est, beta.est, y, w1, tau, sigma.sq.u)
  bb <- matrix(aa[1, 1:2], 2, 1)
  vv.2 <- t(bb)%*%vv%*%bb

  var.N <- (sum((1 - pbar1)/(pbar1*pbar1)) + vv.2)

  var.p <- c(diag(vv), var.N)

  return(var.p)
}

# One (quadratic) covariate predict function (output is the RMSE).

predict_B_glm <- function(mod, datTest, sigma.sq.uTest, B, eta_true = NA) {
  n <- nrow(datTest)

  sigma.sq.u1 <- sigma.sq.uTest

  K1 <- length(mod$coef)

  U1_j <- rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1, B)))
  X1_j <- rep(datTest[, 2], B) - U1_j

  if (K1 == 2) X <- data.frame(cbind(rep(1, B*n), X1_j))
  if (K1 == 3) X <- data.frame(cbind(rep(1, B*n), X1_j, (X1_j)^2))

  colnames(X) <- names(coef(mod))

  etaPred <- as.matrix(X)%*%coef(mod)

  RMSE <- sqrt(mean((rep(eta_true, B) - etaPred)^2, na.rm = T))

  values <- list(RMSE = RMSE)

  return(values)
}

predict_B_gam <- function(mod, datTest, sigma.sq.uTest, B, eta_true = NA) {
  n <- nrow(datTest)

  sigma.sq.u1 <- sigma.sq.uTest

  U1_j <- rnorm(n*B, 0, sd = sqrt(rep(sigma.sq.u1, B)))
  X1_j <- rep(datTest[, 2], B) - U1_j

  datTest_new <- data.frame(cbind(X1_j))

  colnames(datTest_new) <- "x1"

  etaPred <- predict(mod, newdata = datTest_new, type = "link")

  RMSE <- sqrt(mean((rep(eta_true, B) - etaPred)^2, na.rm = T))

  values <- list(RMSE = RMSE)

  return(values)
}

#------------------------------------------------------------------------------------------------------

# Multiple plot function for ggplot.

# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects).
#
# - cols: Number of columns in layout.
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist.

  plots <- c(list(...), plotlist)
  numPlots <- length(plots)

  # If layout is NULL, then use 'cols' to determine layout.

  # Make the panel
  # ncol: Number of columns of plots
  # nrow: Number of rows needed, calculated from # of cols

  if (is.null(layout)) layout <- matrix(seq(1, cols*ceiling(numPlots/cols)), ncol = cols, nrow = ceiling(numPlots/cols))
  if (numPlots == 1) print(plots[[1]])

  # Set up the page.

  else{
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location.

    # Get the i, j matrix positions of the regions that contain this subplot.

    for(i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    }
  }
}

#------------------------------------------------------------------------------------------------------

# Simulation functions (Section 5).

# A coverage probability function.

cov.pc <- function(true, est, esd, N.sim) {
  p <- length(true)
  count <- rep(0, p)

  for(i in 1:p) {
    count[i] <- length(which(((est[, i] - 1.96*esd[, i]) < true[i] & true[i] < (est[, i] + 1.96*esd[, i])) == T))/N.sim
  }

  count
}

# A coverage probability and interval score function.

cov.pc_N <- function(true, est, esd, N.sim) {
  count1 <- length(which(((est - 1.96*esd) < true & true < (est + 1.96*esd)) == T))/N.sim

  list(cp = count1)
}

# Simulation function for a logistic regression ME and missing CR sim.

genCRdata <- function(N, beta, tau, sigma.sq.u, par.int) {
  x1 <- rnorm(N, 0, sd = 1)

  if (par.int == 2) X <- cbind(rep(1, N), x1)
  if (par.int == 3) X <- cbind(rep(1, N), x1, x1^2)
  if (par.int == "spline") X <- NA

  w1 <- x1 + rnorm(N, 0, sd = sqrt(sigma.sq.u))  # Add/contaminate the covariate with error.

  if (par.int == 2) {
    W <- cbind(rep(1, N), w1)
    colnames(W) <- c("Intercept", "w1")
  }
  if (par.int == 3) {
    W <- cbind(rep(1, N), w1, w1^2)
    colnames(W) <- c("Intercept", "w1", "I(w1^2)")
  }

  if (par.int == "spline") W <- NA

  if (par.int != "spline") eta <- X%*%beta
  if (par.int == "spline") eta <- sin(0.85*x1)
  if (par.int == "spline") eta <- cos(2*x1) + 2

  p <- 1/(1 + exp(-drop(eta)))

  y1 <- matrix(runif(N*tau) < p, N, tau)
  y <- apply(y1, 1, sum)
  y2 <- y1 + 0

  hist <- y2[y > 0, ]
  X <- X[y > 0, ]
  W <- W[y > 0, ]

  list(hist = hist, y = apply(hist, 1, sum), W = W, X = X)
}
