#------------------------------------------------------------------------------
# JSS_Replication_Code_Examples.R
#
# This file reproduces all the results (Tables 1-2, Figure 3-4 and 10-11) for
# Examples 1 to 4 of Section 4 and Appendix B of the manuscript.
#------------------------------------------------------------------------------

# Example 1: A simple GLM example taken from Carroll et al. (2006) - The Framingham heart study data set.

library(refitME)
library(simex)

set.seed(2020)

data(Framinghamdata)

glm_naiv <- glm(Y ~ w1 + z1 + z2 + z3, x = TRUE, family = binomial, data = Framinghamdata)

sigma.sq.u <- 0.006295 # Measurement error variance, obtained from Carroll et al. (2006) monograph.

rel.rat <- (1 - sigma.sq.u/var(Framinghamdata$w1))*100

rel.rat

B <- 100 # No. Monte Carlo replication values/SIMEX simulations.

start <- Sys.time()

glm_simex <- simex(glm_naiv, SIMEXvariable = c("w1"), measurement.error = cbind(sqrt(sigma.sq.u)), B = B) # SIMEX.

end <- Sys.time()
t1 <- difftime(end, start, units = "secs")
comp.time <- c(t1)

start <- Sys.time()

W <- as.matrix(Framinghamdata$w1) # Matrix of error-contaminated covariates.

glm_MCEM <- refitME(glm_naiv, sigma.sq.u, B)

end <- Sys.time()
t2 <- difftime(end, start, units = "secs")
comp.time <- c(comp.time, t2)

est.beta <- rbind(coef(glm_naiv), coef(glm_simex), coef(glm_MCEM))
est.beta.se <- rbind(sqrt(diag(vcov(glm_naiv))),
                     sqrt(diag(glm_simex$variance.jackknife)),
                     glm_MCEM$se)

# Parameter estimates.

row.names(est.beta) = row.names(est.beta.se) <- c("Naive GLM", "SIMEX", "MCEM")
colnames(est.beta) = colnames(est.beta.se) <- c("(Intercept)", "SBP", "chol. level", "age", "smoke")

round(est.beta, digits = 3)

# Standard errors.

round(est.beta.se, digits = 3)

# Computational times.

names(comp.time) <- c("SIMEX", "MCEM")
round(comp.time, digits = 3)

# Check a generic functions, e.g., summary, ANOVA, etc.

summary(glm_naiv)
summary(glm_MCEM)

anova(glm_naiv)
anova(glm_MCEM)

#----------------------------------------------------------------------------------------

# Example 2: A GAM example taken from Ganguli et al. (2005) - The air pollution data set.

library(refitME)
library(SemiPar)
library(mgcv)

set.seed(2020)

data(milan.mort)

dat.air <- milan.mort

Y <- dat.air[, 6]

n <- length(Y)

w1 <- log(dat.air[, 9])

z1 <- (dat.air[, 1])
z2 <- (dat.air[, 4])
z3 <- (dat.air[, 5])

dat <- data.frame(cbind(Y, z1, z2, z3, w1))

# Poisson models (we fitted negative binomial models but found no over-dispersion).

# Naive GAM.

gam_naiv <- gam(Y ~ s(w1) + s(z1, k = 25) + s(z2) + s(z3),
                  family = "poisson", data = dat)

# MCEM GAM.

sigma.sq.u <- 0.09154219 # Rel. ratio of 70%.

rel.rat <- (1 - sigma.sq.u/var(dat$w1))*100
rel.rat

gam_MCEM <- refitME(gam_naiv, sigma.sq.u)

# Plots (examine all smooth terms against covariates).

xlab.names <- c("log(TSP)", "Day", "Temp", "Humidity")

plot_gam_naiv <- plot(gam_naiv, select = 1)

op <- par(mfrow = c(2, 2), las = 1)

for(i in 1:4) {
  if (i == 1) {
    plot(gam_MCEM, select = i, ylim = c(-0.35, 0.1), xlim = range(plot_gam_naiv[[1]]$x), rug = FALSE, col = "blue", all.terms = TRUE,
         xlab = xlab.names[i], ylab = "s(Mortality counts)", lwd = 2, cex.lab = 1.3, cex.axis = 1.3,
         cex.main = 2, font.lab = 1.1, cex = 1.4, shade = T)
    lines(plot_gam_naiv[[1]]$x, plot_gam_naiv[[1]]$fit, type= "l", col = "red", lwd = 2, lty = 2)
    title(main = bquote("Reliability ratio of predictor is" ~ .(rel.rat) ~ "%"), outer = F, line = 1, cex = 1.4)
    legend("bottomright", c("Naive GAM", "MCEM GAM"), col = c("red", "blue"), lty = c(2, 1), lwd = 2, bty = "n")
    for(j in 1:2) {
      axis(j, labels = FALSE)
    }
  }
  if (i == 2) {
    plot(gam_MCEM, select = i, ylim = c(-0.25, 0.3), rug = FALSE, col = "blue", all.terms = TRUE,
         xlab = xlab.names[i], ylab = "s(Mortality counts)", lwd = 2, cex.lab = 1.3, cex.axis = 1.3,
         cex.main = 2, font.lab = 1.1, cex = 1.4, shade = T)
    lines(plot_gam_naiv[[2]]$x, plot_gam_naiv[[2]]$fit, type= "l", col = "red", lwd = 2, lty = 2)
    for(j in 1:2) {
      axis(j, labels = FALSE)
    }
  }
  if (i == 3) {
    plot(gam_MCEM, select = i, ylim = c(-0.2, 0.4), rug = FALSE, col = "blue", all.terms = TRUE,
         xlab = xlab.names[i], ylab = "s(Mortality counts)", lwd = 2, cex.lab = 1.3, cex.axis = 1.3,
         cex.main = 2, font.lab = 1.1, cex = 1.4, shade = T)
    lines(plot_gam_naiv[[3]]$x, plot_gam_naiv[[3]]$fit, type= "l", col = "red", lwd = 2, lty = 2)
    for(j in 1:2) {
      axis(j, labels = FALSE)
    }
  }
  if (i == 4) {
    plot(gam_MCEM, select = i, ylim = c(-0.06, 0.08), rug = FALSE, col = "blue", all.terms = TRUE,
         xlab = xlab.names[i], ylab = "s(Mortality counts)", lwd = 2, cex.lab = 1.1, cex.axis = 1.1,
         cex.main = 2, font.lab = 1.1, cex = 1.4, shade = T)
    lines(plot_gam_naiv[[4]]$x, plot_gam_naiv[[4]]$fit, type= "l", col = "red", lwd = 2, lty = 2)
    for(j in 1:2) {
      axis(j, labels = FALSE)
    }
  }
}

title(main = "MCEM (Poisson GAM) fitted to the air pollution data.", outer = T, line = -2)

par(op)

detach(package:mgcv)

#----------------------------------------------------------------------------------------

# Example 3: A Point Process model (PPM) example using Corymbia eximia presence-only data.

library(refitME)
library(caret)

set.seed(2020)

data(Corymbiaeximiadata)

dat <- Corymbiaeximiadata

coord.dat <- cbind(dat$X, dat$Y)
colnames(coord.dat) <- c("Longitude", "Latitude")

scen_par <- 1

Y1 <- dat$Y.obs

n <- length(Y1)

Rain <- dat$Rain
D.Main <- dat$D.Main
MNT <- dat$MNT

# PPM - using a Poisson GLM.

p.wt <- rep(1.e-6, length(Y1))
p.wt[Y1 == 0] <- 1

Y <- Y1/p.wt

X <- cbind(MNT, Rain, sqrt(D.Main))
colnames(X) <- c("w1", "z1", "z2")

# Training/test data.

train <- (1:n)
n.train <- n.test <- n
test <- (1:n)

data.tr <- data.frame(cbind(Y[train], X[train, ], p.wt[train]))
colnames(data.tr)[c(1, ncol(data.tr))] <- c("Y", "p.wt")

PPM_naiv <- glm(Y ~ poly(w1, degree = 2, raw = TRUE) + poly(z1, degree = 2, raw = TRUE) + poly(z2, degree = 2, raw = TRUE), family = "poisson", weights = p.wt, data = data.tr)

# PPM - using MCEM model.

# Set the measurement error variance here (we considered three values in this example).

#sigma.sq.u <- 0.25
#sigma.sq.u <- 0.5
sigma.sq.u <- 1

W <- MNT # Max. temp (measured with error).
W.train <- W[train]
W.test <- W[test]

rel.rat <- (1 - sigma.sq.u/var(W.train))*100
print(rel.rat)  # Express as a relative percentage (%).

PPM_MCEM <- refitME(PPM_naiv, sigma.sq.u)

# Prediction on test and training data.

X.f <- cbind(rep(1, length(MNT)), poly((MNT), degree = 2, raw = TRUE), poly(Rain, degree = 2, raw = TRUE), poly(sqrt(D.Main), degree = 2, raw = TRUE))
W.s <- MNT + scen_par
X.s <- cbind(rep(1, length(MNT)), poly(W.s, degree = 2, raw = TRUE), poly(Rain, degree = 2, raw = TRUE), poly(sqrt(D.Main), degree = 2, raw = TRUE))

X.f.test <- X.f[test, ]
X.s.test <- X.s[test, ]

preds1 <- exp(X.f.test%*%coef(PPM_naiv))
preds2 <- exp(X.s.test%*%coef(PPM_naiv))

preds3 <- exp(X.f.test%*%coef(PPM_MCEM))
preds4 <- exp(X.s.test%*%coef(PPM_MCEM))

coord.dat1 <- coord.dat[test, ]

preds1a <- preds1
preds2a <- preds2
preds3a <- preds3
preds4a <- preds4

p1 <- cbind(coord.dat1[, 1], coord.dat1[, 2], preds3a, rep(1, length(preds3)))
p2 <- cbind(coord.dat1[, 1], coord.dat1[, 2], preds4a, rep(2, length(preds4)))
p3 <- cbind(coord.dat1[, 1], coord.dat1[, 2], preds1a, rep(3, length(preds1)))
p4 <- cbind(coord.dat1[, 1], coord.dat1[, 2], preds2a, rep(4, length(preds2)))
pred.dats <- rbind(p1, p2, p3, p4)

op <- par(mfrow = c(2, 2), las = 1)

plots <- factor(pred.dats[, 4], labels = c("(c) PPM MCEM", "(d) PPM MCEM (+1 degrees)", "(a) PPM", "(b) PPM (+1 degrees)"))
pred.dats <- as.data.frame(pred.dats)
colnames(pred.dats) <- c("x", "y", "preds", "plot")
levelplot(preds ~ x + y | plots, cex = 1, data = pred.dats, asp = "iso", ylab = "Latitude", xlab = "Longitude", col.regions = heat.colors(1024)[900:1], cuts = 900, main = list("", cex = 5), scales = list(y = list(draw = FALSE), x = list(draw = FALSE)), colorkey = list(labels = list(cex = 0.8)))

par(op)

weighted.mean(coord.dat1[, 2], w = preds1)
weighted.mean(coord.dat1[, 2], w = preds2)
weighted.mean(coord.dat1[, 2], w = preds3)
weighted.mean(coord.dat1[, 2], w = preds4)

#-----------------------------------------------------------------------------

# Example 4: A capture-recapture example using the Hong Kong prinia bird data.

library(refitME)
library(VGAM)

setwd(system.file(package = 'refitME'))
source('MCEM_prog.R')

set.seed(2020)

data(Priniadata)

tau <- 17

# Standardize the bird wing length predictor.

w2 <- c(scale(Priniadata$w1))
Priniadata$w2 <- w2

y <- Priniadata$cap

# Start fitting all models here (model M_h).

# Naive model.

CR_naiv1 <- vglm(cbind(cap, noncap) ~ w2, VGAM::posbinomial(omit.constant = TRUE, parallel = TRUE ~ w2), data = Priniadata, trace = FALSE)

# Conditional score (CS) model.

sigma.sq.u <- 0.37/var(Priniadata$w1) # ME variance.

CS_beta.est <- nleqslv(coef(CR_naiv1), est.cs, y = y, w1 = w2, tau = tau, sigma.sq.u = sigma.sq.u, method = c("Newton"))$x
CS_N.est <- N.CS.est(CS_beta.est, y, w2, tau, sigma.sq.u)

var.ests1 <- var.CS(CS_beta.est, y, w2, tau, sigma.sq.u)
CS_beta.est.se <- sqrt(var.ests1)[1:2]

# MCEM model.

CR_naiv2 <- vgam(cbind(cap, noncap) ~ s(w2, df = 2), VGAM::posbinomial(omit.constant = TRUE, parallel = TRUE ~ s(w2, df = 2)), data = Priniadata, trace = FALSE)

B <- 100

CR_MCEM <- refitME(CR_naiv2, sigma.sq.u, B)

N.hat_MCEM <- CR_MCEM$N.est

N.hat.se_MCEM <- CR_MCEM$N.est.se

# Model selection using AIC.

AIC.mod <- c(round(AIC(CR_naiv1), digits = 2), NA, round(AIC(CR_naiv2), digits = 2), round(CR_MCEM$aic.value, digits = 2))

# Combine all results.

Nhat.mod <- c(round(CR_naiv1@extra$N.hat, digits = 2), round(CS_N.est, digits = 2), round(CR_naiv2@extra$N.hat, digits = 2), round(N.hat_MCEM, digits = 2))

Nhat.se.mod <- c(round(CR_naiv1@extra$SE.N.hat, digits = 2), round(sqrt(var.ests1)[3], digits = 2), round(CR_naiv2@extra$SE.N.hat, digits = 2), round(N.hat.se_MCEM, digits = 2))

est <- cbind(AIC.mod, Nhat.mod, Nhat.se.mod)

row.names(est) <- c("Naive GLM", "Conditional Score", "Naive GAM", "MCEM")
colnames(est) <- c("AIC", "N hat", "SE(N hat)")

est

detach(package:VGAM)
