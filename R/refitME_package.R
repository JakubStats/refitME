########################################################################
## refitME_package.r
##
## Source file for the "refitME" R-package.
########################################################################

## Required R-packages.

library(mvtnorm);
library(MASS);
library(SDMTools);
library(mgcv);
library(Matrix);
library(sandwich);

#' MCEMfit_glm
#'
#' Function for wrapping the MCEM algorithm on GLMs where covariates are subject to measurement error/error-in-varaibles.
#' @name MCEMfit_glm
#' @param mod : a glm object (this is the naive fitted model). Make sure the first stored variable is the contaminated covariate (W).
#' @param family : a specified family/distribution.
#' @param sigma.sq.u : measurement error variance.
#' @param sigma.sq.e : variance of the true covariate (X).
#' @param B : the number of Monte Carlo replication values.
#' @param epsilon : a set convergence threshold.
#' @param theta.est : an initial value for the dispersion parameter (this is required for fitting negative binomial models).
#' @param shape.est : an initial value for the shape parameter (this is required for fitting gamma models).
#' @return \code{MCEMfit_glm} returns model coef estimates with standard errors.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J. and Warton, D.I. (2019). A general algorithm for error-in-variables using Monte Carlo expectation maximization.
#' @export
#' @seealso \code{\link{MCEMfit_gam}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown tutorial with examples.
MCEMfit_glm<-function(mod,family,sigma.sq.u,sigma.sq.e,B,epsilon=0.00001,theta.est=1,shape.est=1)
  {
  options(warn=-1);

  reps<-0;
  cond<-TRUE;

  Y<-mod$model[,1];
  W<-mod$model[,-1];
  w1<-mod$model[,2];
  n<-length(Y);

  mod.terms<-attr(mod$terms,"term.labels");

  d<-length(mod.terms);

  if(d==1){p1<-1;}

  if(d>1)
    {
    if(sum((w1)^2!=(W[,2]))==n){p1<-1;}
    if(sum((w1)^2==(W[,2]))==n){p1<-2;}
    }

  muPred<-rep(stats::predict(mod,type="response"),B);

  beta.est<-stats::coef(mod);

  sigma.sq.u1<-sigma.sq.u;
  sigma.sq.e1<-sigma.sq.e;

  U1_j<-stats::rnorm(n*B,0,sd=sqrt(rep(sigma.sq.u1,B)));
  X1_j<-rep(w1,B)-U1_j;

  X<-cbind(rep(1,B*n),X1_j);
  if(p1==2){X<-cbind(rep(1,B*n),X1_j,(X1_j)^2);}

  if(!is.null(ncol(W)))
    {
    if(d>2 & p1==2){X<-cbind(X,do.call(rbind,replicate(B,W[,-c(1:p1)],simplify=FALSE)));}
    if(d>2 & p1==1){X<-cbind(X,do.call(rbind,replicate(B,W[,-c(1:p1)],simplify=FALSE)));}
    if(d==2 & p1==2){X<-X;}
    }

  X<-as.matrix(X);

  colnames(X)<-names(stats::coef(mod));

  colnames(X)[2]<-"x1";

  if(p1==2){colnames(X)[2:3]<-c("x1","I(x1^2)");}

  mu.e1<-mean(X1_j);

  bigY<-rep(Y,B);

  while(cond)
    {
    ## MC and E-step.

    prX<-stats::dnorm(X1_j,mu.e1,sd=sqrt(sigma.sq.e1));

    if(family=="gaussian"){prY<-stats::dnorm(bigY,muPred,1);}
    if(family=="binomial"){prY<-stats::dbinom(bigY,1,muPred);}
    if(family=="poisson"){prY<-stats::dpois(bigY,muPred);}
    #if(family=="Gamma"){prY<-stats::dgamma(bigY,shape=shape.est,rate=shape.est/muPred);}
    if(family=="Gamma"){prY<-stats::dgamma(bigY,shape=shape.est,scale=muPred/shape.est);}
    if(family=="negbin"){prY<-stats::dnbinom(bigY,size=theta.est,mu=muPred);}

    ## M-step (updates).

    bigW<-matrix(prY*prX,n,B);

    sumW<-rep(apply(bigW,1,sum,na.rm=T),B);

    weights1<-as.vector(bigW)/sumW;
    weights1[is.nan(weights1)]<-0;

    if(family=="gaussian")
      {
      mod<-stats::lm(bigY~X-1,weights=weights1);
      sigma.sq.est<-(summary(mod)$sigma)^2;
      }

    if(family=="binomial"){mod<-stats::glm(bigY~X-1,weights=weights1,family="binomial");}
    if(family=="poisson"){mod<-stats::glm(bigY~X-1,weights=weights1,family="poisson");}
    if(family=="Gamma"){mod<-stats::glm(bigY~X-1,weights=weights1,family=Gamma(link="log"));}
    if(family=="negbin"){mod<-MASS::glm.nb(bigY~X-1,weights=weights1,init.theta=theta.est);}

    beta.update<-stats::coef(mod);
    if(family=="negbin"){theta.update<-mod$theta;}
    muPred<-stats::predict(mod,type="response");
    #if(family=="Gamma"){shape.update<-summary(mod)[14]$dispersion;}
    if(family=="Gamma")
      {
      #shape.update<-MASS::gamma.shape(mod,it.lim=10000,eps.max=0.0000001,verbose =T)$alpha;
      shape.update<-summary(mod)[14]$dispersion
      muPred<-stats::predict(mod,type="response",dispersion=shape.update);
      }

    sigma.sq.e1.update<-SDMTools::wt.var(X[,2],w=weights1);

    mu.e1.update<-stats::weighted.mean(X[,2],w=weights1);

    ## Convergence monitoring.

    beta.norm<-sum((beta.est-beta.update)^2);
    if(family=="negbin"){theta.norm<-sum((theta.est-theta.update)^2);}
    if(family=="Gamma"){shape.norm<-sum((shape.est-shape.update)^2);}
    #if(family=="Gamma"){shape.norm<-abs(shape.est-shape.update);}

    diff.sig_e<-abs(sigma.sq.e1.update-sigma.sq.e1);
    diff.mu_e<-sum((mu.e1.update-mu.e1)^2);

    reps<-reps+1;   # Keeps track of number of iterations.

    if(family=="binomial" | family=="poisson" | family=="gaussian")
      {
      if(diff.mu_e<epsilon && diff.sig_e<epsilon && beta.norm<epsilon)
        {
        cond<-FALSE;
        print("convergence :-)");
        print(reps);
        break;
        }
      }

    if(family=="negbin")
      {
      if(diff.mu_e<epsilon && diff.sig_e<epsilon && beta.norm<epsilon && theta.norm<epsilon)
        {
        cond<-FALSE;
        print("convergence :-)");
        print(reps);
        break;
        }
      }

    if(family=="Gamma")
      {
      if(diff.mu_e<epsilon && diff.sig_e<epsilon && beta.norm<epsilon && shape.norm<epsilon)
        {
        cond<-FALSE;
        print("convergence :-)");
        print(reps);
        break;
        }
      }

    ## Update parameters.

    beta.est<-beta.update;
    if(family=="negbin"){theta.est<-theta.update;}
    if(family=="Gamma"){shape.est<-shape.update;}
    sigma.sq.e1<-sigma.sq.e1.update;
    mu.e1<-mu.e1.update;
    }

  eff.samp.size<-1/sum(weights1^2);

  ## Standard error calculations start here.

  beta.est.se1<-sqrt(diag(stats::vcov(mod))); # Naive SE estimator.

  #if(family=="gaussian"){beta.est.se1<-sqrt(diag(stats::vcov(mod)/sigma.sq.est));}
  if(family=="gaussian"){beta.est.se1<-sqrt(diag(stats::vcov(mod)*B));}

  estfun_mat<-sandwich::estfun(mod);
  #if(family=="gaussian"){estfun_mat<-sandwich::estfun(mod)/sigma.sq.est;}
  if(family=="gaussian"){estfun_mat<-sandwich::estfun(mod)*B;}

  K1<-length(beta.update);

  S_1<-matrix(0,nrow=K1,ncol=K1);
  SS_1<-matrix(0,nrow=K1,ncol=K1);

  ind_mat<-matrix(1:(n*B),ncol=n,byrow=T);

  iii1<-1;

  for(iii in 1:n)
    {
    index_vec<-ind_mat[,iii];

    S_1<-S_1+(apply(estfun_mat[index_vec,],2,sum))%*%t(apply(estfun_mat[index_vec,],2,sum));
    }

  if(family=="gaussian")
    {
    #SS_1<-t(sandwich::estfun(mod)/sigma.sq.est)%*%(sandwich::estfun(mod)/sigma.sq.est/mod$weights);
    #u.bar<-solve(stats::vcov(mod)/sigma.sq.est);
    #beta.est.se2<-sqrt(diag(solve(u.bar-SS_1*sigma.sq.est^2+S_1*sigma.sq.est^2)));
    SS_1<-t(sandwich::estfun(mod)*B)%*%(sandwich::estfun(mod)*B/mod$weights);
    u.bar<-solve(stats::vcov(mod)*B);
    beta.est.se2<-sqrt(diag(solve(u.bar-SS_1/B^2+S_1/B^2)));
    }
  if(family!="gaussian")
    {
    SS_1<-t(sandwich::estfun(mod))%*%(sandwich::estfun(mod)/mod$prior);
    u.bar<-solve(stats::vcov(mod));
    beta.est.se2<-sqrt(diag(solve(u.bar-SS_1+S_1)));
    }

  if(length(which(is.nan(beta.est.se2)))>0){beta.est.se2<-c(rep(NA,K1));}

  values<-list(beta=beta.est,beta.se1=beta.est.se1,beta.se2=beta.est.se2,mod=mod,eff.samp.size=eff.samp.size);

  return(values)
  }

#' MCEMfit_gam
#'
#' Function for wrapping the MCEM algorithm on GAMs where covariates are subject to measurement error/error-in-varaibles.
#' @name MCEMfit_gam
#' @param mod : a gam object (this is the naive fitted model). Make sure the first stored variable is the contaminated covariate (W).
#' @param family : a specified family/distribution.
#' @param sigma.sq.u : measurement error variance.
#' @param sigma.sq.e : variance of the true covariate.
#' @param B : the number of Monte Carlo replication values.
#' @param epsilon : a set convergence threshold.
#' @param theta.est : an initial value for the dispersion parameter (this is required for fitting negative binomial models).
#' @param shape.est : an initial value for the shape parameter (this is required for fitting gamma models).
#' @return \code{MCEMfit_glm} returns model coef estimates with standard errors.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J. and Warton, D.I. (2019). A general algorithm for error-in-variables using Monte Carlo expectation maximization.
#' @export
#' @seealso \code{\link{MCEMfit_glm}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown tutorial with examples.
MCEMfit_gam<-function(mod,family,sigma.sq.u,sigma.sq.e,B,epsilon=0.00001,theta.est=1,shape.est=10)
  {
  options(warn=-1);

  reps<-0;
  cond<-TRUE;

  Y<-mod$model[,1];
  W<-mod$model[,-1];
  w1<-mod$model[,2];
  n<-length(Y);

  mod.terms<-attr(mod$terms,"term.labels");

  d<-length(mod.terms);

  form.name<-stats::formula(mod);

  beta.est<-stats::coef(mod);

  muPred<-rep(stats::predict(mod,type="response"),B);

  sigma.sq.u1<-sigma.sq.u;
  sigma.sq.e1<-sigma.sq.e;

  U1_j<-stats::rnorm(n*B,0,sd=sqrt(rep(sigma.sq.u1,B)));
  X1_j<-rep(w1,B)-U1_j;

  mu.e1<-mean(X1_j);

  x1<-cbind(X1_j);

  bigY<-rep(Y,B);

  dat_new<-cbind(bigY,x1);

  form.name<-stats::as.formula(stats::update(stats::formula(mod),bigY~.+s(x1,k=5)-s(w1,k=5)));

  if(d>1)
    {
    dat_new<-cbind(bigY,x1);

    for(jj in 2:d)
      {
      dat_new<-cbind(dat_new,rep(W[,jj],B));
      }
    }

  colnames(dat_new)<-c("bigY",mod.terms);
  colnames(dat_new)[2]<-"x1";

  dat_new<-as.data.frame(dat_new);

  while(cond)
    {
    ## MC and E-step.

    prX<-stats::dnorm(X1_j,mu.e1,sigma.sq.e1);

    if(family=="gaussian"){prY<-stats::dnorm(bigY,muPred,1);}
    if(family=="binomial"){prY<-stats::dbinom(bigY,1,muPred);}
    if(family=="poisson"){prY<-stats::dpois(bigY,muPred);}
    if(family=="Gamma"){prY<-stats::dgamma(bigY,rate=shape.est/muPred,shape=shape.est);}
    if(family=="negbin"){prY<-stats::dnbinom(bigY,size=theta.est,mu=muPred);}

    ## M-step (updates).

    bigW<-matrix(prY*prX,n,B);

    sumW<-rep(apply(bigW,1,sum),B);

    #dat_new$bigW<-as.vector(bigW);
    #dat_new$sumW<-sumW;

    weights1<-as.vector(bigW)/sumW;
    weights1[is.nan(weights1)]<-0;

    dat_new$weights1<-weights1;

    if(family=="gaussian")
      {
      mod<-mgcv::gam(formula=form.name,weights=weights1,data=dat_new);
      sigma.sq.est<-(summary(mod)$dispersion);
      }
    if(family=="binomial"){mod<-mgcv::gam(formula=form.name,weights=weights1,family="binomial",gamma=1.4,data=dat_new);}
    if(family=="poisson"){mod<-mgcv::gam(formula=form.name,weights=weights1,family="poisson",data=dat_new);}
    if(family=="Gamma"){mod<-mgcv::gam(formula=form.name,weights=weights1,family=Gamma(link="log"),data=dat_new);}
    if(family=="negbin"){mod<-mgcv::gam(form.name,family=nb(),weights=weights1,data=dat_new);}

    beta.update<-stats::coef(mod);
    if(family=="negbin"){theta.update<-mod$family$getTheta(TRUE);}
    if(family=="Gamma"){shape.update<-MASS::gamma.shape(mod)[1]$alpha;}
    muPred<-stats::predict(mod,type="response");

    sigma.sq.e1.update<-SDMTools::wt.var(X1_j,w=weights1);

    mu.e1.update<-stats::weighted.mean(x1,w=weights1);

    ## Convergence monitoring.

    #gcv.norm<-sum((gcv.est-gcv.update)^2);
    beta.norm<-sum((beta.est-beta.update)^2);
    if(family=="negbin"){theta.norm<-sum((theta.est-theta.update)^2);}
    if(family=="Gamma"){shape.norm<-sum((shape.est-shape.update)^2);}
    diff.sig_e<-abs(sigma.sq.e1.update-sigma.sq.e1);
    diff.mu_e<-sum((mu.e1.update-mu.e1)^2);

    reps<-reps+1;   # Keeps track of number of iterations.

    if(family=="binomial" | family=="poisson" | family=="gaussian")
      {
      if((diff.mu_e<epsilon && diff.sig_e<epsilon && beta.norm<epsilon) | reps>50)
        {
        cond<-FALSE;
        print("convergence :-)");
        print(reps);
        break;
        }
      }

    if(family=="negbin")
      {
      if((diff.mu_e<epsilon && diff.sig_e<epsilon && beta.norm<epsilon && theta.norm<epsilon) | reps>50)
        {
        cond<-FALSE;
        print("convergence :-)");
        print(reps);
        break;
        }
      }

    if(family=="Gamma")
      {
      if((diff.mu_e<epsilon && diff.sig_e<epsilon && beta.norm<epsilon && shape.norm<epsilon) | reps>50)
        {
        cond<-FALSE;
        print("convergence :-)");
        print(reps);
        break;
        }
      }

    ## Update parameters.

    #gcv.est<-gcv.update;
    beta.est<-beta.update;
    if(family=="negbin"){theta.est<-theta.update;}
    if(family=="Gamma"){shape.est<-shape.update;}
    sigma.sq.e1<-sigma.sq.e1.update;
    mu.e1<-mu.e1.update;
    }

  eff.samp.size<-1/sum(weights1^2);

  ## Standard error calculations start here.

  beta.est<-stats::coef(mod);
  beta.est.se1<-sqrt(diag(stats::vcov(mod))); # Naive SE estimator.
  if(family=="gaussian"){beta.est.se1<-sqrt(diag(stats::vcov(mod)*B));}

  #if(family=="gaussian"){beta.est.se1<-sqrt(diag(stats::vcov(mod)/sigma.sq.est));}
  if(family=="gaussian"){beta.est.se1<-sqrt(diag(stats::vcov(mod)*B));}

  estfun_mat<-sandwich::estfun(mod);
  #if(family=="gaussian"){estfun_mat<-sandwich::estfun(mod)/sigma.sq.est;}
  if(family=="gaussian"){estfun_mat<-sandwich::estfun(mod)*B;}

  X<-stats::model.matrix(mod);
  K1<-ncol(X);

  ind_mat<-matrix(1:(n*B),ncol=n,byrow=T);

  S_1<-matrix(0,nrow=K1,ncol=K1);
  SS_1<-matrix(0,nrow=K1,ncol=K1);

  iii1<-1;

  for(iii in 1:n)
    {
    index_vec<-ind_mat[,iii];

    S_1<-S_1+(apply(estfun_mat[index_vec,],2,sum))%*%t(apply(estfun_mat[index_vec,],2,sum));
    }

  if(family=="gaussian")
    {
    u.bar<-solve(stats::vcov(mod)*B);
    beta.est.se2<-sqrt(diag(solve(u.bar-SS_1/B^2+S_1/B^2)));
    }
  if(family=="binomial" | family=="poisson" | family=="Gamma" | family=="negbin")
    {
    SS_1<-t(sandwich::estfun(mod))%*%(sandwich::estfun(mod)/mod$prior);
    u.bar<-solve(stats::vcov(mod));
    beta.est.se2<-sqrt(diag(solve(u.bar-SS_1+S_1)));
    }

  if(length(which(is.nan(beta.est.se2)))>0){beta.est.se2<-c(rep(NA,K1));}

  values<-list(beta=beta.est,beta.se1=beta.est.se1,beta.se2=beta.est.se2,mod=mod,eff.samp.size=eff.samp.size);

  return(values);
  }

#' refitME
#'
#' Function that extracts the fitted (naive) model object and wraps the MCEM algorithm to correct for measurement error/error-in-varaibles (currently available for lm, glm and gam, excludes quassi models).
#' @name refitME
#' @param mod : a gam object (this is the naive fitted model). Make sure the first stored variable is the contaminated one (W).
#' @param sigma.sq.u : measurement error variance.
#' @param B : the number of Monte Carlo replication values.
#' @param epsilon : convergence threshold.
#' @return \code{refitME} returns model coef estimates with standard errors.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J. and Warton, D.I. (2019). A general algorithm for error-in-variables using Monte Carlo expectation maximization.
#' @export
#' @seealso \code{\link{MCEMfit_glm}} and \code{\link{MCEMfit_gam}}
#' @source See \url{https://github.com/JakubStats/refitME} for an RMarkdown tutorial with examples.
refitME<-function(mod,sigma.sq.u,B,epsilon=0.00001)
  {
  sigma.sq.e<-stats::var(mod$model[,2])-sigma.sq.u;

  ob.type<-attr(mod,"class")[1];

  if(ob.type=="lm" | ob.type=="glm" | ob.type=="negbin")
    {
    if(ob.type=="lm"){family<-"gaussian";}
    if(ob.type=="glm"){family<-mod$family$family;}
    if(ob.type=="negbin"){family<-"negbin";}

    return(MCEMfit_glm(mod,family,sigma.sq.u,sigma.sq.e,B,epsilon))
    }
  if(ob.type=="gam")
    {
    family<-mod$family$family;

    return(MCEMfit_gam(mod,family,sigma.sq.u,sigma.sq.e,B,epsilon))
    }
  }



