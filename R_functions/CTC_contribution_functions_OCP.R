## Causal tail coefficient estimators and causality permutation Monte Carlo test.
## Code for the methodology described in the article:
## Olivier C. Pasche, Valerie Chavez-Demoulin, and Anthony C. Davison. 2021.
## Causal Modelling of Heavy-Tailed Variables and Confounders with Application to River Flow
## https://arxiv.org/abs/2110.06686
## Created by Olivier Colin PASCHE, EPFL Lausanne (CH), Spring and Autumn 2020.

# libraries
library("evd")
library("ismev")
library(causalXtreme) # code by Gnecco et al. (2019)

library(maxLik)

library(boot)

# Script Dependencies
source("../R_functions/utils.R")
source("../R_functions/Modified_gpd_fitting.R")
source("../R_functions/Bootstrap_helpers.R")
suppressWarnings(try(source("../R_functions/deprecated/Partial_CIs_diagnostics.R"), silent=TRUE))# Deprecated

## Peaks Over Threshold Generalized Pareto Distribution based Causal Tail Coefficient estimators

GPD_causal_tail_coeff <- function(X1, X2, k = NULL, threshold_q = 0.95, parametric_F1=TRUE, both_tails = FALSE){
  ## GPD causal tail coefficient estimator
  n <- NROW(X1) # number of observations
  if(is.null(k)){
    k<-floor(n ^ 0.4)
  }
  #Checks
  if (threshold_q>=1 | threshold_q<=0) {
    stop("threshold_q must be between 0 and 1.")
  }
  if (k <= 1 | k>floor((1-threshold_q)*n) ) {
    stop("k must be greater than 1 and smaller than the number of threshold exceedences.")
  }
  if (both_tails) {
    stop("both_tails not yet possible.")
  }
  
  r1 <- rank(X1, ties.method = "first") # ranks of X1
  r2 <- rank(X2, ties.method = "first") # ranks of X2
  
  #Fit POT GPD model and compute the fitted cdf in the data points for X2
  u2 <- quantile(X2,threshold_q)
  fit2<-fpot(X2,threshold=u2)
  FX2 <- pgpd(X2, loc=u2, scale=fit2$estimate[1], shape=fit2$estimate[2])
  FX2 <- (FX2*(1-threshold_q)+threshold_q)*(X2>=u2) + r2/n * (X2<u2)
  
  #Compute the parametric causal tail coefficient
  if(parametric_F1){
    #Fit POT GPD model and compute the fitted cdf in the data points for X1
    u1 <- quantile(X1,threshold_q)
    fit1<-fpot(X1,threshold=u1)
    FX1 <- pgpd(X1, loc=u1, scale=fit1$estimate[1], shape=fit1$estimate[2])
    FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
    #Compute the parametric causal tail coefficient
    #not always exactly k excesses for parametric approach
    ctc <- mean(FX2[FX1 > (1-k/n)])
  }else{
    ctc <- 1/k * sum(FX2[r1 > n - k])
  }
  return(ctc)
}

LGPD_causal_tail_coeff <- function(X1, X2, H, k = NULL, threshold_q = 0.95, parametric_F1=TRUE, return_MLEs=FALSE, both_tails = FALSE){
  ## Post-fit corrected H-conditional LGPD causal tail coefficient estimator
  n <- NROW(X1) # number of observations
  if(is.null(k)){
    k<-floor(n ^ 0.4)
  }
  #Checks
  if (threshold_q>=1 | threshold_q<=0) {
    stop("threshold_q must be between 0 and 1.")
  }
  if (k <= 1 | k>floor((1-threshold_q)*n) ) {
    stop("k must be greater than 1 and smaller than the number of threshold exceedences.")
  }
  if (both_tails) {
    stop("both_tails not yet possible.")
  }
  
  Hcs <- scale(H) #Centering and scaling H
  
  r1 <- rank(X1, ties.method = "first") # ranks of X1
  r2 <- rank(X2, ties.method = "first") # ranks of X2
  
  #Fit POT GPD model and compute the fitted cdf in the data points for X2
  u2 <- quantile(X2,threshold_q)
  fit2 <- gpd.fit(X2,threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE)
  FX2 <- pgpd(X2, loc=u2, scale=pmax(fit2$mle[1]+fit2$mle[2]*Hcs,0.0001)[1:length(Hcs)], shape=fit2$mle[3])
  FX2 <- (FX2*(1-threshold_q)+threshold_q)*(X2>=u2) + r2/n * (X2<u2)
  
  #Compute the parametric causal tail coefficient
  if(parametric_F1){
    #Fit POT GPD model and compute the fitted cdf in the data points for X1
    u1 <- quantile(X1,threshold_q)
    fit1 <- gpd.fit(X1,threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE)
    FX1 <- pgpd(X1, loc=u1, scale=pmax(fit1$mle[1]+fit1$mle[2]*Hcs,0.0001)[1:length(Hcs)], shape=fit1$mle[3])
    FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
    #Compute the parametric causal tail coefficient
    ctc <- mean(FX2[FX1 > (1-k/n)]) #not always exactly k excesses for parametric approach
    if(return_MLEs){
      return(list(ctc=ctc, F2_scale0=fit2$mle[1], F2_scale1=fit2$mle[2], F2_shape=fit2$mle[3],
                  se_F2_scale0=fit2$se[1], se_F2_scale1=fit2$se[2], se_F2_shape=fit2$se[3],
                  F1_scale0=fit1$mle[1], F1_scale1=fit1$mle[2], F1_shape=fit1$mle[3],
                  se_F1_scale0=fit1$se[1], se_F1_scale1=fit1$se[2], se_F1_shape=fit1$se[3]))
    }
  }else{
    ctc <- 1/k * sum(FX2[r1 > n - k])
    if(return_MLEs){
      return(list(ctc=ctc, F2_scale0=fit2$mle[1], F2_scale1=fit2$mle[2], F2_shape=fit2$mle[3],
                  se_F2_scale0=fit2$se[1], se_F2_scale1=fit2$se[2], se_F2_shape=fit2$se[3],
                  F1_scale0=NULL, F1_scale1=NULL, F1_shape=NULL,
                  se_F1_scale0=NULL, se_F1_scale1=NULL, se_F1_shape=NULL))
    }
  }
  return(ctc)
}

constrained_LGPD_CTC <- function(X1, X2, H, k = NULL, threshold_q = 0.95, parametric_F1=TRUE, return_MLEs=FALSE, both_tails = FALSE){
  ## Constrained fit H-conditional LGPD causal tail coefficient estimator
  n <- NROW(X1) # number of observations
  if(is.null(k)){
    k<-floor(n ^ 0.4)
  }
  #Checks
  if (threshold_q>=1 | threshold_q<=0) {
    stop("threshold_q must be between 0 and 1.")
  }
  if (k <= 1 | k>floor((1-threshold_q)*n) ) {
    stop("k must be greater than 1 and smaller than the number of threshold exceedences.")
  }
  if (both_tails) {
    stop("both_tails not yet possible.")
  }
  
  Hcs <- scale(H) #Centering and scaling H
  
  r1 <- rank(X1, ties.method = "first") # ranks of X1
  r2 <- rank(X2, ties.method = "first") # ranks of X2
  
  #Fit POT GPD model and compute the fitted cdf in the data points for X2
  u2 <- quantile(X2,threshold_q)
  fit2 <- gpd_constraints_fit_maxLik(X2,threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE,
                                     constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
  FX2 <- pgpd(X2, loc=u2, scale=fit2$mle[1]+fit2$mle[2]*Hcs, shape=fit2$mle[3])
  FX2 <- (FX2*(1-threshold_q)+threshold_q)*(X2>=u2) + r2/n * (X2<u2)
  
  #Compute the parametric causal tail coefficient
  if(parametric_F1){
    #Fit POT GPD model and compute the fitted cdf in the data points for X1
    u1 <- quantile(X1,threshold_q)
    fit1 <- gpd_constraints_fit_maxLik(X1,threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE,
                                       constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
    FX1 <- pgpd(X1, loc=u1, scale=fit1$mle[1]+fit1$mle[2]*Hcs, shape=fit1$mle[3])
    FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
    #Compute the parametric causal tail coefficient
    ctc <- mean(FX2[FX1 > (1-k/n)])
    if(return_MLEs){
      return(list(ctc=ctc, F2_scale0=fit2$mle[1], F2_scale1=fit2$mle[2], F2_shape=fit2$mle[3],
                  se_F2_scale0=fit2$se[1], se_F2_scale1=fit2$se[2], se_F2_shape=fit2$se[3],
                  F1_scale0=fit1$mle[1], F1_scale1=fit1$mle[2], F1_shape=fit1$mle[3],
                  se_F1_scale0=fit1$se[1], se_F1_scale1=fit1$se[2], se_F1_shape=fit1$se[3]))
    }
  }else{
    ctc <- 1/k * sum(FX2[r1 > n - k])
    if(return_MLEs){
      return(list(ctc=ctc, F2_scale0=fit2$mle[1], F2_scale1=fit2$mle[2], F2_shape=fit2$mle[3],
                  se_F2_scale0=fit2$se[1], se_F2_scale1=fit2$se[2], se_F2_shape=fit2$se[3],
                  F1_scale0=NULL, F1_scale1=NULL, F1_shape=NULL,
                  se_F1_scale0=NULL, se_F1_scale1=NULL, se_F1_shape=NULL))
    }
  }
  return(ctc)
}

student_constrained_LGPD_CTC <- function(X1, X2, H, k = NULL, threshold_q = 0.95, t_df=t_df, parametric_F1=TRUE, return_MLEs=FALSE, both_tails = FALSE){
  ## Constrained fit H-conditional LGPD causal tail coefficient estimator for known Student t data distribution
  n <- NROW(X1) # number of observations
  if(is.null(k)){
    k<-floor(n ^ 0.4)
  }
  #Checks
  if (threshold_q>=1 | threshold_q<=0) {
    stop("threshold_q must be between 0 and 1.")
  }
  if (k <= 1 | k>floor((1-threshold_q)*n) ) {
    stop("k must be greater than 1 and smaller than the number of threshold exceedences.")
  }
  if (both_tails) {
    stop("both_tails not yet possible.")
  }
  
  Hcs <- scale(H) #Centering and scaling H
  
  r1 <- rank(X1, ties.method = "first") # ranks of X1
  r2 <- rank(X2, ties.method = "first") # ranks of X2
  
  #Fit POT GPD model and compute the fitted cdf in the data points for X2
  u2 <- quantile(X2,threshold_q)
  fit2 <- gpd_constraints_fit_optim(X2,threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE, method = "L-BFGS-B",
                                    lower=c(-Inf,-u2/(t_df*max(Hcs)),-Inf), upper=c(Inf,-u2/(t_df*min(Hcs)),Inf))
  FX2 <- pgpd(X2, loc=u2, scale=fit2$mle[1]+fit2$mle[2]*Hcs, shape=fit2$mle[3])
  FX2 <- (FX2*(1-threshold_q)+threshold_q)*(X2>=u2) + r2/n * (X2<u2)
  
  #Compute the parametric causal tail coefficient
  if(parametric_F1){
    #Fit POT GPD model and compute the fitted cdf in the data points for X1
    u1 <- quantile(X1,threshold_q)
    fit1 <- gpd_constraints_fit_optim(X1,threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE, method = "L-BFGS-B",
                                      lower=c(-Inf,-u1/(t_df*max(Hcs)),-Inf), upper=c(Inf,-u1/(t_df*min(Hcs)),Inf))
    FX1 <- pgpd(X1, loc=u1, scale=fit1$mle[1]+fit1$mle[2]*Hcs, shape=fit1$mle[3])
    FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
    #Compute the parametric causal tail coefficient
    ctc <- mean(FX2[FX1 > (1-k/n)])
    if(return_MLEs){
      return(list(ctc=ctc, F2_scale0=fit2$mle[1], F2_scale1=fit2$mle[2], F2_shape=fit2$mle[3],
                  F1_scale0=fit1$mle[1], F1_scale1=fit1$mle[2], F1_shape=fit1$mle[3]))
    }
  }else{
    ctc <- 1/k * sum(FX2[r1 > n - k])
    if(return_MLEs){
      return(list(ctc=ctc, F2_scale0=fit2$mle[1], F2_scale1=fit2$mle[2], F2_shape=fit2$mle[3],
                  F1_scale0=NULL, F1_scale1=NULL, F1_shape=NULL))
    }
  }
  return(ctc)
}

# Non-linear GPD CTC estimators

expGPD_causal_tail_coeff <- function(X1, X2, H, k = NULL, threshold_q = 0.95, parametric_F1=TRUE, both_tails = FALSE){
  ## H-conditional GPD causal tail coefficient estimator with exponential link function
  n <- NROW(X1) # number of observations
  if(is.null(k)){
    k<-floor(n ^ 0.4)
  }
  #Checks
  if (threshold_q>=1 | threshold_q<=0) {
    stop("threshold_q must be between 0 and 1.")
  }
  if (k <= 1 | k>floor((1-threshold_q)*n) ) {
    stop("k must be greater than 1 and smaller than the number of threshold exceedences.")
  }
  if (both_tails) {
    stop("both_tails not yet possible.")
  }
  
  Hcs <- scale(H) #Centering and scaling H
  
  r1 <- rank(X1, ties.method = "first") # ranks of X1
  r2 <- rank(X2, ties.method = "first") # ranks of X2
  
  #Fit POT GPD model and compute the fitted cdf in the data points for X2
  u2 <- quantile(X2,threshold_q)
  fit2 <- gpd.fit(X2,threshold=u2,ydat=as.matrix(Hcs),sigl=1, siglink=exp, show=FALSE)
  FX2 <- pgpd(X2, loc=u2, scale=exp(fit2$mle[1]+fit2$mle[2]*Hcs), shape=fit2$mle[3])
  FX2 <- (FX2*(1-threshold_q)+threshold_q)*(X2>=u2) + r2/n * (X2<u2)
  
  #Compute the parametric causal tail coefficient
  if(parametric_F1){
    #Fit POT GPD model and compute the fitted cdf in the data points for X1
    u1 <- quantile(X1,threshold_q)
    fit1 <- gpd.fit(X1,threshold=u1,ydat=as.matrix(Hcs),sigl=1, siglink=exp, show=FALSE)
    FX1 <- pgpd(X1, loc=u1, scale=exp(fit1$mle[1]+fit1$mle[2]*Hcs), shape=fit1$mle[3])
    FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
    #Compute the parametric causal tail coefficient
    ctc <- mean(FX2[FX1 > (1-k/n)])
  }else{
    ctc <- 1/k * sum(FX2[r1 > n - k])
  }
  return(ctc)
}


## CAUSALITY PERMUTATION MONTE CARLO TEST

CTC_causality_permutation_test <- function(X1, X2, H=NULL, k=floor(NROW(X1)^0.4), R=10000, constrained_fit=FALSE, threshold_q = 0.95, both_tails=FALSE, parametric_F1=TRUE){
  n <- NROW(X1) # number of observations
  if(is.null(k)){
    k<-floor(n ^ 0.4)
  }
  #Checks
  if (threshold_q>=1 | threshold_q<=0) {
    stop("threshold_q must be between 0 and 1.")
  }
  if (k <= 1 | k>floor((1-threshold_q)*n) ) {
    stop("k must be greater than 1 and smaller than the number of threshold exceedences.")
  }
  if (both_tails) {
    stop("Both tails CTC not possible in CTC_causality_permutation_test.")
  }
  if(!parametric_F1){stop("Non-Parametric F1 not possible for the LGPD CTC in CTC_causality_permutation_test.")}
  
  r1 <- rank(X1, ties.method = "first") # ranks of X1
  r2 <- rank(X2, ties.method = "first") # ranks of X2
  
  # Fit and transform F1(X1) F2(X2)
  
  if(is.null(H)){# Non-parametric CTC
    FX1 <- r1/n
    FX2 <- r2/n
    
  } else {# LGPD CTC
    Hcs <- scale(H) #Centering and scaling H
    u2 <- quantile(X2,threshold_q)
    u1 <- quantile(X1,threshold_q)
    
    if(constrained_fit){
      #Fit POT GPD model and compute the fitted cdf in the data points for X2
      fit2 <- gpd_constraints_fit_maxLik(X2,threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE,
                                         constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
      FX2 <- pgpd(X2, loc=u2, scale=fit2$mle[1]+fit2$mle[2]*Hcs, shape=fit2$mle[3])
      FX2 <- (FX2*(1-threshold_q)+threshold_q)*(X2>=u2) + r2/n * (X2<u2)
      #Fit POT GPD model and compute the fitted cdf in the data points for X1
      fit1 <- gpd_constraints_fit_maxLik(X1,threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE,
                                         constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
      FX1 <- pgpd(X1, loc=u1, scale=fit1$mle[1]+fit1$mle[2]*Hcs, shape=fit1$mle[3])
      FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
      
    } else {# Post-fit corrected LGPD CTC
      #Fit POT GPD model and compute the fitted cdf in the data points for X2
      fit2 <- gpd.fit(X2,threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE)
      FX2 <- pgpd(X2, loc=u2, scale=pmax(fit2$mle[1]+fit2$mle[2]*Hcs,0.0001)[1:length(Hcs)], shape=fit2$mle[3])
      FX2 <- (FX2*(1-threshold_q)+threshold_q)*(X2>=u2) + r2/n * (X2<u2)
      #Fit POT GPD model and compute the fitted cdf in the data points for X1
      fit1 <- gpd.fit(X1,threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE)
      FX1 <- pgpd(X1, loc=u1, scale=pmax(fit1$mle[1]+fit1$mle[2]*Hcs,0.0001)[1:length(Hcs)], shape=fit1$mle[3])
      FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
    }
  }
  # Compute the coefficients using the F1(X1) and F2(X2)
  ctc12 <- mean(FX2[FX1 > (1-k/n)])
  ctc21 <- mean(FX1[FX2 > (1-k/n)])
  diff12 <- ctc12 - ctc21
  
  # Create R datasets with random indices permutations
  # and Compute the coefficients using the F1(X1) and F2(X2) for each sample
  ctc12p <- rep(as.double(NA), R)
  ctc21p <- rep(as.double(NA), R)
  for(i in 1:R){
    FX1p <- rep(as.double(NA), n)
    FX2p <- rep(as.double(NA), n)
    permbool <- sample(x=c(TRUE,FALSE), size=n, replace=TRUE)
    FX1p[!permbool] <- FX1[!permbool]
    FX1p[permbool] <- FX2[permbool]
    FX2p[!permbool] <- FX2[!permbool]
    FX2p[permbool] <- FX1[permbool]
    ctc12p[i] <- mean(FX2p[FX1p > (1-k/n)])
    ctc21p[i] <- mean(FX1p[FX2p > (1-k/n)])
  }
  diff12p <- ctc12p - ctc21p
  
  # Monte Carlo test (p-value estimate)
  Pmc <- (1+sum(diff12p>=diff12))/(R+1)
  
  return(list(Pmc=Pmc, ctc12=ctc12, ctc21=ctc21, diff12=diff12, ctc12p=ctc12p, ctc21p=ctc21p, diff12p=diff12p))
}

expGPD_CTC_causality_permutation_test <- function(X1, X2, H, k=floor(NROW(X1)^0.4), R=10000, threshold_q = 0.95, both_tails=FALSE, parametric_F1=TRUE){
  n <- NROW(X1) # number of observations
  if(is.null(k)){
    k<-floor(n ^ 0.4)
  }
  #Checks
  if (threshold_q>=1 | threshold_q<=0) {
    stop("threshold_q must be between 0 and 1.")
  }
  if (k <= 1 | k>floor((1-threshold_q)*n) ) {
    stop("k must be greater than 1 and smaller than the number of threshold exceedences.")
  }
  if (both_tails) {
    stop("Both tails CTC not possible in CTC_causality_permutation_test.")
  }
  if(!parametric_F1){stop("Non-Parametric F1 not possible for the LGPD CTC in CTC_causality_permutation_test.")}
  
  r1 <- rank(X1, ties.method = "first") # ranks of X1
  r2 <- rank(X2, ties.method = "first") # ranks of X2
  
  # Fit and transform F1(X1) F2(X2)
  
  Hcs <- scale(H) #Centering and scaling H
  u2 <- quantile(X2,threshold_q)
  u1 <- quantile(X1,threshold_q)
  
  #Fit POT GPD model and compute the fitted cdf in the data points for X2
  fit2 <- gpd.fit(X2,threshold=u2,ydat=as.matrix(Hcs),sigl=1, siglink=exp, show=FALSE)
  FX2 <- pgpd(X2, loc=u2, scale=exp(fit2$mle[1]+fit2$mle[2]*Hcs), shape=fit2$mle[3])
  FX2 <- (FX2*(1-threshold_q)+threshold_q)*(X2>=u2) + r2/n * (X2<u2)
  #Fit POT GPD model and compute the fitted cdf in the data points for X1
  fit1 <- gpd.fit(X1,threshold=u1,ydat=as.matrix(Hcs),sigl=1, siglink=exp, show=FALSE)
  FX1 <- pgpd(X1, loc=u1, scale=exp(fit1$mle[1]+fit1$mle[2]*Hcs), shape=fit1$mle[3])
  FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
  
  # Compute the coefficients using the F1(X1) and F2(X2)
  ctc12 <- mean(FX2[FX1 > (1-k/n)])
  ctc21 <- mean(FX1[FX2 > (1-k/n)])
  diff12 <- ctc12 - ctc21
  
  # Create R datasets with random indices permutations
  # and Compute the coefficients using the F1(X1) and F2(X2) for each sample
  ctc12p <- rep(as.double(NA), R)
  ctc21p <- rep(as.double(NA), R)
  for(i in 1:R){
    FX1p <- rep(as.double(NA), n)
    FX2p <- rep(as.double(NA), n)
    permbool <- sample(x=c(TRUE,FALSE), size=n, replace=TRUE)
    FX1p[!permbool] <- FX1[!permbool]
    FX1p[permbool] <- FX2[permbool]
    FX2p[!permbool] <- FX2[!permbool]
    FX2p[permbool] <- FX1[permbool]
    ctc12p[i] <- mean(FX2p[FX1p > (1-k/n)])
    ctc21p[i] <- mean(FX1p[FX2p > (1-k/n)])
  }
  diff12p <- ctc12p - ctc21p
  
  # Monte Carlo test
  Pmc <- (1+sum(diff12p>=diff12))/(R+1)
  
  return(list(Pmc=Pmc, ctc12=ctc12, ctc21=ctc21, diff12=diff12, ctc12p=ctc12p, ctc21p=ctc21p, diff12p=diff12p))
}



