
## DEPRECATED

## Partial Confidence Intervals and diagnostics for Causal Tail Coefficients

causal_tail_coeff_CI <- function(v1, v2, k = floor(n ^ 0.4), alpha_CI = 0.05, to_rank = TRUE, both_tails = FALSE){
  ## Partially using the original 'causal_tail_coeff' function implementation from:
  ## [Nicola Gnecco, Nicolai Meinshausen, Jonas Peters, and Sebastian Engelke. Causal discovery in heavy-tailed models, 2019.]
  # number of observations
  n <- NROW(v1)
  
  # check k
  if (k <= 1 | k >= n) {
    stop("k must be greater than 1 and smaller than n.")
  }
  
  # rank variables?
  if (to_rank){
    r1 <- rank(v1, ties.method = "first") # ranks of v1
    r2 <- rank(v2, ties.method = "first") # ranks of v2
  } else{
    r1 <- v1
    r2 <- v2
  }
  
  #Confidence interval preparation
  DKW_eps <- sqrt(log(2/alpha_CI)/(2*n))
  
  # compute causal tail coefficient
  if (both_tails){
    k <- (k %/% 2) * 2
    ctc <- 1 / (k * n) * sum(2 * abs(r2[r1 > n - k / 2 | r1 <= k / 2] - (n + 1) / 2))
    warning("CI for \"both_tails\" not yet implemented, returning NAs")
    ctc_up <- NA #not 1 / k * sum(2 * abs(pmin(r2/n+DKW_eps,1)[r1 > n - k / 2 | r1 <= k / 2] - (n + 1) / (2*n)))
    ctc_down <- NA #not 1 / k * sum(2 * abs(pmax(r2/n-DKW_eps,0)[r1 > n - k / 2 | r1 <= k / 2] - (n + 1) / (2*n)))
  } else{
    ctc <- 1 / (k * n) * sum(r2[r1 > n - k])
    ctc_down <- 1 / k * sum(pmax(r2/n-DKW_eps,0)[r1 > n - k])
    ctc_up <- 1 / k * sum(pmin(r2/n+DKW_eps,1)[r1 > n - k])
  }
  return(c(CI_down=ctc_down, estim=ctc, CI_up=ctc_up))
}

GPD_causal_tail_coeff_CI <- function(X1, X2, k = NULL, threshold_q = 0.95, parametric_F1=TRUE, alpha_CI = 0.05, both_tails = FALSE){
  n <- NROW(X1) # number of observations
  if(is.null(k)){
    k<-floor(n ^ 0.4)
  }
  # checks
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
  
  #Confidence interval preparation
  DKW_eps <- sqrt(log(2/alpha_CI)/(2*n))
  sigm <- fit2$estimate[1]
  xi <- fit2$estimate[2]
  CI_GPD2 <- function(y){
    if(y<=u2){return(0)}
    dF1 <- -xi*(y-u2)*(1+xi*(y-u2)/sigm)^(-(1+xi)/xi)/(xi*sigm^2)
    dF2 <- ((1+xi*(y-u2)/sigm)^(-(1+xi)/xi)*(xi*(y-u2)/sigm-(1+xi*(y-u2)/sigm)*log(1+xi*(y-u2)/sigm)))/(xi^2)
    return(qnorm(1-alpha_CI/2)*sqrt((matrix(c(dF1, dF2), nrow = 1, ncol = 2) %*% fit2$var.cov %*% matrix(c(dF1, dF2), nrow = 2, ncol = 1))[1,1]))
  }
  CI_F2_eps <- function(x){(CI_GPD2(x)*(1-threshold_q))*(x>=u2) + DKW_eps * (x<u2)}
  CI_eps2 <- NULL
  for(i in 1:n){
    CI_eps2[i] <- CI_F2_eps(X2[i])
  }
  CI_F2_up <- pmin(FX2 + CI_eps2,1)
  CI_F2_down <- pmax(FX2 - CI_eps2,0)
  
  #Compute the parametric causal tail coefficient and the CIs
  if(parametric_F1){
    #Fit POT GPD model and compute the fitted cdf in the data points for X1
    u1 <- quantile(X1,threshold_q)
    fit1<-fpot(X1,threshold=u1)
    FX1 <- pgpd(X1, loc=u1, scale=fit1$estimate[1], shape=fit1$estimate[2])
    FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
    #Compute the parametric causal tail coefficient
    ctc <- mean(FX2[FX1 > (1-k/n)])
    ctc_down <- mean(CI_F2_down[FX1 > (1-k/n)])
    ctc_up <- mean(CI_F2_up[FX1 > (1-k/n)])
  }else{
    ctc <- 1/k * sum(FX2[r1 > n - k])
    ctc_down <- 1/k * sum(CI_F2_down[r1 > n - k])
    ctc_up <- 1/k * sum(CI_F2_up[r1 > n - k])
  }
  return(c(CI_down=ctc_down, estim=ctc, CI_up=ctc_up))
}

LGPD_causal_tail_coeff_CI <- function(X1, X2, H, k = NULL, threshold_q = 0.95, parametric_F1=TRUE, alpha_CI = 0.05, both_tails = FALSE){
  n <- NROW(X1) # number of observations
  if(is.null(k)){
    k<-floor(n ^ 0.4)
  }
  # checks
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
  
  #Confidence interval preparation
  DKW_eps <- sqrt(log(2/alpha_CI)/(2*n))
  xi <- fit2$mle[3]
  CI_GPD2 <- function(y,h){
    if(y<=u2){return(0)}
    sigm <- max(fit2$mle[1]+fit2$mle[2]*h,0.0001) #max(fit2$mle[1]+fit2$mle[2]*h,0.0001)
    dF1 <- -xi*(y-u2)*(1+xi*(y-u2)/sigm)^(-(1+xi)/xi)/(xi*sigm^2)
    dF2 <- h*dF1
    dF3 <- ((1+xi*(y-u2)/sigm)^(-(1+xi)/xi)*(xi*(y-u2)/sigm-(1+xi*(y-u2)/sigm)*log(1+xi*(y-u2)/sigm)))/(xi^2) #max(sigm,0.0001)
    return(qnorm(1-alpha_CI/2)*sqrt((matrix(c(dF1, dF2, dF3), nrow = 1, ncol = 3) %*% as.matrix(fit2$cov) %*% matrix(c(dF1, dF2, dF3), nrow = 3, ncol = 1))[1,1]))
  }
  CI_F2_eps <- function(x,h){(CI_GPD2(x,h)*(1-threshold_q))*(x>=u2) + DKW_eps * (x<u2)}
  CI_eps2 <- NULL
  for(i in 1:n){
    CI_eps2[i] <- CI_F2_eps(X2[i],Hcs[i])
  }
  CI_F2_up <- pmin(FX2 + CI_eps2,1)
  CI_F2_down <- pmax(FX2 - CI_eps2,0)
  
  #Compute the parametric causal tail coefficient and the CIs
  if(parametric_F1){
    #Fit POT GPD model and compute the fitted cdf in the data points for X1
    u1 <- quantile(X1,threshold_q)
    fit1 <- gpd.fit(X1,threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE)
    FX1 <- pgpd(X1, loc=u1, scale=pmax(fit1$mle[1]+fit1$mle[2]*Hcs,0.0001)[1:length(Hcs)], shape=fit1$mle[3])
    FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
    #Compute the parametric causal tail coefficient
    ctc <- mean(FX2[FX1 > (1-k/n)])
    ctc_down <- mean(CI_F2_down[FX1 > (1-k/n)])
    ctc_up <- mean(CI_F2_up[FX1 > (1-k/n)])
  }else{
    ctc <- 1/k * sum(FX2[r1 > n - k])
    ctc_down <- 1/k * sum(CI_F2_down[r1 > n - k])
    ctc_up <- 1/k * sum(CI_F2_up[r1 > n - k])
  }
  return(c(CI_down=ctc_down, estim=ctc, CI_up=ctc_up))
}


LGPD_causal_tail_coeff_diagnostic <- function(X1, X2, H, k = NULL, threshold_q = 0.95, parametric_F1=TRUE, alpha_CI = 0.05, both_tails = FALSE){
  n <- NROW(X1) # number of observations
  if(is.null(k)){
    k<-floor(n ^ 0.4)
  }
  # checks
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
  
  #Confidence interval preparation
  DKW_eps <- sqrt(log(2/alpha_CI)/(2*n))
  xi <- fit2$mle[3]
  CI_GPD2 <- function(y,h){
    if(y<=u2){return(0)}
    sigm <- max(fit2$mle[1]+fit2$mle[2]*h,0.0001) #max(fit2$mle[1]+fit2$mle[2]*h,0.0001)
    dF1 <- -xi*(y-u2)*(1+xi*(y-u2)/sigm)^(-(1+xi)/xi)/(xi*sigm^2)
    dF2 <- h*dF1
    dF3 <- ((1+xi*(y-u2)/sigm)^(-(1+xi)/xi)*(xi*(y-u2)/sigm-(1+xi*(y-u2)/sigm)*log(1+xi*(y-u2)/sigm)))/(xi^2) #max(sigm,0.0001)
    return(qnorm(1-alpha_CI/2)*sqrt((matrix(c(dF1, dF2, dF3), nrow = 1, ncol = 3) %*% as.matrix(fit2$cov) %*% matrix(c(dF1, dF2, dF3), nrow = 3, ncol = 1))[1,1]))
  }
  CI_F2_eps <- function(x,h){(CI_GPD2(x,h)*(1-threshold_q))*(x>=u2) + DKW_eps * (x<u2)}
  CI_eps2 <- NULL
  for(i in 1:n){
    CI_eps2[i] <- CI_F2_eps(X2[i],Hcs[i])
  }
  CI_F2_up <- pmin(FX2 + CI_eps2,1)
  CI_F2_down <- pmax(FX2 - CI_eps2,0)
  
  #Compute the parametric causal tail coefficient and the CIs
  if(parametric_F1){
    #Fit POT GPD model and compute the fitted cdf in the data points for X1
    u1 <- quantile(X1,threshold_q)
    fit1 <- gpd.fit(X1,threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE)#,siglink=exp)
    FX1 <- pgpd(X1, loc=u1, scale=pmax(fit1$mle[1]+fit1$mle[2]*Hcs,0.0001)[1:length(Hcs)], shape=fit1$mle[3])
    FX1 <- (FX1*(1-threshold_q)+threshold_q)*(X1>=u1) + r1/n * (X1<u1)
    #Compute the parametric causal tail coefficient
    ctc <- mean(FX2[FX1 > (1-k/n)])
    ctc_down <- mean(CI_F2_down[FX1 > (1-k/n)])
    ctc_up <- mean(CI_F2_up[FX1 > (1-k/n)])
    return(list(estim=ctc, CI_down=ctc_down, CI_up=ctc_up,
                F2_scale0=fit2$mle[1], F2_scale0_se=fit2$se[1], F2_scale1=fit2$mle[2], F2_scale1_se=fit2$se[2], F2_shape=fit2$mle[3], F2_shape_se=fit2$se[3],
                F1_scale0=fit1$mle[1], F1_scale0_se=fit1$se[1], F1_scale1=fit1$mle[2], F1_scale1_se=fit1$se[2], F1_shape=fit1$mle[3], F1_shape_se=fit1$se[3],
                u2=u2, u1=u1, conv2=fit2$conv, conv1=fit1$conv))
  }else{
    u1 <- NULL
    fit1 <- NULL
    ctc <- 1/k * sum(FX2[r1 > n - k])
    ctc_down <- 1/k * sum(CI_F2_down[r1 > n - k])
    ctc_up <- 1/k * sum(CI_F2_up[r1 > n - k])
    return(list(estim=ctc, CI_down=ctc_down, CI_up=ctc_up,
                F2_scale0=fit2$mle[1], F2_scale0_se=fit2$se[1], F2_scale1=fit2$mle[2], F2_scale1_se=fit2$se[2], F2_shape=fit2$mle[3], F2_shape_se=fit2$se[3],
                u2=u2, conv2=fit2$conv))
  }
}

