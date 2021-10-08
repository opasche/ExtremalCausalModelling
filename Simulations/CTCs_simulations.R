#library(tidyverse)
library("evd")
library("ismev")
#library(latex2exp)
library(causalXtreme)
library(EnvStats)
#library(boot)
#library("ggpubr")
#theme_set(theme_bw() + theme(legend.position = "top"))

source("../R_functions/CTC_contribution_functions_OCP.R")

set.seed(1)

for(sims_round in 1:1){
  
  start_time <- Sys.time()
  
  ## Simulations
  cat("\n=================== SIMULATIONS: CTCs distributions (",sims_round,") ===================\n", sep="")
  t_df <- 4
  n_sim <- 1e3
  size_sim <- 1e5
  thresh_q <- 0.9
  param.k <- 2*floor(size_sim^0.4) #floor(size_sim^0.4)
  
  save_path <- paste0("CTCs_simulations/last/CTCs_t4_H3_q90_100k2k1000", "_", sims_round, ".RData")
  check_directory(dirname(save_path), recursive=TRUE)
  
  H <- matrix(rt(n_sim*size_sim, df=3), c(n_sim,size_sim))
  X1.i <- matrix(rt(n_sim*size_sim, df=4), c(n_sim,size_sim))
  X2.i <- matrix(rt(n_sim*size_sim, df=4), c(n_sim,size_sim))
  #rlnorm(n_sim*size_sim, meanlog=0, sdlog=1)#rt(n_sim*size_sim, df=t_df)#rpareto(n_sim*size_sim, location=1, shape=2)
  
  cat("X1 indep. X2:\n")
  ctc12.i <- rep(as.double(NA), n_sim)# Pre-attribute memory space for faster loop
  ctc21.i <- rep(as.double(NA), n_sim)
  gpdctc12.i <- rep(as.double(NA), n_sim)
  gpdctc21.i <- rep(as.double(NA), n_sim)
  pfcorr_lgpdctc12.i <- rep(as.double(NA), n_sim)
  pfcorr_lgpdctc21.i <- rep(as.double(NA), n_sim)
  constr_lgpdctc12.i <- rep(as.double(NA), n_sim)
  constr_lgpdctc21.i <- rep(as.double(NA), n_sim)
  for (i in 1:n_sim){
    ctc12.i[i] <- causal_tail_coeff(X1.i[i,], X2.i[i,], k=param.k, both_tails = FALSE)
    ctc21.i[i] <- causal_tail_coeff(X2.i[i,], X1.i[i,], k=param.k, both_tails = FALSE)
    gpdctc12.i[i] <- GPD_causal_tail_coeff(X1.i[i,], X2.i[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    gpdctc21.i[i] <- GPD_causal_tail_coeff(X2.i[i,], X1.i[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    pfcorr_lgpdctc12.i[i] <- LGPD_causal_tail_coeff(X1.i[i,], X2.i[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    pfcorr_lgpdctc21.i[i] <- LGPD_causal_tail_coeff(X2.i[i,], X1.i[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    constr_lgpdctc12.i[i] <- constrained_LGPD_CTC(X1.i[i,], X2.i[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    constr_lgpdctc21.i[i] <- constrained_LGPD_CTC(X2.i[i,], X1.i[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
  }
  cat("X1 X2 with confounder:\n")
  X1.ih <- H + X1.i
  X2.ih <- H + X2.i
  ctc12.ih <- rep(as.double(NA), n_sim)
  ctc21.ih <- rep(as.double(NA), n_sim)
  gpdctc12.ih <- rep(as.double(NA), n_sim)
  gpdctc21.ih <- rep(as.double(NA), n_sim)
  pfcorr_lgpdctc12.ih <- rep(as.double(NA), n_sim)
  pfcorr_lgpdctc21.ih <- rep(as.double(NA), n_sim)
  constr_lgpdctc12.ih <- rep(as.double(NA), n_sim)
  constr_lgpdctc21.ih <- rep(as.double(NA), n_sim)
  for (i in 1:n_sim){
    ctc12.ih[i] <- causal_tail_coeff(X1.ih[i,], X2.ih[i,], k=param.k, both_tails = FALSE)
    ctc21.ih[i] <- causal_tail_coeff(X2.ih[i,], X1.ih[i,], k=param.k, both_tails = FALSE)
    gpdctc12.ih[i] <- GPD_causal_tail_coeff(X1.ih[i,], X2.ih[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    gpdctc21.ih[i] <- GPD_causal_tail_coeff(X2.ih[i,], X1.ih[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    pfcorr_lgpdctc12.ih[i] <- LGPD_causal_tail_coeff(X1.ih[i,], X2.ih[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    pfcorr_lgpdctc21.ih[i] <- LGPD_causal_tail_coeff(X2.ih[i,], X1.ih[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    constr_lgpdctc12.ih[i] <- constrained_LGPD_CTC(X1.ih[i,], X2.ih[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    constr_lgpdctc21.ih[i] <- constrained_LGPD_CTC(X2.ih[i,], X1.ih[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
  }
  cat("X1 -> X2:\n")
  rm(X2.ih)# save space
  X1.c <- X1.i# no copy (same tracemem())
  X2.c <- X1.c + X2.i
  ctc12.c <- rep(as.double(NA), n_sim)
  ctc21.c <- rep(as.double(NA), n_sim)
  gpdctc12.c <- rep(as.double(NA), n_sim)
  gpdctc21.c <- rep(as.double(NA), n_sim)
  pfcorr_lgpdctc12.c <- rep(as.double(NA), n_sim)
  pfcorr_lgpdctc21.c <- rep(as.double(NA), n_sim)
  constr_lgpdctc12.c <- rep(as.double(NA), n_sim)
  constr_lgpdctc21.c <- rep(as.double(NA), n_sim)
  for (i in 1:n_sim){
    ctc12.c[i] <- causal_tail_coeff(X1.c[i,], X2.c[i,], k=param.k, both_tails = FALSE)
    ctc21.c[i] <- causal_tail_coeff(X2.c[i,], X1.c[i,], k=param.k, both_tails = FALSE)
    gpdctc12.c[i] <- GPD_causal_tail_coeff(X1.c[i,], X2.c[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    gpdctc21.c[i] <- GPD_causal_tail_coeff(X2.c[i,], X1.c[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    pfcorr_lgpdctc12.c[i] <- LGPD_causal_tail_coeff(X1.c[i,], X2.c[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    pfcorr_lgpdctc21.c[i] <- LGPD_causal_tail_coeff(X2.c[i,], X1.c[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    constr_lgpdctc12.c[i] <- constrained_LGPD_CTC(X1.c[i,], X2.c[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    constr_lgpdctc21.c[i] <- constrained_LGPD_CTC(X2.c[i,], X1.c[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
  }
  cat("X1 -> X2 with confounder:\n")
  rm(X1.c,X2.c)# save space
  X1.ch <- X1.ih# no copy (same tracemem())
  X2.ch <- X1.ch + H + X2.i
  rm(X1.i,X2.i)
  ctc12.ch <- rep(as.double(NA), n_sim)
  ctc21.ch <- rep(as.double(NA), n_sim)
  gpdctc12.ch <- rep(as.double(NA), n_sim)
  gpdctc21.ch <- rep(as.double(NA), n_sim)
  pfcorr_lgpdctc12.ch <- rep(as.double(NA), n_sim)
  pfcorr_lgpdctc21.ch <- rep(as.double(NA), n_sim)
  constr_lgpdctc12.ch <- rep(as.double(NA), n_sim)
  constr_lgpdctc21.ch <- rep(as.double(NA), n_sim)
  for (i in 1:n_sim){
    ctc12.ch[i] <- causal_tail_coeff(X1.ch[i,], X2.ch[i,], k=param.k, both_tails = FALSE)
    ctc21.ch[i] <- causal_tail_coeff(X2.ch[i,], X1.ch[i,], k=param.k, both_tails = FALSE)
    gpdctc12.ch[i] <- GPD_causal_tail_coeff(X1.ch[i,], X2.ch[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    gpdctc21.ch[i] <- GPD_causal_tail_coeff(X2.ch[i,], X1.ch[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    pfcorr_lgpdctc12.ch[i] <- LGPD_causal_tail_coeff(X1.ch[i,], X2.ch[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    pfcorr_lgpdctc21.ch[i] <- LGPD_causal_tail_coeff(X2.ch[i,], X1.ch[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    constr_lgpdctc12.ch[i] <- constrained_LGPD_CTC(X1.ch[i,], X2.ch[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
    constr_lgpdctc21.ch[i] <- constrained_LGPD_CTC(X2.ch[i,], X1.ch[i,], H[i,], k=param.k, threshold_q=thresh_q, parametric_F1=TRUE)
  }
  
  if(!is.null(save_path)){
    save(ctc12.i,ctc21.i,ctc12.ih,ctc21.ih,ctc12.c,ctc21.c,ctc12.ch,ctc21.ch,
         gpdctc12.i,gpdctc21.i,gpdctc12.ih,gpdctc21.ih,gpdctc12.c,gpdctc21.c,gpdctc12.ch,gpdctc21.ch,
         pfcorr_lgpdctc12.i,pfcorr_lgpdctc21.i,pfcorr_lgpdctc12.ih,pfcorr_lgpdctc21.ih,
         pfcorr_lgpdctc12.c,pfcorr_lgpdctc21.c,pfcorr_lgpdctc12.ch,pfcorr_lgpdctc21.ch,
         constr_lgpdctc12.i,constr_lgpdctc21.i,constr_lgpdctc12.ih,constr_lgpdctc21.ih,
         constr_lgpdctc12.c,constr_lgpdctc21.c,constr_lgpdctc12.ch,constr_lgpdctc21.ch,
         file=save_path)
  }
  
  end_time <- Sys.time()
  cat("\nRun time:\n")
  print(end_time - start_time)
  
  rm(list = setdiff(ls(), union(lsf.str(), "sims_round")))
}

#Sys.sleep(300)
#system('shutdown -t 30 -s')

