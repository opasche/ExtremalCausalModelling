library(tidyverse)
library("evd")
library("ismev")
library(causalXtreme)
library(EnvStats)
library(parallel)
library(doParallel)

start_time <- Sys.time()

#sink("log_last.txt")
cat("\n=================== SIMULATIONS: CTC Permutation Causality MC Test ===================\n", sep="")

parallel_loop <- TRUE

betas <- c(0, 0.01, 0.05, 0.1, 0.2) #sequence of causal strengths of X1 on X2

set.seed(1)

t_df <- 4
n_sim <- 1e3
size_sim <- 1e4
thresh_q <- 0.9
param.k <- 2*floor(size_sim^0.4)
param.R <- 1000

H <- matrix(rt(n_sim*size_sim, df=4), c(n_sim,size_sim))
X1.i <- matrix(rt(n_sim*size_sim, df=4), c(n_sim,size_sim))
X2.i <- matrix(rt(n_sim*size_sim, df=4), c(n_sim,size_sim))
#rlnorm(n_sim*size_sim, meanlog=0, sdlog=1)#rt(n_sim*size_sim, df=t_df)#rpareto(n_sim*size_sim, location=1, shape=2)

nb_cores_total <- detectCores()#12
nb_cores_used <- min(nb_cores_total-1, length(betas))

if(parallel_loop){
  cl <- makeCluster(nb_cores_used)
  registerDoParallel(cl)
  `%fun%` <- `%dopar%`
} else {
  `%fun%` <- `%do%`
}

packages_to_pass <- c("tidyverse","evd","ismev","causalXtreme","EnvStats")
objects_to_pass <- c("t_df", "n_sim", "size_sim", "thresh_q", "param.k", "param.R", "H", "X1.i", "X2.i")

res_tbl <- foreach(beta21=betas, .packages=packages_to_pass, .export=objects_to_pass, .errorhandling="stop", .combine=rbind) %fun% {
  
  set.seed(2)
  source("../R_functions/CTC_contribution_functions_OCP.R")
  
  save_path <- paste0("CTCs_simulations/Permutation_tests/",
                      "CTCs_Permutation_test_power_t4_q90_10k2k_R1k_beta", str_replace(toString(beta21),"([.])","p"), "_1000", ".RData")
  check_directory(dirname(save_path), recursive=TRUE)
  
  beta1H <- 1
  beta2H <- 1
  
  
  cat("X1 -> X2:\n")
  X1.c <- X1.i# no copy (same tracemem())
  X2.c <- beta21*X1.c + X2.i
  res_ctc.c <- list()
  Pmc_ctc.c <- rep(as.double(NA), n_sim)
  res_pfcorr_lgpdctc.c <- list()
  Pmc_pfcorr_lgpdctc.c <- rep(as.double(NA), n_sim)
  res_constr_lgpdctc.c <- list()
  Pmc_constr_lgpdctc.c <- rep(as.double(NA), n_sim)
  for (i in 1:n_sim){
    res_ctc.c[[i]] <- CTC_causality_permutation_test(X1.c[i,], X2.c[i,], H=NULL, k=param.k, R=param.R, constrained_fit=FALSE, threshold_q=thresh_q)
    Pmc_ctc.c[i] <- res_ctc.c[[i]]$Pmc
    res_pfcorr_lgpdctc.c[[i]] <- CTC_causality_permutation_test(X1.c[i,], X2.c[i,], H=H[i,], k=param.k, R=param.R, constrained_fit=FALSE, threshold_q=thresh_q)
    Pmc_pfcorr_lgpdctc.c[i] <- res_pfcorr_lgpdctc.c[[i]]$Pmc
    res_constr_lgpdctc.c[[i]] <- CTC_causality_permutation_test(X1.c[i,], X2.c[i,], H=H[i,], k=param.k, R=param.R, constrained_fit=TRUE, threshold_q=thresh_q)
    Pmc_constr_lgpdctc.c[i] <- res_constr_lgpdctc.c[[i]]$Pmc
  }
  cat("X1 -> X2 with confounder:\n")
  rm(X1.c,X2.c)# save space
  X1.ch <- X1.i + beta1H*H
  X2.ch <- beta21*X1.ch + beta2H*H + X2.i
  #rm(X1.i,X2.i)
  res_ctc.ch <- list()
  Pmc_ctc.ch <- rep(as.double(NA), n_sim)
  res_pfcorr_lgpdctc.ch <- list()
  Pmc_pfcorr_lgpdctc.ch <- rep(as.double(NA), n_sim)
  res_constr_lgpdctc.ch <- list()
  Pmc_constr_lgpdctc.ch <- rep(as.double(NA), n_sim)
  for (i in 1:n_sim){
    res_ctc.ch[[i]] <- CTC_causality_permutation_test(X1.ch[i,], X2.ch[i,], H=NULL, k=param.k, R=param.R, constrained_fit=FALSE, threshold_q=thresh_q)
    Pmc_ctc.ch[i] <- res_ctc.ch[[i]]$Pmc
    res_pfcorr_lgpdctc.ch[[i]] <- CTC_causality_permutation_test(X1.ch[i,], X2.ch[i,], H=H[i,], k=param.k, R=param.R, constrained_fit=FALSE, threshold_q=thresh_q)
    Pmc_pfcorr_lgpdctc.ch[i] <- res_pfcorr_lgpdctc.ch[[i]]$Pmc
    res_constr_lgpdctc.ch[[i]] <- CTC_causality_permutation_test(X1.ch[i,], X2.ch[i,], H=H[i,], k=param.k, R=param.R, constrained_fit=TRUE, threshold_q=thresh_q)
    Pmc_constr_lgpdctc.ch[i] <- res_constr_lgpdctc.ch[[i]]$Pmc
  }
  
  if(!is.null(save_path)){
    save(Pmc_ctc.c,Pmc_ctc.ch,Pmc_pfcorr_lgpdctc.c,Pmc_pfcorr_lgpdctc.ch,Pmc_constr_lgpdctc.c,Pmc_constr_lgpdctc.ch,
         #res_ctc.c,res_ctc.ch,res_pfcorr_lgpdctc.c,res_pfcorr_lgpdctc.ch,res_constr_lgpdctc.c,res_constr_lgpdctc.ch,
         file=save_path)
  }
  tibble(beta21=beta21, Pmc_ctc.c=Pmc_ctc.c, Pmc_ctc.ch=Pmc_ctc.ch,
         Pmc_pfcorr_lgpdctc.c=Pmc_pfcorr_lgpdctc.c, Pmc_pfcorr_lgpdctc.ch=Pmc_pfcorr_lgpdctc.ch,
         Pmc_constr_lgpdctc.c=Pmc_constr_lgpdctc.c, Pmc_constr_lgpdctc.ch=Pmc_constr_lgpdctc.ch)
}

if(parallel_loop){
  stopCluster(cl)
}

#write_csv(res_tbl,"CTCs_simulations/Permutation_tests/CTCs_Permutation_test_power_last.csv")

end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)
#sink()

#rm(list = setdiff(ls(), union(lsf.str(), "sims_round")))


#Sys.sleep(300)
#system('shutdown -t 30 -s')

