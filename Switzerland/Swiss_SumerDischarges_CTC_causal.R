library(tidyverse)
library(tsibble)
library(feasts)
library(fable)
library("evd")
library("evdbayes")
library("ismev")
library(latex2exp)
library(causalXtreme)
library(boot)
library("ggpubr")
#theme_set(theme_bw() + theme(legend.position = "top"))

Swiss_SummerDischarges_CTCs_causal <- function(causal_pair=c("station_42","station_63"), precip_causal=c("ENT","FLU","EIT","LUZ"), precip_keyw="-loc",
                                               n_bootstrap_sims=0, R_test=10000, thresh_q=0.9, mult_k=2, seed=NULL, save_path="Results/causal_pairs/"){
  #-upstr -loc -all -same
  #source("../R_functions/CTC_contribution_functions_OCP.R")
  
  if(!is.null(seed)){set.seed(seed)}
  
  start_time <- Sys.time()
  
  ## This is all summer data as a tibble object
  load("Data/river_dat.Rdata")
  Disch <- river_dat %>% arrange(Date)
  rm(river_dat)
  
  ## This is the meta-info for all stations, some scraped from the admin DAFU website
  info_disch <- read_csv("data_wrangled/Discharges_info_merged.csv",guess_max = 10000)
  
  Precipitations <- read_csv("Data/precip.csv",guess_max = 10000)
  names(Precipitations)[1] <- "Date"
  
  Wplt <- 220#200 300
  Hplt <- floor(0.4 * Wplt)
  
  #Plot Colors
  col_dif <- rgb(83/255,81/255,84/255)# Greys: "grey30" rgb(128/255,133/255,133/255) rgb(83/255,81/255,84/255)
  col_v12 <- "#69b3a2"# Greens: "#69b3a2" rgb(132/255,186/255,91/255) rgb(62/255,150/255,81/255)
  col_v21 <- rgb(57/255,106/255,177/255)# Blues: "#404080" rgb(114/255,147/255,203/255) rgb(57/255,106/255,177/255)
  
  
  ## Estimating point process of threshold excesses (to compare shape parameter)
  cat("\n=================== GPD FIT OF 95% EXCESSES (All sations) ===================\n")
  ctrl <- list(maxit = 100000)#names(ctrl) <- c("maxit")
  EVresults <- tibble(Station = character(), scale = double(), se_scale = double(), shape = double(), se_shape = double(), Convergence = character())
  for (stati in names(select(Disch, -Date))){
    disch_vect <- drop_na(Disch[stati])[[stati]]
    u<-quantile(disch_vect,0.95)
    EVresult <- fpot(disch_vect, threshold=u, control=ctrl)#, model="pp", method="Nelder-Mead", std.err = FALSE)
    EVresults <- EVresults %>% add_row(Station=stati, scale=EVresult$estimate[1], se_scale=EVresult$std.err[1],
                                       shape=EVresult$estimate[2], se_shape=EVresult$std.err[2], Convergence=EVresult$convergence)
  }
  rm(disch_vect)
  rm(EVresult)
  #print(as.data.frame(EVresults))
  
  #cat("\n=================== GPD FIT OF 95% EXCESSES (By Elevation) ===================\n")
  elevations <- select(info_disch, id, Elevation)
  EVresults_elev <- bind_cols(EVresults, elevations) %>% arrange(desc(Elevation))
  
  #cat("\n=================== GPD FIT OF 95% EXCESSES (By Average Volume) ===================\n")
  ave_vols <- select(info_disch, id, Ave)
  EVresults_vol <- bind_cols(EVresults, ave_vols) %>% arrange(desc(Ave))
  
  ## Causal orders and matrices of the stations (for upper tail only and both tails)
  #cat("\n=================== CAUSAL ORDER RETRIEVAL (SPACE) ===================\n")
  common_Disch <- drop_na(select(Disch, -Date))
  #cat("\nLength of common discharges series:",nrow(common_Disch),"\n")
  
  
  
  ## ====================================== DISCHARGES BOOTSTRAPPING PARAMETERS AND DATA CONSTRUCTION ======================================
  if(precip_causal=="all"){precip_causal <- names(select(Precipitations, -Date))}
  
  causal_str <- paste(substring(causal_pair[1],9,11),substring(causal_pair[2],9,11), sep="-")
  name_plot_str <- paste("c",causal_str, precip_keyw, sep="")#-upstr -loc -all -same
  
  size_max <- nrow(Disch)
  k_max <- mult_k*floor(size_max^0.4)
  Causal_data <- drop_na(Disch[causal_pair])
  size_causal <- nrow(Causal_data)
  k_causal <- mult_k*floor(size_causal^0.4)
  
  Precip_H <- Precipitations %>% mutate(H_causal=rowMeans(select(., all_of(precip_causal))),#na.rm=TRUE
  ) %>% select(Date,H_causal)
  Disch_Precip <-  full_join(Disch,Precip_H,by="Date") %>% arrange(Date)
  Causal_data_h <- drop_na(Disch_Precip[c(causal_pair,"H_causal")])
  size_causal_h <- nrow(Causal_data_h)
  k_causal_h <- mult_k*floor(size_causal_h^0.4)
  
  
  # Shapes of variables
  EVprecipc <- tibble(scale = double(), se_scale = double(), shape = double(), se_shape = double(), Convergence = character())
  precip_vectc <- drop_na(Precip_H["H_causal"])[["H_causal"]]
  u<-quantile(precip_vectc,0.95)
  EVresultc <- fpot(precip_vectc, threshold=u, control=ctrl)#, model="pp", method="Nelder-Mead", std.err = FALSE)
  EVprecipc <- EVprecipc %>% add_row(scale=EVresultc$estimate[1], se_scale=EVresultc$std.err[1],
                                     shape=EVresultc$estimate[2], se_shape=EVresultc$std.err[2], Convergence=EVresultc$convergence)
  
  elevations <- elevations %>% arrange(id)
  ave_vols <- ave_vols %>% arrange(id)
  EVall <- bind_cols(EVresults, select(elevations,-id), select(ave_vols,-id))
  
  sink(paste(save_path,name_plot_str,".txt", sep=""))
  cat("\n=================== Stations Estimated GPD Parameters ===================\n")
  cat("\nCausal pair:\n", sep="")
  print(as.data.frame(EVall %>% filter(Station==causal_pair[1] | Station==causal_pair[2])))
  cat("With Precipitation:\n", sep="")
  print(as.data.frame(EVprecipc))
  
  
  cat("\n")
  cat("\n=================== Bootstrapping Discharge data (NP CTC) ===================\n")
  
  # bootstrapping with R replications for tail coefficients
  cat("\nNon-parametric Causal tail coeffs for causal pair ",paste(causal_pair)," (n=",size_causal,", k=",k_causal,"):\n", sep="")
  if(n_bootstrap_sims>0){
  bs_causal <- boot(data=Causal_data, statistic=CTC_boot_stats, R=n_bootstrap_sims, k=k_causal, both_tails=FALSE)
  # view results
  print(bs_causal$t0)
  #plot(bs_causal, index=1)
  #plot(bs_causal, index=2)
  cat("\n")
  # get 95% confidence intervals
  print_CI(boot.ci(bs_causal, conf = 0.95, type=c("basic","perc"), index=1), "CTC(1,2)")#type=c("basic","bca","perc")
  print_CI(boot.ci(bs_causal, conf = 0.95, type=c("basic","perc"), index=2), "CTC(2,1)")
  print_CI(boot.ci(bs_causal, conf = 0.95, type=c("basic","perc"), index=3), "CTC(1,2) - CTC(2,1)")
  cat("\n")
  
  plot_xlims <- c(min(bs_causal$t[,3], na.rm=TRUE),max(bs_causal$t[,3], na.rm=TRUE))
  hcausal <- ggplot(data.frame(diff12=bs_causal$t[,3]), aes(x=x) ) +
    geom_histogram( aes(x = diff12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="\"Causal\" pair of stations", x="Non-parametric CTC difference", y="Density") +
    expand_limits(x=plot_xlims)
  #plot(ggarrange(hindep, hcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  
  plot_xlims <- c(min(bs_causal$t[,1:2], na.rm=TRUE), 1)
  pcausal <- ggplot(data.frame(Gamma12=bs_causal$t[,1], Gamma21=bs_causal$t[,2]), aes(x=x) ) +
    geom_histogram( aes(x = Gamma12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12 ) +
    geom_histogram( aes(x = Gamma21, y = -..density..), binwidth=0.005,boundary=1,closed="right", fill= col_v21) + labs(title="\"Causal\" pair of stations", x="Non-parametric causal tail coefficients", y="Density") +
    expand_limits(x=plot_xlims)
  #plot(ggarrange(pindep, pcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  pstati <- ggarrange(hcausal, pcausal, labels = c("A", "B"), ncol = 1, nrow = 2)
  ggsave(paste(name_plot_str, "_npctc.png", sep=""), plot=pstati, device="png", path=save_path, width=Wplt/2, height=Hplt, units="mm", dpi=300)
  }
  
  causal_test <- CTC_causality_permutation_test(Causal_data[[1]], Causal_data[[2]], H=NULL, k=k_causal, R=R_test,
                                                constrained_fit=FALSE, threshold_q=thresh_q)
  Pmc_npctc <- causal_test$Pmc
  cat("Monte Carlo Permutation test (R = ",R_test,"): Pmc = ",Pmc_npctc,"\n\n", sep="")
  
  
  
  
  cat("\n")
  cat("\n=================== Bootstrapping Discharge data (unconstrained LGPD CTC) ===================\n")
  
  # bootstrapping with R replications for tail coefficients
  cat("\nLGPD Causal tail coeffs for causal pair ",paste(causal_pair)," (n=",size_causal_h,", k=",k_causal_h,"):\n", sep="")
  if(n_bootstrap_sims>0){
  bs_causal_h <- boot(data=Causal_data_h, statistic=LGPDCTC_boot_stats, R=n_bootstrap_sims, k=k_causal_h, constrained_fit=FALSE, threshold_q=thresh_q, parametric_F1=TRUE, both_tails=FALSE)
  # view results
  print(bs_causal_h$t0)
  #plot(bs_causal_h, index=1)
  #plot(bs_causal_h, index=2)
  cat("\n")
  # get 95% confidence intervals
  print_CI(boot.ci(bs_causal_h, conf = 0.95, type=c("basic","perc"), index=1), "LGPD CTC(1,2)")#type=c("basic","bca","perc")
  print_CI(boot.ci(bs_causal_h, conf = 0.95, type=c("basic","perc"), index=2), "LGPD CTC(2,1)")
  print_CI(boot.ci(bs_causal_h, conf = 0.95, type=c("basic","perc"), index=3), "LGPD CTC(1,2) - LGPD CTC(2,1)")
  cat("\n")
  
  plot_xlims <- c(min(bs_causal_h$t[,3], na.rm=TRUE),max(bs_causal_h$t[,3], na.rm=TRUE))
  hcausal <- ggplot(data.frame(diff12=bs_causal_h$t[,3]), aes(x=x) ) +
    geom_histogram( aes(x = diff12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="\"Causal\" pair of stations", x="Post-fit corrected LGPD CTC difference", y="Density") +
    expand_limits(x=plot_xlims)
  #plot(ggarrange(hindep, hcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  
  plot_xlims <- c(min(bs_causal_h$t[,1:2], na.rm=TRUE), 1)
  pcausal <- ggplot(data.frame(Gamma12=bs_causal_h$t[,1], Gamma21=bs_causal_h$t[,2]), aes(x=x) ) +
    geom_histogram( aes(x = Gamma12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram( aes(x = Gamma21, y = -..density..), binwidth=0.005,boundary=1,closed="right", fill= col_v21) + labs(title="\"Causal\" pair of stations", x="Post-fit corrected LGPD CTCs", y="Density") +
    expand_limits(x=plot_xlims)
  #plot(ggarrange(pindep, pcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  pstati <- ggarrange(hcausal, pcausal, labels = c("A", "B"), ncol = 1, nrow = 2)
  ggsave(paste(name_plot_str, "_pfcorr_lgpdctc.png", sep=""), plot=pstati, device="png", path=save_path, width=Wplt/2, height=Hplt, units="mm", dpi=300)
  }
  
  MLEs_unconstr <- LGPD_causal_tail_coeff(Causal_data_h[[1]], Causal_data_h[[2]], H=Causal_data_h[[3]], k=k_causal_h,
                                          threshold_q=thresh_q, return_MLEs=TRUE)
  
  causal_test <- CTC_causality_permutation_test(Causal_data_h[[1]], Causal_data_h[[2]], H=Causal_data_h[[3]], k=k_causal_h, R=R_test,
                                                constrained_fit=FALSE, threshold_q=thresh_q)
  Pmc_pfcorr_lgpdctc <- causal_test$Pmc
  cat("Monte Carlo Permutation test (R = ",R_test,"): Pmc = ",Pmc_pfcorr_lgpdctc,"\n\n", sep="")
  
  
  cat("\n")
  cat("\n=================== Bootstrapping Discharge data (constrained LGPD CTC) ===================\n")
  
  # bootstrapping with R replications for tail coefficients
  cat("\nLGPD Causal tail coeffs for causal pair ",paste(causal_pair)," (n=",size_causal_h,", k=",k_causal_h,"):\n", sep="")
  if(n_bootstrap_sims>0){
  bs_causal_h <- boot(data=Causal_data_h, statistic=LGPDCTC_boot_stats, R=n_bootstrap_sims, k=k_causal_h, constrained_fit=TRUE, threshold_q=thresh_q, parametric_F1=TRUE, both_tails=FALSE)
  # view results
  print(bs_causal_h$t0)
  #plot(bs_causal_h, index=1)
  #plot(bs_causal_h, index=2)
  cat("\n")
  # get 95% confidence intervals
  print_CI(boot.ci(bs_causal_h, conf = 0.95, type=c("basic","perc"), index=1), "LGPD CTC(1,2)")#type=c("basic","bca","perc")
  print_CI(boot.ci(bs_causal_h, conf = 0.95, type=c("basic","perc"), index=2), "LGPD CTC(2,1)")
  print_CI(boot.ci(bs_causal_h, conf = 0.95, type=c("basic","perc"), index=3), "LGPD CTC(1,2) - LGPD CTC(2,1)")
  cat("\n")
  
  plot_xlims <- c(min(bs_causal_h$t[,3], na.rm=TRUE),max(bs_causal_h$t[,3], na.rm=TRUE))
  hcausal <- ggplot(data.frame(diff12=bs_causal_h$t[,3]), aes(x=x) ) +
    geom_histogram( aes(x = diff12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="\"Causal\" pair of stations", x="Constrained fit LGPD CTC difference", y="Density") +
    expand_limits(x=plot_xlims)
  #plot(ggarrange(hindep, hcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  
  plot_xlims <- c(min(bs_causal_h$t[,1:2], na.rm=TRUE), 1)
  pcausal <- ggplot(data.frame(Gamma12=bs_causal_h$t[,1], Gamma21=bs_causal_h$t[,2]), aes(x=x) ) +
    geom_histogram( aes(x = Gamma12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram( aes(x = Gamma21, y = -..density..), binwidth=0.005,boundary=1,closed="right", fill= col_v21) + labs(title="\"Causal\" pair of stations", x="Constrained fit LGPD CTCs", y="Density") +
    expand_limits(x=plot_xlims)
  #plot(ggarrange(pindep, pcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  pstati <- ggarrange(hcausal, pcausal, labels = c("A", "B"), ncol = 1, nrow = 2)
  ggsave(paste(name_plot_str, "_constr_lgpdctc.png", sep=""), plot=pstati, device="png", path=save_path, width=Wplt/2, height=Hplt, units="mm", dpi=300)
  }
  
  MLEs_constr <- constrained_LGPD_CTC(Causal_data_h[[1]], Causal_data_h[[2]], H=Causal_data_h[[3]], k=k_causal_h,
                                      threshold_q=thresh_q, return_MLEs=TRUE)
  
  causal_test <- CTC_causality_permutation_test(Causal_data_h[[1]], Causal_data_h[[2]], H=Causal_data_h[[3]], k=k_causal_h, R=R_test,
                                                constrained_fit=TRUE, threshold_q=thresh_q)
  Pmc_constr_lgpdctc <- causal_test$Pmc
  cat("Monte Carlo Permutation test (R = ",R_test,"): Pmc = ",Pmc_constr_lgpdctc,"\n\n", sep="")
  
  
  did_exp <- FALSE
  
  try({
  cat("\n")
  cat("\n=================== Bootstrapping Discharge data (exp GPD CTC) ===================\n")
  
  # bootstrapping with R replications for tail coefficients
  cat("\nGPD Causal tail coeffs for causal pair ",paste(causal_pair)," (n=",size_causal_h,", k=",k_causal_h,"):\n", sep="")
  causal_test <- expGPD_CTC_causality_permutation_test(Causal_data_h[[1]], Causal_data_h[[2]], H=Causal_data_h[[3]], k=k_causal_h, R=R_test,
                                                       threshold_q=thresh_q)
  Pmc_exp_gpdctc <- causal_test$Pmc
  did_exp <- TRUE
  
  if(n_bootstrap_sims>0){
  bs_causal_h <- boot(data=Causal_data_h, statistic=expGPDCTC_boot_stats, R=n_bootstrap_sims, k=k_causal_h, threshold_q=thresh_q, parametric_F1=TRUE, both_tails=FALSE)
  # view results
  print(bs_causal_h$t0)
  #plot(bs_causal_h, index=1)
  #plot(bs_causal_h, index=2)
  cat("\n")
  # get 95% confidence intervals
  print_CI(boot.ci(bs_causal_h, conf = 0.95, type=c("basic","perc"), index=1), "GPD CTC(1,2)")#type=c("basic","bca","perc")
  print_CI(boot.ci(bs_causal_h, conf = 0.95, type=c("basic","perc"), index=2), "GPD CTC(2,1)")
  print_CI(boot.ci(bs_causal_h, conf = 0.95, type=c("basic","perc"), index=3), "GPD CTC(1,2) - GPD CTC(2,1)")
  cat("\n")
  
  plot_xlims <- c(min(bs_causal_h$t[,3], na.rm=TRUE),max(bs_causal_h$t[,3], na.rm=TRUE))
  hcausal <- ggplot(data.frame(diff12=bs_causal_h$t[,3]), aes(x=x) ) +
    geom_histogram( aes(x = diff12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="\"Causal\" pair of stations", x="Exponential GPD CTC difference", y="Density") +
    expand_limits(x=plot_xlims)
  #plot(ggarrange(hindep, hcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  
  plot_xlims <- c(min(bs_causal_h$t[,1:2], na.rm=TRUE), 1)
  pcausal <- ggplot(data.frame(Gamma12=bs_causal_h$t[,1], Gamma21=bs_causal_h$t[,2]), aes(x=x) ) +
    geom_histogram( aes(x = Gamma12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram( aes(x = Gamma21, y = -..density..), binwidth=0.005,boundary=1,closed="right", fill= col_v21) + labs(title="\"Causal\" pair of stations", x="Exponential GPD CTCs", y="Density") +
    expand_limits(x=plot_xlims)
  #plot(ggarrange(pindep, pcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  pstati <- ggarrange(hcausal, pcausal, labels = c("A", "B"), ncol = 1, nrow = 2)
  ggsave(paste(name_plot_str, "_exp_lgpdctc.png", sep=""), plot=pstati, device="png", path=save_path, width=Wplt/2, height=Hplt, units="mm", dpi=300)
  }
  cat("Monte Carlo Permutation test (R = ",R_test,"): Pmc = ",Pmc_exp_gpdctc,"\n", sep="")
  })
  
  output <- list(Pmc_npctc=Pmc_npctc, Pmc_pfcorr_lgpdctc=Pmc_pfcorr_lgpdctc, Pmc_constr_lgpdctc=Pmc_constr_lgpdctc)
  if(did_exp){
    output$Pmc_exp_gpdctc <- Pmc_exp_gpdctc
  } else {
    output$Pmc_exp_gpdctc <- NaN
  }
  output <- c(output, list(MLEs_unconstr=MLEs_unconstr, MLEs_constr=MLEs_constr, shape_H=EVprecipc$shape, se_shape_H=EVprecipc$se_shape,
                           size_causal=size_causal, k_causal=k_causal, size_causal_h=size_causal_h, k_causal_h=k_causal_h))
  
  end_time <- Sys.time()
  cat("\nRun time:\n")
  print(end_time - start_time)
  sink()
  return(output)
}

