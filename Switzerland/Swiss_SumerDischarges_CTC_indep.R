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

Swiss_SummerDischarges_CTCs_indep <- function(indep_pair=c("station_30","station_45"), precip_indep="all", precip_keyw="-all",
                                              n_bootstrap_sims=0, R_test=10000, thresh_q=0.9, mult_k=2, seed=NULL, save_path="Results/indep_pairs/"){
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
  if(precip_indep=="all"){precip_indep <- names(select(Precipitations, -Date))}
  
  indep_str <- paste(substring(indep_pair[1],9,11),substring(indep_pair[2],9,11), sep="-")
  name_plot_str <- paste("i",indep_str, precip_keyw, sep="")#-upstr -loc -all -same
  
  size_max <- nrow(Disch)
  k_max <- mult_k*floor(size_max^0.4)
  Indep_data <- drop_na(Disch[indep_pair])
  size_indep <- nrow(Indep_data)
  k_indep <- mult_k*floor(size_indep^0.4)
  
  Precip_H <- Precipitations %>% mutate(H_indep=rowMeans(select(., all_of(precip_indep)))) %>% select(Date,H_indep)
  Disch_Precip <-  full_join(Disch,Precip_H,by="Date") %>% arrange(Date)
  Indep_data_h <- drop_na(Disch_Precip[c(indep_pair,"H_indep")])
  size_indep_h <- nrow(Indep_data_h)
  k_indep_h <- mult_k*floor(size_indep_h^0.4)
  
  
  # Shapes of variables
  EVprecipi <- tibble(scale = double(), se_scale = double(), shape = double(), se_shape = double(), Convergence = character())
  precip_vecti <- drop_na(Precip_H["H_indep"])[["H_indep"]]
  u<-quantile(precip_vecti,0.95)
  EVresulti <- fpot(precip_vecti, threshold=u, control=ctrl)#, model="pp", method="Nelder-Mead", std.err = FALSE)
  EVprecipi <- EVprecipi %>% add_row(scale=EVresulti$estimate[1], se_scale=EVresulti$std.err[1],
                                     shape=EVresulti$estimate[2], se_shape=EVresulti$std.err[2], Convergence=EVresulti$convergence)
  
  elevations <- elevations %>% arrange(id)
  ave_vols <- ave_vols %>% arrange(id)
  EVall <- bind_cols(EVresults, select(elevations,-id), select(ave_vols,-id))
  
  sink(paste(save_path,name_plot_str,".txt", sep=""))
  cat("\n=================== Stations Estimated GPD Parameters ===================\n")
  cat("\nIndependent pair:\n", sep="")
  print(as.data.frame(EVall %>% filter(Station==indep_pair[1] | Station==indep_pair[2])))
  cat("With Precipitation:\n", sep="")
  print(as.data.frame(EVprecipi))
  
  
  cat("\n")
  cat("\n=================== Bootstrapping Discharge data (NP CTC) ===================\n")
  
  # bootstrapping with R replications for tail coefficients
  cat("\nNon-parametric Causal tail coeffs for \"independent\" pair ",paste(indep_pair)," (n=",size_indep,", k=",k_indep,"):\n", sep="")
  if(n_bootstrap_sims>0){
  bs_indep <- boot(data=Indep_data, statistic=CTC_boot_stats, R=n_bootstrap_sims, k=k_indep, both_tails=FALSE)
  # view results
  print(bs_indep$t0)
  #plot(bs_indep, index=1)
  #plot(bs_indep, index=2)
  cat("\n")
  # get 95% confidence intervals
  print_CI(boot.ci(bs_indep, conf = 0.95, type=c("basic","perc"), index=1), "CTC(1,2)")#type=c("basic","bca","perc")
  print_CI(boot.ci(bs_indep, conf = 0.95, type=c("basic","perc"), index=2), "CTC(2,1)")
  print_CI(boot.ci(bs_indep, conf = 0.95, type=c("basic","perc"), index=3), "CTC(1,2) - CTC(2,1)")
  cat("\n")
  
  
  plot_xlims <- c(min(bs_indep$t[,3], na.rm=TRUE),max(bs_indep$t[,3], na.rm=TRUE))
  hindep <- ggplot(data.frame(diff12=bs_indep$t[,3]), aes(x=x) ) +
    geom_histogram( aes(x = diff12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="\"Independent\" stations", x="Non-parametric CTC difference", y="Density") +
    expand_limits(x=plot_xlims)#xlim(plot_xlims) scale_x_continuous(limits=plot_xlims) coord_cartesian(xlim=plot_xlims) expand_limits(x=plot_xlims)
  #plot(ggarrange(hindep, hcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  
  plot_xlims <- c(min(bs_indep$t[,1:2], na.rm=TRUE), 1)
  pindep <- ggplot(data.frame(Gamma12=bs_indep$t[,1], Gamma21=bs_indep$t[,2]), aes(x=x) ) +
    geom_histogram( aes(x = Gamma12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12 ) +
    geom_histogram( aes(x = Gamma21, y = -..density..), binwidth=0.005,boundary=1,closed="right", fill= col_v21) + labs(title="\"Independent\" stations", x="Non-parametric causal tail coefficients", y="Density") +
    expand_limits(x=plot_xlims)#xlim(plot_xlims) scale_x_continuous(limits=plot_xlims) coord_cartesian(xlim=plot_xlims) expand_limits(x=plot_xlims)
  #plot(ggarrange(pindep, pcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  pstati <- ggarrange(hindep, pindep, labels = c("A", "B"), ncol = 1, nrow = 2)
  ggsave(paste(name_plot_str, "_npctc.png", sep=""), plot=pstati, device="png", path=save_path, width=Wplt/2, height=Hplt, units="mm", dpi=300)
  }
  
  indep_test <- CTC_causality_permutation_test(Indep_data[[1]], Indep_data[[2]], H=NULL, k=k_indep, R=R_test,
                                               constrained_fit=FALSE, threshold_q=thresh_q)
  Pmc_npctc <- indep_test$Pmc
  cat("Monte Carlo Permutation test (R = ",R_test,"): Pmc = ",Pmc_npctc,"\n\n", sep="")
  
  
  
  
  cat("\n")
  cat("\n=================== Bootstrapping Discharge data (unconstrained LGPD CTC) ===================\n")
  
  # bootstrapping with R replications for tail coefficients
  cat("\nLGPD Causal tail coeffs for \"independent\" pair ",paste(indep_pair)," (n=",size_indep_h,", k=",k_indep_h,"):\n", sep="")
  if(n_bootstrap_sims>0){
  bs_indep_h <- boot(data=Indep_data_h, statistic=LGPDCTC_boot_stats, R=n_bootstrap_sims, k=k_indep_h, constrained_fit=FALSE, threshold_q=thresh_q, parametric_F1=TRUE, both_tails=FALSE)
  # view results
  print(bs_indep_h$t0)
  #plot(bs_indep_h, index=1)
  #plot(bs_indep_h, index=2)
  cat("\n")
  # get 95% confidence intervals
  print_CI(boot.ci(bs_indep_h, conf = 0.95, type=c("basic","perc"), index=1), "LGPD CTC(1,2)")#type=c("basic","bca","perc")
  print_CI(boot.ci(bs_indep_h, conf = 0.95, type=c("basic","perc"), index=2), "LGPD CTC(2,1)")
  print_CI(boot.ci(bs_indep_h, conf = 0.95, type=c("basic","perc"), index=3), "LGPD CTC(1,2) - LGPD CTC(2,1)")
  cat("\n")
  
  
  plot_xlims <- c(min(bs_indep_h$t[,3], na.rm=TRUE),max(bs_indep_h$t[,3], na.rm=TRUE))
  hindep <- ggplot(data.frame(diff12=bs_indep_h$t[,3]), aes(x=x) ) +
    geom_histogram( aes(x = diff12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="\"Independent\" stations", x="Post-fit corrected LGPD CTC difference", y="Density") +
    expand_limits(x=plot_xlims)#xlim(plot_xlims) scale_x_continuous(limits=plot_xlims) coord_cartesian(xlim=plot_xlims) expand_limits(x=plot_xlims)
  #plot(ggarrange(hindep, hcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  
  plot_xlims <- c(min(bs_indep_h$t[,1:2], na.rm=TRUE), 1)
  pindep <- ggplot(data.frame(Gamma12=bs_indep_h$t[,1], Gamma21=bs_indep_h$t[,2]), aes(x=x) ) +
    geom_histogram( aes(x = Gamma12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram( aes(x = Gamma21, y = -..density..), binwidth=0.005,boundary=1,closed="right", fill= col_v21) + labs(title="\"Independent\" stations", x="Post-fit corrected LGPD CTCs", y="Density") +
    expand_limits(x=plot_xlims)#xlim(plot_xlims) scale_x_continuous(limits=plot_xlims) coord_cartesian(xlim=plot_xlims) expand_limits(x=plot_xlims)
  #plot(ggarrange(pindep, pcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  pstati <- ggarrange(hindep, pindep, labels = c("A", "B"), ncol = 1, nrow = 2)
  ggsave(paste(name_plot_str, "_pfcorr_lgpdctc.png", sep=""), plot=pstati, device="png", path=save_path, width=Wplt/2, height=Hplt, units="mm", dpi=300)
  }
  
  MLEs_unconstr <- LGPD_causal_tail_coeff(Indep_data_h[[1]], Indep_data_h[[2]], H=Indep_data_h[[3]], k=k_indep_h,
                                          threshold_q=thresh_q, return_MLEs=TRUE)
  
  indep_test <- CTC_causality_permutation_test(Indep_data_h[[1]], Indep_data_h[[2]], H=Indep_data_h[[3]], k=k_indep_h, R=R_test,
                                               constrained_fit=FALSE, threshold_q=thresh_q)
  Pmc_pfcorr_lgpdctc <- indep_test$Pmc
  cat("Monte Carlo Permutation test (R = ",R_test,"): Pmc = ",Pmc_pfcorr_lgpdctc,"\n\n", sep="")
  
  
  cat("\n")
  cat("\n=================== Bootstrapping Discharge data (constrained LGPD CTC) ===================\n")
  
  # bootstrapping with R replications for tail coefficients
  cat("\nLGPD Causal tail coeffs for \"independent\" pair ",paste(indep_pair)," (n=",size_indep_h,", k=",k_indep_h,"):\n", sep="")
  if(n_bootstrap_sims>0){
  bs_indep_h <- boot(data=Indep_data_h, statistic=LGPDCTC_boot_stats, R=n_bootstrap_sims, k=k_indep_h, constrained_fit=TRUE, threshold_q=thresh_q, parametric_F1=TRUE, both_tails=FALSE)
  # view results
  print(bs_indep_h$t0)
  #plot(bs_indep_h, index=1)
  #plot(bs_indep_h, index=2)
  cat("\n")
  # get 95% confidence intervals
  print_CI(boot.ci(bs_indep_h, conf = 0.95, type=c("basic","perc"), index=1), "LGPD CTC(1,2)")#type=c("basic","bca","perc")
  print_CI(boot.ci(bs_indep_h, conf = 0.95, type=c("basic","perc"), index=2), "LGPD CTC(2,1)")
  print_CI(boot.ci(bs_indep_h, conf = 0.95, type=c("basic","perc"), index=3), "LGPD CTC(1,2) - LGPD CTC(2,1)")
  cat("\n")
  
  
  plot_xlims <- c(min(bs_indep_h$t[,3], na.rm=TRUE),max(bs_indep_h$t[,3], na.rm=TRUE))
  hindep <- ggplot(data.frame(diff12=bs_indep_h$t[,3]), aes(x=x) ) +
    geom_histogram( aes(x = diff12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="\"Independent\" stations", x="Constrained fit LGPD CTC difference", y="Density") +
    expand_limits(x=plot_xlims)#xlim(plot_xlims) scale_x_continuous(limits=plot_xlims) coord_cartesian(xlim=plot_xlims) expand_limits(x=plot_xlims)
  #plot(ggarrange(hindep, hcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  
  plot_xlims <- c(min(bs_indep_h$t[,1:2], na.rm=TRUE), 1)
  pindep <- ggplot(data.frame(Gamma12=bs_indep_h$t[,1], Gamma21=bs_indep_h$t[,2]), aes(x=x) ) +
    geom_histogram( aes(x = Gamma12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram( aes(x = Gamma21, y = -..density..), binwidth=0.005,boundary=1,closed="right", fill= col_v21) + labs(title="\"Independent\" stations", x="Constrained fit LGPD CTCs", y="Density") +
    expand_limits(x=plot_xlims)#xlim(plot_xlims) scale_x_continuous(limits=plot_xlims) coord_cartesian(xlim=plot_xlims) expand_limits(x=plot_xlims)
  #plot(ggarrange(pindep, pcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
  pstati <- ggarrange(hindep, pindep, labels = c("A", "B"), ncol = 1, nrow = 2)
  ggsave(paste(name_plot_str, "_constr_lgpdctc.png", sep=""), plot=pstati, device="png", path=save_path, width=Wplt/2, height=Hplt, units="mm", dpi=300)
  }
  
  MLEs_constr <- constrained_LGPD_CTC(Indep_data_h[[1]], Indep_data_h[[2]], H=Indep_data_h[[3]], k=k_indep_h,
                                      threshold_q=thresh_q, return_MLEs=TRUE)
  
  indep_test <- CTC_causality_permutation_test(Indep_data_h[[1]], Indep_data_h[[2]], H=Indep_data_h[[3]], k=k_indep_h, R=R_test,
                                               constrained_fit=TRUE, threshold_q=thresh_q)
  Pmc_constr_lgpdctc <- indep_test$Pmc
  cat("Monte Carlo Permutation test (R = ",R_test,"): Pmc = ",Pmc_constr_lgpdctc,"\n\n", sep="")
  
  
  did_exp <- FALSE
  
  try({
    cat("\n")
    cat("\n=================== Bootstrapping Discharge data (exp GPD CTC) ===================\n")
    
    # bootstrapping with R replications for tail coefficients
    cat("\nGPD Causal tail coeffs for \"independent\" pair ",paste(indep_pair)," (n=",size_indep_h,", k=",k_indep_h,"):\n", sep="")
    indep_test <- expGPD_CTC_causality_permutation_test(Indep_data_h[[1]], Indep_data_h[[2]], H=Indep_data_h[[3]], k=k_indep_h, R=R_test,
                                                        threshold_q=thresh_q)
    Pmc_exp_gpdctc <- indep_test$Pmc
    did_exp <- TRUE
    
    if(n_bootstrap_sims>0){
    bs_indep_h <- boot(data=Indep_data_h, statistic=expGPDCTC_boot_stats, R=n_bootstrap_sims, k=k_indep_h, threshold_q=thresh_q, parametric_F1=TRUE, both_tails=FALSE)
    # view results
    print(bs_indep_h$t0)
    #plot(bs_indep_h, index=1)
    #plot(bs_indep_h, index=2)
    cat("\n")
    # get 95% confidence intervals
    print_CI(boot.ci(bs_indep_h, conf = 0.95, type=c("basic","perc"), index=1), "GPD CTC(1,2)")#type=c("basic","bca","perc")
    print_CI(boot.ci(bs_indep_h, conf = 0.95, type=c("basic","perc"), index=2), "GPD CTC(2,1)")
    print_CI(boot.ci(bs_indep_h, conf = 0.95, type=c("basic","perc"), index=3), "GPD CTC(1,2) - GPD CTC(2,1)")
    cat("\n")
    
    
    plot_xlims <- c(min(bs_indep_h$t[,3], na.rm=TRUE),max(bs_indep_h$t[,3], na.rm=TRUE))
    hindep <- ggplot(data.frame(diff12=bs_indep_h$t[,3]), aes(x=x) ) +
      geom_histogram( aes(x = diff12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="\"Independent\" stations", x="Exponential GPD CTC difference", y="Density") +
      expand_limits(x=plot_xlims)#xlim(plot_xlims) scale_x_continuous(limits=plot_xlims) coord_cartesian(xlim=plot_xlims) expand_limits(x=plot_xlims)
    #plot(ggarrange(hindep, hcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
    
    plot_xlims <- c(min(bs_indep_h$t[,1:2], na.rm=TRUE), 1)
    pindep <- ggplot(data.frame(Gamma12=bs_indep_h$t[,1], Gamma21=bs_indep_h$t[,2]), aes(x=x) ) +
      geom_histogram( aes(x = Gamma12, y = ..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
      geom_histogram( aes(x = Gamma21, y = -..density..), binwidth=0.005,boundary=1,closed="right", fill= col_v21) + labs(title="\"Independent\" stations", x="Exponential GPD CTCs", y="Density") +
      expand_limits(x=plot_xlims)#xlim(plot_xlims) scale_x_continuous(limits=plot_xlims) coord_cartesian(xlim=plot_xlims) expand_limits(x=plot_xlims)
    #plot(ggarrange(pindep, pcausal, labels = c("A", "B"), ncol = 2, nrow = 1))
    pstati <- ggarrange(hindep, pindep, labels = c("A", "B"), ncol = 1, nrow = 2)
    ggsave(paste(name_plot_str, "_exp_lgpdctc.png", sep=""), plot=pstati, device="png", path=save_path, width=Wplt/2, height=Hplt, units="mm", dpi=300)
    }
    cat("Monte Carlo Permutation test (R = ",R_test,"): Pmc = ",Pmc_exp_gpdctc,"\n\n", sep="")
  })
  
  output <- list(Pmc_npctc=Pmc_npctc, Pmc_pfcorr_lgpdctc=Pmc_pfcorr_lgpdctc, Pmc_constr_lgpdctc=Pmc_constr_lgpdctc)
  if(did_exp){
    output$Pmc_exp_gpdctc <- Pmc_exp_gpdctc
  } else {
    output$Pmc_exp_gpdctc <- NaN
  }
  output <- c(output, list(MLEs_unconstr=MLEs_unconstr, MLEs_constr=MLEs_constr, shape_H=EVprecipi$shape, se_shape_H=EVprecipi$se_shape,
                           size_indep=size_indep, k_indep=k_indep, size_indep_h=size_indep_h, k_indep_h=k_indep_h))
  
  end_time <- Sys.time()
  cat("\nRun time:\n")
  print(end_time - start_time)
  sink()
  return(output)
}

