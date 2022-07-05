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
source("../R_functions/CTC_contribution_functions_OCP.R")
source("competitors/competitor_helpers.R")

Swiss_SummerDischarges_competitors_routine <- function(station_pair=c("station_42","station_63"), precip_stations=c("ENT","FLU","EIT","LUZ"), precip_keyw="-loc",
                                                       alpha_pc=0.01, seed=NULL, save_path="Results/new_pairs/",
                                                       pair_type=c("unknown","causal","non-causal","independent"), ...){
  #-upstr -loc -all -same
  #source("../R_functions/CTC_contribution_functions_OCP.R")
  
  pair_type <- match.arg(pair_type)
  Pair_type <- str_to_sentence(pair_type)#str_to_title str_to_sentence
  
  check_directory(save_path, recursive=TRUE)
  
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
  if(precip_stations=="all"){precip_stations <- names(select(Precipitations, -Date))}
  
  stations_str <- paste(substring(station_pair[1],9,11),substring(station_pair[2],9,11), sep="-")
  caus_init <- if(pair_type=="causal"){"c"}else if(pair_type=="unknown"){"u"}else{"i"}
  name_plot_str <- paste(caus_init, stations_str, precip_keyw, sep="")#-upstr -loc -all -same
  
  n_max <- nrow(Disch)
  Stations_data <- drop_na(Disch[station_pair])
  n_obs <- nrow(Stations_data)
  
  Precip_H <- Precipitations %>% mutate(H_av=rowMeans(select(., all_of(precip_stations)))) %>% select(Date,H_av)#, na.rm=TRUE
  Disch_Precip <-  full_join(Disch,Precip_H,by="Date") %>% arrange(Date)
  Data_h <- drop_na(Disch_Precip[c(station_pair,"H_av")])
  n_obs_h <- nrow(Data_h)
  
  
  # Shapes of variables
  EVprecip <- tibble(scale = double(), se_scale = double(), shape = double(), se_shape = double(), Convergence = character())
  precip_vect <- drop_na(Precip_H["H_av"])[["H_av"]]
  u<-quantile(precip_vect,0.95)
  EVresultp <- fpot(precip_vect, threshold=u, control=ctrl)#, model="pp", method="Nelder-Mead", std.err = FALSE)
  EVprecip <- EVprecip %>% add_row(scale=EVresultp$estimate[1], se_scale=EVresultp$std.err[1],
                                   shape=EVresultp$estimate[2], se_shape=EVresultp$std.err[2], Convergence=EVresultp$convergence)
  
  elevations <- elevations %>% arrange(id)
  ave_vols <- ave_vols %>% arrange(id)
  EVall <- bind_cols(EVresults, select(elevations,-id), select(ave_vols,-id))
  
  sink(paste(save_path,name_plot_str,".txt", sep=""))
  cat("\n=================== Stations Estimated GPD Parameters ===================\n")
  cat("\n", Pair_type, " pair:\n", sep="")
  print(as.data.frame(EVall %>% filter(Station==station_pair[1] | Station==station_pair[2])))
  cat("With Precipitation:\n", sep="")
  print(as.data.frame(EVprecip))
  
  
  cat("\n")
  cat("\n=================== LiNGAM ===================\n")
  
  # 
  cat("\nLiNGAM for ", pair_type, " pair ",paste(station_pair)," only (n=",n_obs,"):\n", sep="")
  
  lingam_p <- LiNGAM(Stations_data, labels=c("X1", "X2"), verbose=FALSE)
  
  lingam_p_coef12 <- lingam_p$Bpruned[2,1]
  lingam_p_caus12 <- lingam_p$adj[1,2]
  
  cat("LiNGAM (1,2): B1->2 = ", lingam_p_coef12, "\n", sep="")
  print(lingam_p$Bpruned)
  # print(lingam_p$adj)
  print(lingam_p$effects)
  
  
  cat("\nLiNGAM with precipitation for ", name_plot_str," (n=",n_obs_h,"):\n", sep="")
  
  lingam_h <- LiNGAM(Data_h, labels=c("X1", "X2", "H"), verbose=FALSE)
  
  lingam_h_coef12 <- lingam_h$Bpruned[2,1]
  lingam_h_caus12 <- lingam_h$adj[1,2]
  
  cat("LiNGAM (1,2,H): B1->2 = ", lingam_h_coef12, "\n", sep="")
  print(lingam_h$Bpruned)
  # print(lingam_h$adj)
  print(lingam_h$effects)
  
  
  
  cat("\n")
  cat("\n=================== PC Algorithm ===================\n")
  
  # 
  cat("\nPC Alg. for ", pair_type, " pair ",paste(station_pair)," only (n=",n_obs,"):\n", sep="")
  
  pc_p <- pc_algorithm(Stations_data, rank_version=FALSE, alpha=alpha_pc, labels=c("X1", "X2"),
                       skel.method="stable.fast", ...)#"stable" "stable.fast"
  
  pc_p_caus12 <- causal_cpdag_helper(pc_p$cpdag, 1, 2)
  
  cat("PC Alg. (1,2): 1-?->2 = ", pc_p_caus12, "\n", sep="")
  print(pc_p$cpdag)
  #print(pc_p$effects)
  print(pc_p$fit)
  
  
  cat("\nPC Alg. with precipitation for ", name_plot_str," (n=",n_obs_h,"):\n", sep="")
  
  pc_h <- pc_algorithm(Data_h, rank_version=FALSE, alpha=alpha_pc, labels=c("X1", "X2", "H"),
                       skel.method="stable.fast", ...)#"stable" "stable.fast"
  
  pc_h_caus12 <- causal_cpdag_helper(pc_h$cpdag, 1, 2)
  
  cat("PC Alg. (1,2,H): 1-?->2 = ", pc_h_caus12, "\n", sep="")
  print(pc_h$cpdag)
  #print(pc_h$effects)
  print(pc_h$fit)
  
  
  cat("\n")
  cat("\n=================== PC (rank) Algorithm ===================\n")
  
  # 
  cat("\nPC rank Alg. for ", pair_type, " pair ",paste(station_pair)," only (n=",n_obs,"):\n", sep="")
  
  pcr_p <- pc_algorithm(Stations_data, rank_version=TRUE, alpha=alpha_pc, labels=c("X1", "X2"),
                        skel.method="stable.fast", ...)#"stable" "stable.fast"
  
  pcr_p_caus12 <- causal_cpdag_helper(pcr_p$cpdag, 1, 2)
  
  cat("PC rank Alg. (1,2): 1-?->2 = ", pcr_p_caus12, "\n", sep="")
  print(pcr_p$cpdag)
  #print(pcr_p$effects)
  print(pcr_p$fit)
  
  
  cat("\nPC rank Alg. with precipitation for ", name_plot_str," (n=",n_obs_h,"):\n", sep="")
  
  pcr_h <- pc_algorithm(Data_h, rank_version=TRUE, alpha=alpha_pc, labels=c("X1", "X2", "H"),
                        skel.method="stable.fast", ...)#"stable" "stable.fast"
  
  pcr_h_caus12 <- causal_cpdag_helper(pcr_h$cpdag, 1, 2)
  
  cat("PC rank Alg. (1,2,H): 1-?->2 = ", pcr_h_caus12, "\n", sep="")
  print(pcr_h$cpdag)
  #print(pcr_h$effects)
  print(pcr_h$fit)
  
  
  
  output <- c(lingam_p_coef12=lingam_p_coef12, lingam_h_coef12=lingam_h_coef12,
              #lingam_p_caus12=lingam_p_caus12, lingam_h_caus12=lingam_h_caus12,
              pc_p_caus12=pc_p_caus12, pc_h_caus12=pc_h_caus12,
              pcr_p_caus12=pcr_p_caus12, pcr_h_caus12=pcr_h_caus12,
              alpha_pc=alpha_pc, n_obs=n_obs, n_obs_h=n_obs_h)
  
  Results <- list(output=output,
                  lingam_p=lingam_p, lingam_h=lingam_h, pc_p=pc_p, pc_h=pc_h, pcr_p=pcr_p, pcr_h=pcr_h)
  safe_save_rds(Results, paste0(save_path, "RESULTS_", name_plot_str, ".rds"))
  
  end_time <- Sys.time()
  cat("\nRun time:\n")
  print(end_time - start_time)
  sink()
  return(Results)
}

## Proxies for labelling plots and outputs nicely differenciating causal/non-causal stations

Swiss_SummerDischarges_competitors_causal <- function(causal_pair=c("station_42","station_63"), precip_causal=c("ENT","FLU","EIT","LUZ"), precip_keyw="-loc",
                                                      alpha_pc=0.01, seed=NULL, save_path="Results/causal_pairs/"){
  #-upstr -loc -all -same
  
  Results <- Swiss_SummerDischarges_competitors_routine(station_pair=causal_pair, precip_stations=precip_causal, precip_keyw=precip_keyw,
                                                       alpha_pc=alpha_pc, seed=seed, save_path=save_path, pair_type="causal")
  # output <- out_in
  # output$size_causal <- out_in$n_obs
  # output$size_causal_h <- out_in$n_obs_h
  return(Results)
}


Swiss_SummerDischarges_competitors_indep <- function(indep_pair=c("station_30","station_45"), precip_indep="all", precip_keyw="-all",
                                                     alpha_pc=0.01, seed=NULL, save_path="Results/indep_pairs/"){
  #-upstr -loc -all -same
  
  Results <- Swiss_SummerDischarges_competitors_routine(station_pair=indep_pair, precip_stations=precip_indep, precip_keyw=precip_keyw,
                                                       alpha_pc=alpha_pc, seed=seed, save_path=save_path, pair_type="non-causal")
  # output <- out_in
  # output$size_indep <- out_in$n_obs
  # output$size_indep_h <- out_in$n_obs_h
  return(Results)
}






