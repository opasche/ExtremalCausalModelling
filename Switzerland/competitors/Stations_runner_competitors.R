library(parallel)
library(doParallel)

library(tidyverse)
library(tsibble)
library(feasts)
library(fable)
library(evd)
library(ismev)
library(latex2exp)
library(causalXtreme)
library(boot)
library(ggpubr)

source("../R_functions/CTC_contribution_functions_OCP.R")
source("competitors/Swiss_SumerDischarges_competitors_routine.R")

start_time_all <- Sys.time()

nb_cores_total <- detectCores()#12
nb_cores_used <- nb_cores_total-1
cl <- makeCluster(nb_cores_used)
registerDoParallel(cl)

alpha_pc <- 0.01 #0.01 5e-2
seed <- 1
save_path_causal <- "Results/competitors/causal_pairs/test/"# "Results/competitors/causal_pairs/"
save_path_indep <- "Results/competitors/indep_pairs/test/"# "Results/competitors/indep_pairs/"
check_directory(save_path_causal, recursive=TRUE)
check_directory(save_path_indep, recursive=TRUE)

params_list_causal <- list(
  list(causal_pair = c("station_43","station_62"),
       precip_causal = c("INT","LTB","BRZ","MER"),
       precip_kw = "-loc"),
  
  list(causal_pair = c("station_42","station_63"),
       precip_causal = c("ENT","FLU","EIT","LUZ"),
       precip_kw = "-loc"),
  
  list(causal_pair = c("station_36","station_63"),
       precip_causal = c("ENT","FLU","EIT","LUZ"),
       precip_kw = "-loc"),
  
  list(causal_pair = c("station_24","station_61"),
       precip_causal = c("ZWE","BOT","WIS","THU"),
       precip_kw = "-loc"),
  
  list(causal_pair = c("station_44","station_61"),
       precip_causal = c("FRU","KIE","ABO","KAS","THU"),#,"WIS"
       precip_kw = "-loc"),
  
  list(causal_pair = c("station_22","station_38"),
       precip_causal = c("LAG"),
       precip_kw = "-loc"),
  
  list(causal_pair = c("station_22","station_35"),
       precip_causal = c("LAG","WIE","AIE","BUD","WIE","AIE","BUD"),
       precip_kw = "-loc")
)


params_list_indep <- list(
  list(indep_pair = c("station_30","station_45"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_36","station_39"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_42","station_34"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_42","station_34"),
       precip_indep = c("LUZ","EIT","SRN","ENG"), precip_kw = "-loc"),
  
  list(indep_pair = c("station_32","station_33"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_62","station_63"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_60","station_57"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_13","station_14"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_17","station_22"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_12","station_21"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_26","station_28"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_27","station_31"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_23","station_39"),
       precip_indep = "all", precip_kw = "-all"),
  
  list(indep_pair = c("station_23","station_35"),
       precip_indep = "all", precip_kw = "-all")
)


packages_to_pass <- c("tidyverse","tsibble","feasts","fable","evd","ismev","latex2exp","causalXtreme","boot","ggpubr")
variables_to_pass <- c("alpha_pc", "seed", "save_path_causal", "save_path_indep")


results_causal <- foreach(params=params_list_causal, .packages=packages_to_pass, .export=variables_to_pass,
                          .errorhandling="stop", .combine=rbind) %dopar% {
  source("../R_functions/CTC_contribution_functions_OCP.R")
  source("competitors/Swiss_SumerDischarges_competitors_routine.R")
  Results <- Swiss_SummerDischarges_competitors_causal(causal_pair=params$causal_pair, precip_causal=params$precip_causal, precip_keyw=params$precip_kw,
                                                       alpha_pc=alpha_pc, seed=seed, save_path=save_path_causal)
  Results$output
}
#save(results_causal, params_list_causal, file="Results/results_Station_runner_last.RData")

results_indep <- foreach(params=params_list_indep, .packages=packages_to_pass, .export=variables_to_pass,
                         .errorhandling="stop", .combine=rbind) %dopar% {
  source("../R_functions/CTC_contribution_functions_OCP.R")
  source("competitors/Swiss_SumerDischarges_competitors_routine.R")
  Results <- Swiss_SummerDischarges_competitors_indep(indep_pair=params$indep_pair, precip_indep=params$precip_indep, precip_keyw=params$precip_kw,
                                                      alpha_pc=alpha_pc, seed=seed, save_path=save_path_indep)
  Results$output
}
#save(results_causal, results_indep, params_list_causal, params_list_indep, file="Results/results_Station_runner_last.RData")

stopCluster(cl)

results_causal_tibble0 <- tibble::tibble(
  id=1:length(params_list_causal),
  station_pair=sapply(params_list_causal, function(x) (paste(substring(x$causal_pair[1],9,11),substring(x$causal_pair[2],9,11), sep="-"))),
  precipitation=sapply(params_list_causal, function(x) (paste(x$precip_causal, sep="-", collapse = "-"))),
  precip_kw=sapply(params_list_causal, function(x) (x$precip_kw))
)

results_indep_tibble0 <- tibble::tibble(
  id=1:length(params_list_indep),
  station_pair=sapply(params_list_indep, function(x) (paste(substring(x$indep_pair[1],9,11),substring(x$indep_pair[2],9,11), sep="-"))),
  precipitation=sapply(params_list_indep, function(x) (paste(x$precip_indep, sep="-", collapse = "-"))),
  precip_kw=sapply(params_list_indep, function(x) (x$precip_kw))
)

results_causal_tibble <- bind_cols(results_causal_tibble0, tibble::as_tibble(results_causal))
results_indep_tibble <- bind_cols(results_indep_tibble0, tibble::as_tibble(results_indep))

#save(results_causal_tibble, results_indep_tibble, file="Results/results_Station_runner_last_tibbles.RData")

sink("Results_station_runner.txt")

cat("\nResults causal stations:\n")
print(as.data.frame(results_causal_tibble))

cat("\nResults indep stations:\n")
print(as.data.frame(results_indep_tibble))

end_time_results <- Sys.time()
cat("\nRun time (Results):\n")
print(end_time_results - start_time_all)

sink()

results_causal_tibble <- bind_cols(results_causal_tibble0, pair_type="causal", tibble::as_tibble(results_causal))
results_indep_tibble <- bind_cols(results_indep_tibble0, pair_type="non-causal", tibble::as_tibble(results_indep))
results_CH_competitors_all <- bind_rows(results_causal_tibble,results_indep_tibble)

filename <- paste0("Results/competitors/results_Station_runner_",format(Sys.time(),'%Y%m%d_%H%M%S'),".csv")
write_csv(results_CH_competitors_all, file=filename)


end_time_all <- Sys.time()
cat("\nRun time:\n")
print(end_time_all - start_time_all)


#Sys.sleep(300)
#system('shutdown -t 30 -s')
