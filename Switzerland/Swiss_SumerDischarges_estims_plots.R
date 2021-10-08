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

source("../R_functions/CTC_contribution_functions_OCP.R")

set.seed(1)

start_time <- Sys.time()

## This is all summer data as a tibble object
load("Data/river_dat.Rdata")
Disch <- river_dat %>% arrange(Date)
rm(river_dat)

## This is the meta-info for all stations, some scraped from the admin DAFU website
info_disch <- read_csv("data_wrangled/Discharges_info_merged.csv",guess_max = 10000)

Precipitations <- read_csv("Data/precip.csv",guess_max = 10000)
names(Precipitations)[1] <- "Date"

save_path <- "Results/causal_indep_pairs/"#"Results/" or NULL
Wplt <- 220#200 300
Hplt <- floor(0.4 * Wplt)

#Plot Colors
col_dif <- rgb(83/255,81/255,84/255)# Greys: "grey30" rgb(128/255,133/255,133/255) rgb(83/255,81/255,84/255)
col_v12 <- "#69b3a2"# Greens: "#69b3a2" rgb(132/255,186/255,91/255) rgb(62/255,150/255,81/255)
col_v21 <- rgb(57/255,106/255,177/255)# Blues: "#404080" rgb(114/255,147/255,203/255) rgb(57/255,106/255,177/255)

## Distribution of discharges among stations:
p1 <- select(Disch, -Date) %>% gather(key="Station", value="Value", factor_key=TRUE) %>% drop_na() %>% ggplot(aes(x=Value)) +
  geom_density(aes(color=Station), alpha=0.8, size=1 ) +
  coord_cartesian(xlim=c(0,50), ylim=c(0, 3)) + # zooms in
  labs(title="Density plot", subtitle="Danube river discharges at the different stations", x="Discharge")#,caption="Source: mpg"#,fill="# Cylinders"
#plot(p1)


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
print(as.data.frame(EVresults))

cat("\n=================== GPD FIT OF 95% EXCESSES (By Elevation) ===================\n")
elevations <- select(info_disch, id, Elevation)
EVresults_elev <- bind_cols(EVresults, elevations) %>% arrange(desc(Elevation))
print(as.data.frame(EVresults_elev))

par(mfrow=c(1,2))
plot(EVresults$scale, EVresults$shape, main="GPD scale-shape relationship for exteme discharges")
plot(EVresults_elev$Elevation, EVresults_elev$shape, main="GPD shape by station elevation for extreme discharges")
par(mfrow=c(1,1))

plot_scale <- EVresults_elev %>% ggplot(aes(x=scale, y=shape)) +
  geom_errorbar(aes(ymin=shape-se_shape, ymax=shape+se_shape), color="grey") +#, position=position_dodge(0.05), width=2
  geom_errorbar(aes(xmin=scale-se_scale, xmax=scale+se_scale), color="grey") +#darkgrey
  geom_point(size=-1) + geom_text(aes(label=id), nudge_x = 0, nudge_y = 0) +
  labs(title="GPD scale-shape relationship for exteme discharges", x=bquote("Estimated scale parameter ["*m^3*s^-1*"]"), y="Estimated shape parameter")
plot_elev <- EVresults_elev %>% ggplot(aes(x=Elevation, y=shape)) +
  geom_errorbar(aes(ymin=shape-se_shape, ymax=shape+se_shape), color="grey") +
  geom_point(size=-1) + geom_text(aes(label=id), nudge_x = 0, nudge_y = 0) +
  labs(title="GPD shape by station elevation for extreme discharges", x=bquote("Elevation [m]"), y="Estimated shape parameter")
pshape <- ggarrange(plot_scale, plot_elev, labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("Summer_shapes_ids.png", plot=pshape, device="png", path=save_path, width=Wplt, height=Wplt/2, units="mm", dpi=300)
plot(pshape)

cat("\n=================== GPD FIT OF 95% EXCESSES (By Average Volume) ===================\n")
ave_vols <- select(info_disch, id, Ave)
EVresults_vol <- bind_cols(EVresults, ave_vols) %>% arrange(desc(Ave))
print(as.data.frame(EVresults_vol))

plot_vol <- EVresults_vol %>% ggplot(aes(x=Ave, y=shape)) +
  geom_errorbar(aes(ymin=shape-se_shape, ymax=shape+se_shape), color="grey") +
  geom_point(size=-1) + geom_text(aes(label=id), nudge_x = 0, nudge_y = 0) +
  labs(title="GPD shape by station average discharges", x=bquote("Average Discharge ["*m^3*s^-1*"]"), y="Estimated shape parameter")
ggsave("Summer_shapeVol_ids.png", plot=plot_vol, device="png", path=save_path, width=Wplt/2, height=Wplt/2, units="mm", dpi=300)
plot(plot_vol)

plot_scale_log <- plot_scale + scale_x_log10()# + scale_x_continuous(trans="log")
plot_vol_log <- plot_vol + scale_x_log10()# + scale_x_continuous(trans="log")
pshape_log <- ggarrange(plot_scale_log, plot_elev, plot_vol_log, labels="AUTO", ncol=3, nrow=1)
ggsave("Summer_shapes_all_log.png", plot=pshape_log, device="png", path=save_path, width=3*Wplt/2, height=Wplt/2, units="mm", dpi=300)
plot(pshape_log)


## Causal orders and matrices of the stations (for upper tail only and both tails)
cat("\n=================== CAUSAL ORDER RETRIEVAL (SPACE) ===================\n")
common_Disch <- drop_na(select(Disch, -Date))
cat("\nLength of common discharges series:",nrow(common_Disch),"\n")

# print("Causal matrix (upper tails):")
# print(causal_tail_matrix(common_Disch, both_tails = FALSE))
# print("Causal matrix (both tails):")
# print(causal_tail_matrix(common_Disch, both_tails = TRUE))

cat("\nCausal Order estimation using upper tails only:\n")
causalorder_uppertail <- ease(common_Disch, both_tails = FALSE)
#print(names(common_Disch)[causalorder_uppertail])
print(causalorder_uppertail)

#cat("\nCausal Order estimation using both tails:")
#causalorder_bothtails <- ease(common_Disch, both_tails = TRUE)
#print(names(common_Disch)[causalorder_bothtails])




## ====================================== DISCHARGES BOOTSTRAPPING PARAMETERS AND DATA CONSTRUCTION ======================================
n_sim <- 1000#n_sim
R_test <- 10000
thresh_q <- 0.9#0.95
causal_pair <- c("station_42","station_63")
precip_causal <- c("ENT","FLU","EIT","LUZ")
indep_pair <- c("station_30","station_45")
precip_indep <- names(select(Precipitations, -Date))
mult_k <- 2#1

causal_str <- paste(substring(causal_pair[1],9,11),substring(causal_pair[2],9,11), sep="-")
indep_str <- paste(substring(indep_pair[1],9,11),substring(indep_pair[2],9,11), sep="-")
name_plot_str <- paste("c",causal_str,"-loc_i",indep_str, "-all", sep="")#-upstr -loc -all -same

size_max <- nrow(Disch)
k_max <- mult_k*floor(size_max^0.4)
Causal_data <- drop_na(Disch[causal_pair])
Indep_data <- drop_na(Disch[indep_pair])
size_causal <- nrow(Causal_data)
size_indep <- nrow(Indep_data)
k_causal <- mult_k*floor(size_causal^0.4)
k_indep <- mult_k*floor(size_indep^0.4)

Precip_H <- Precipitations %>% mutate(H_causal=rowMeans(select(., all_of(precip_causal))),#na.rm=TRUE
                                      H_indep=rowMeans(select(., all_of(precip_indep)))) %>% select(Date,H_causal,H_indep)
Disch_Precip <-  full_join(Disch,Precip_H,by="Date") %>% arrange(Date)
Causal_data_h <- drop_na(Disch_Precip[c(causal_pair,"H_causal")])
Indep_data_h <- drop_na(Disch_Precip[c(indep_pair,"H_indep")])
size_causal_h <- nrow(Causal_data_h)
size_indep_h <- nrow(Indep_data_h)
k_causal_h <- mult_k*floor(size_causal_h^0.4)
k_indep_h <- mult_k*floor(size_indep_h^0.4)

# Shapes of variables
EVprecipc <- tibble(scale = double(), se_scale = double(), shape = double(), se_shape = double(), Convergence = character())
precip_vectc <- drop_na(Precip_H["H_causal"])[["H_causal"]]
u<-quantile(precip_vectc,0.95)
EVresultc <- fpot(precip_vectc, threshold=u, control=ctrl)#, model="pp", method="Nelder-Mead", std.err = FALSE)
EVprecipc <- EVprecipc %>% add_row(scale=EVresultc$estimate[1], se_scale=EVresultc$std.err[1],
                                   shape=EVresultc$estimate[2], se_shape=EVresultc$std.err[2], Convergence=EVresultc$convergence)
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
cat("\nCausal pair:\n", sep="")
print(as.data.frame(EVall %>% filter(Station==causal_pair[1] | Station==causal_pair[2])))
cat("With Precipitation:\n", sep="")
print(as.data.frame(EVprecipc))



end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)
sink()

