library(tidyverse)
#library("evd")
#library("evdbayes")
#library("ismev")
library(latex2exp)
#library(causalXtreme)
#library(EnvStats)
#library(boot)
#library("plyr") library(reshape) library(ggplot2) library(tidyr) library(dplyr)
library("ggpubr")
#theme_set(theme_bw() + theme(legend.position = "top"))

source("../R_functions/CTC_contribution_functions_OCP.R")

save_path <- "Results/CTCs_simulations/q90_1M_2/"#"Results/CTCs_simulations/" or NULL
sims_path <- "CTCs_simulations/q90_1M/"
check_directory(dirname(save_path), recursive=TRUE)

#files_list <- list.files(path="CTCs_simulations/q90_1M/", pattern="CTCs_*_1M2k1000.RData", full.names=FALSE)
general_params_code <- "_q90_1M2k1000"#"_q90_100k2k1000""_q90_10k2k1000""_1M2k1000"
distr_codes <- c("ln",
                 "ln_H0p5",
                 "ln_H1p5",
                 "ln_H2",
                 "Pa1-2",
                 "Pa1-2_H1-1",
                 "Pa2_H1p5",
                 "Pa2_H1p2",
                 "Pa1-2_H1-4",
                 "Pa2p3_H1p8",
                 "Pa3_H1p5",
                 "t4",
                 "t4_H2",
                 "t4_H3",
                 "t4_H6")
codes_list <- paste0(distr_codes, general_params_code)

#Plots widths and heights for the 2 types of plots
Wdif <- 200
Hdif <- floor(1/3 * Wdif)
Wvals <- 200
Hvals <- floor(1/2 * Wvals)

#Plot Colors
col_dif <- rgb(83/255,81/255,84/255)# Greys: "grey30" rgb(128/255,133/255,133/255) rgb(83/255,81/255,84/255)
col_v12 <- "#69b3a2"# Greens: "#69b3a2" rgb(132/255,186/255,91/255) rgb(62/255,150/255,81/255)
col_v21 <- rgb(57/255,106/255,177/255)# Blues: "#404080" rgb(114/255,147/255,203/255) rgb(57/255,106/255,177/255)

set.seed(1)

start_time <- Sys.time()

for(sims_code in codes_list){
  cat("\n================================", sims_code, "================================\n")
  
  load(paste0(sims_path, "CTCs_", sims_code, ".RData"))
  Sys.sleep(1)
  
  diffctc.i <- ctc12.i - ctc21.i
  diffctc.ih <- ctc12.ih - ctc21.ih
  diffctc.c <- ctc12.c - ctc21.c
  diffctc.ch <- ctc12.ch - ctc21.ch
  diffgpdctc.i <- gpdctc12.i - gpdctc21.i
  diffgpdctc.ih <- gpdctc12.ih - gpdctc21.ih
  diffgpdctc.c <- gpdctc12.c - gpdctc21.c
  diffgpdctc.ch <- gpdctc12.ch - gpdctc21.ch
  diffpfcorr_lgpdctc.i <- pfcorr_lgpdctc12.i - pfcorr_lgpdctc21.i
  diffpfcorr_lgpdctc.ih <- pfcorr_lgpdctc12.ih - pfcorr_lgpdctc21.ih
  diffpfcorr_lgpdctc.c <- pfcorr_lgpdctc12.c - pfcorr_lgpdctc21.c
  diffpfcorr_lgpdctc.ch <- pfcorr_lgpdctc12.ch - pfcorr_lgpdctc21.ch
  diffconstr_lgpdctc.i <- constr_lgpdctc12.i - constr_lgpdctc21.i
  diffconstr_lgpdctc.ih <- constr_lgpdctc12.ih - constr_lgpdctc21.ih
  diffconstr_lgpdctc.c <- constr_lgpdctc12.c - constr_lgpdctc21.c
  diffconstr_lgpdctc.ch <- constr_lgpdctc12.ch - constr_lgpdctc21.ch
  
  ## ================ NP CTC ================
  cat("\n================ NP CTC ================\n")
  
  cat("NAs: \tA(",sum(is.na(ctc12.i)),",",sum(is.na(ctc21.i)),") B(",sum(is.na(ctc12.ih)),",",sum(is.na(ctc21.ih)),")
    \tC(",sum(is.na(ctc12.c)),",",sum(is.na(ctc21.c)),") D(",sum(is.na(ctc12.ch)),",",sum(is.na(ctc21.ch)),")\n", sep="")
  
  plot_xlimsi <- c(min(c(diffctc.i,diffctc.ih,-0.1),na.rm=TRUE), max(c(diffctc.i,diffctc.ih,0.1),na.rm=TRUE))
  plot_xlimsc <- c(min(c(diffctc.c,diffctc.ch,0),na.rm=TRUE), max(c(diffctc.c,diffctc.ch,0.3),na.rm=TRUE))
  
  
  d1 <- ggplot(data.frame(diff12=diffctc.i), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="Independent variables", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_vline(xintercept=0, col="black", size=1)
  d2 <- ggplot(data.frame(diff12=diffctc.ih), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="Common hidden confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_vline(xintercept=0, col="black", size=1)
  d3 <- ggplot(data.frame(diff12=diffctc.c), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="X1 -> X2", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_vline(xintercept=0.25, col="black", size=1)
  d4 <- ggplot(data.frame(diff12=diffctc.ch), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="X1 -> X2 with confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_vline(xintercept=1/8, col="black", size=1)
  pdif_ctc <- ggarrange(d1, d2, d3, d4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2) %>%
    annotate_figure(left=text_grob("Density", rot = 90), bottom="Non-parametric CTC difference")
  #plot(pdif_ctc)
  Sys.sleep(0.5)
  ggsave(paste(sims_code, "_dif_npctc.png", sep=""), plot=pdif_ctc, device="png", path=save_path, width=Wdif, height=Hdif, units="mm", dpi=300)
  Sys.sleep(0.5)
  
  plot_xlimsi <- c(min(c(ctc12.i,ctc21.i,ctc12.ih,ctc21.ih, 0.75),na.rm=TRUE), 1)
  plot_xlimsc <- c(min(c(ctc12.c,ctc21.c,ctc12.ch,ctc21.ch, 0.75),na.rm=TRUE), 1)
  pci <- ggplot(data.frame(Gamma12=ctc12.i, Gamma21=ctc21.i), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="Independent variables", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_segment(x=0.5, xend=0.5, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.5, xend=0.5, y=0, yend=-Inf, col=col_v21, size=1)
  pcih <- ggplot(data.frame(Gamma12=ctc12.ih, Gamma21=ctc21.ih), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="Common hidden confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_segment(x=0.75, xend=0.75, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.75, xend=0.75, y=0, yend=-Inf, col=col_v21, size=1)
  pcc <- ggplot(data.frame(Gamma12=ctc12.c, Gamma21=ctc21.c), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="X1 -> X2", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_segment(x=1, xend=1, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.75, xend=0.75, y=0, yend=-Inf, col=col_v21, size=1)
  pcch <- ggplot(data.frame(Gamma12=ctc12.ch, Gamma21=ctc21.ch), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="X1 -> X2 with confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_segment(x=1, xend=1, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=7/8, xend=7/8, y=0, yend=-Inf, col=col_v21, size=1)
  pvals_ctc <- ggarrange(pci,pcih,pcc,pcch, labels=c("A", "B", "C", "D"), ncol=2, nrow=2) %>%
    annotate_figure(left=text_grob("Density", rot = 90), bottom="Non-parametric causal tail coefficients")
  #plot(pvals_ctc)
  Sys.sleep(0.5)
  ggsave(paste(sims_code, "_vals_npctc.png", sep=""), plot=pvals_ctc, device="png", path=save_path, width=Wvals, height=Hvals, units="mm", dpi=300)
  Sys.sleep(0.5)
  
  
  ## ================ GPD CTC ================
  cat("\n================ GPD CTC ================\n")
  
  cat("NAs: \tA(",sum(is.na(gpdctc12.i)),",",sum(is.na(gpdctc21.i)),") B(",sum(is.na(gpdctc12.ih)),",",sum(is.na(gpdctc21.ih)),")
    \tC(",sum(is.na(gpdctc12.c)),",",sum(is.na(gpdctc21.c)),") D(",sum(is.na(gpdctc12.ch)),",",sum(is.na(gpdctc21.ch)),")\n", sep="")
  
  plot_xlimsi <- c(min(c(diffgpdctc.i,diffgpdctc.ih,-0.1),na.rm=TRUE), max(c(diffgpdctc.i,diffgpdctc.ih,0.1),na.rm=TRUE))
  plot_xlimsc <- c(min(c(diffgpdctc.c,diffgpdctc.ch,0),na.rm=TRUE), max(c(diffgpdctc.c,diffgpdctc.ch,0.3),na.rm=TRUE))
  
  
  d1 <- ggplot(data.frame(diff12=diffgpdctc.i), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="Independent variables", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_vline(xintercept=0, col="black", size=1)
  d2 <- ggplot(data.frame(diff12=diffgpdctc.ih), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="Common hidden confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_vline(xintercept=0, col="black", size=1)
  d3 <- ggplot(data.frame(diff12=diffgpdctc.c), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="X1 -> X2", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_vline(xintercept=0.25, col="black", size=1)
  d4 <- ggplot(data.frame(diff12=diffgpdctc.ch), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="X1 -> X2 with confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_vline(xintercept=1/8, col="black", size=1)
  pdif_gpdctc <- ggarrange(d1, d2, d3, d4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2) %>%
    annotate_figure(left=text_grob("Density", rot = 90), bottom="GPD CTC difference")
  #plot(pdif_gpdctc)
  Sys.sleep(0.5)
  ggsave(paste(sims_code, "_dif_gpdctc.png", sep=""), plot=pdif_gpdctc, device="png", path=save_path, width=Wdif, height=Hdif, units="mm", dpi=300)
  Sys.sleep(0.5)
  
  plot_xlimsi <- c(min(c(gpdctc12.i,gpdctc21.i,gpdctc12.ih,gpdctc21.ih, 0.75),na.rm=TRUE), 1)
  plot_xlimsc <- c(min(c(gpdctc12.c,gpdctc21.c,gpdctc12.ch,gpdctc21.ch, 0.75),na.rm=TRUE), 1)
  pci <- ggplot(data.frame(Gamma12=gpdctc12.i, Gamma21=gpdctc21.i), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="Independent variables", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_segment(x=0.5, xend=0.5, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.5, xend=0.5, y=0, yend=-Inf, col=col_v21, size=1)
  pcih <- ggplot(data.frame(Gamma12=gpdctc12.ih, Gamma21=gpdctc21.ih), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="Common hidden confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_segment(x=0.75, xend=0.75, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.75, xend=0.75, y=0, yend=-Inf, col=col_v21, size=1)
  pcc <- ggplot(data.frame(Gamma12=gpdctc12.c, Gamma21=gpdctc21.c), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="X1 -> X2", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_segment(x=1, xend=1, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.75, xend=0.75, y=0, yend=-Inf, col=col_v21, size=1)
  pcch <- ggplot(data.frame(Gamma12=gpdctc12.ch, Gamma21=gpdctc21.ch), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="X1 -> X2 with confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_segment(x=1, xend=1, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=7/8, xend=7/8, y=0, yend=-Inf, col=col_v21, size=1)
  pvals_gpdctc <- ggarrange(pci,pcih,pcc,pcch, labels=c("A", "B", "C", "D"), ncol=2, nrow=2) %>%
    annotate_figure(left=text_grob("Density", rot = 90), bottom="GPD causal tail coefficients")
  #plot(pvals_gpdctc)
  Sys.sleep(0.5)
  ggsave(paste(sims_code, "_vals_gpdctc.png", sep=""), plot=pvals_gpdctc, device="png", path=save_path, width=Wvals, height=Hvals, units="mm", dpi=300)
  Sys.sleep(0.5)
  
  
  ## ================ Post-fit corr. LGPD CTC ================
  cat("\n================ Post-fit corr. LGPD CTC ================\n")
  
  cat("NAs: \tA(",sum(is.na(pfcorr_lgpdctc12.i)),",",sum(is.na(pfcorr_lgpdctc21.i)),") B(",sum(is.na(pfcorr_lgpdctc12.ih)),",",sum(is.na(pfcorr_lgpdctc21.ih)),")
    \tC(",sum(is.na(pfcorr_lgpdctc12.c)),",",sum(is.na(pfcorr_lgpdctc21.c)),") D(",sum(is.na(pfcorr_lgpdctc12.ch)),",",sum(is.na(pfcorr_lgpdctc21.ch)),")\n", sep="")
  
  plot_xlimsi <- c(min(c(diffpfcorr_lgpdctc.i,diffpfcorr_lgpdctc.ih,-0.1),na.rm=TRUE), max(c(diffpfcorr_lgpdctc.i,diffpfcorr_lgpdctc.ih,0.1),na.rm=TRUE))
  plot_xlimsc <- c(min(c(diffpfcorr_lgpdctc.c,diffpfcorr_lgpdctc.ch,0),na.rm=TRUE), max(c(diffpfcorr_lgpdctc.c,diffpfcorr_lgpdctc.ch,0.3),na.rm=TRUE))
  
  
  d1 <- ggplot(data.frame(diff12=diffpfcorr_lgpdctc.i), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="Independent variables", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_vline(xintercept=0, col="black", size=1)
  d2 <- ggplot(data.frame(diff12=diffpfcorr_lgpdctc.ih), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="Common hidden confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_vline(xintercept=0, col="black", size=1)
  d3 <- ggplot(data.frame(diff12=diffpfcorr_lgpdctc.c), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="X1 -> X2", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_vline(xintercept=0.25, col="black", size=1)
  d4 <- ggplot(data.frame(diff12=diffpfcorr_lgpdctc.ch), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="X1 -> X2 with confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_vline(xintercept=1/8, col="black", size=1)# + geom_vline(xintercept=0.25, col="black", size=1)
  pdif_pfclgpdctc <- ggarrange(d1, d2, d3, d4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2) %>%
    annotate_figure(left=text_grob("Density", rot = 90), bottom="Post-fit corrected LGPD CTC difference")
  #plot(pdif_pfclgpdctc)
  Sys.sleep(0.5)
  ggsave(paste(sims_code, "_dif_pfcorr_lgpdctc.png", sep=""), plot=pdif_pfclgpdctc, device="png", path=save_path, width=Wdif, height=Hdif, units="mm", dpi=300)
  Sys.sleep(0.5)
  
  
  plot_xlimsi <- c(min(c(pfcorr_lgpdctc12.i,pfcorr_lgpdctc21.i,pfcorr_lgpdctc12.ih,pfcorr_lgpdctc21.ih, 0.75),na.rm=TRUE), 1)
  plot_xlimsc <- c(min(c(pfcorr_lgpdctc12.c,pfcorr_lgpdctc21.c,pfcorr_lgpdctc12.ch,pfcorr_lgpdctc21.ch, 0.75),na.rm=TRUE), 1)
  pci <- ggplot(data.frame(Gamma12=pfcorr_lgpdctc12.i, Gamma21=pfcorr_lgpdctc21.i), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="Independent variables", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_segment(x=0.5, xend=0.5, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.5, xend=0.5, y=0, yend=-Inf, col=col_v21, size=1)
  pcih <- ggplot(data.frame(Gamma12=pfcorr_lgpdctc12.ih, Gamma21=pfcorr_lgpdctc21.ih), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="Common hidden confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_segment(x=0.75, xend=0.75, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.75, xend=0.75, y=0, yend=-Inf, col=col_v21, size=1)#+
  #geom_segment(x=0.5, xend=0.5, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.5, xend=0.5, y=0, yend=-Inf, col=col_v21, size=1)
  pcc <- ggplot(data.frame(Gamma12=pfcorr_lgpdctc12.c, Gamma21=pfcorr_lgpdctc21.c), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="X1 -> X2", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_segment(x=1, xend=1, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.75, xend=0.75, y=0, yend=-Inf, col=col_v21, size=1)
  pcch <- ggplot(data.frame(Gamma12=pfcorr_lgpdctc12.ch, Gamma21=pfcorr_lgpdctc21.ch), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="X1 -> X2 with confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_segment(x=1, xend=1, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=7/8, xend=7/8, y=0, yend=-Inf, col=col_v21, size=1)#+
  #geom_segment(x=0.75, xend=0.75, y=0, yend=-Inf, col=col_v21, size=1)
  pvals_pfclgpdctc <- ggarrange(pci,pcih,pcc,pcch, labels=c("A", "B", "C", "D"), ncol=2, nrow=2) %>%
    annotate_figure(left=text_grob("Density", rot = 90), bottom="Post-fit corrected LGPD CTCs")
  #plot(pvals_pfclgpdctc)
  Sys.sleep(0.5)
  ggsave(paste(sims_code, "_vals_pfcorr_lgpdctc.png", sep=""), plot=pvals_pfclgpdctc, device="png", path=save_path, width=Wvals, height=Hvals, units="mm", dpi=300)
  Sys.sleep(0.5)
  
  
  ## ================ Constrained fit LGPD CTC ================
  cat("\n================ Constrained fit LGPD CTC ================\n")
  
  cat("NAs: \tA(",sum(is.na(constr_lgpdctc12.i)),",",sum(is.na(constr_lgpdctc21.i)),") B(",sum(is.na(constr_lgpdctc12.ih)),",",sum(is.na(constr_lgpdctc21.ih)),")
    \tC(",sum(is.na(constr_lgpdctc12.c)),",",sum(is.na(constr_lgpdctc21.c)),") D(",sum(is.na(constr_lgpdctc12.ch)),",",sum(is.na(constr_lgpdctc21.ch)),")\n", sep="")
  
  plot_xlimsi <- c(min(c(diffconstr_lgpdctc.i,diffconstr_lgpdctc.ih,-0.1),na.rm=TRUE), max(c(diffconstr_lgpdctc.i,diffconstr_lgpdctc.ih,0.1),na.rm=TRUE))
  plot_xlimsc <- c(min(c(diffconstr_lgpdctc.c,diffconstr_lgpdctc.ch,0),na.rm=TRUE), max(c(diffconstr_lgpdctc.c,diffconstr_lgpdctc.ch,0.3),na.rm=TRUE))
  
  
  d1 <- ggplot(data.frame(diff12=diffconstr_lgpdctc.i), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="Independent variables", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_vline(xintercept=0, col="black", size=1)
  d2 <- ggplot(data.frame(diff12=diffconstr_lgpdctc.ih), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="Common hidden confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_vline(xintercept=0, col="black", size=1)
  d3 <- ggplot(data.frame(diff12=diffconstr_lgpdctc.c), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="X1 -> X2", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_vline(xintercept=0.25, col="black", size=1)
  d4 <- ggplot(data.frame(diff12=diffconstr_lgpdctc.ch), aes(x=x) ) +
    geom_histogram(aes(x=diff12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_dif) + labs(title="X1 -> X2 with confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_vline(xintercept=1/8, col="black", size=1)# + geom_vline(xintercept=0.25, col="black", size=1)
  pdif_constrlgpdctc <- ggarrange(d1, d2, d3, d4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2) %>%
    annotate_figure(left=text_grob("Density", rot = 90), bottom="Constrained fit LGPD CTC difference")
  #plot(pdif_constrlgpdctc)
  Sys.sleep(0.5)
  ggsave(paste(sims_code, "_dif_constr_lgpdctc.png", sep=""), plot=pdif_constrlgpdctc, device="png", path=save_path, width=Wdif, height=Hdif, units="mm", dpi=300)
  Sys.sleep(0.5)
  
  
  plot_xlimsi <- c(min(c(constr_lgpdctc12.i,constr_lgpdctc21.i,constr_lgpdctc12.ih,constr_lgpdctc21.ih, 0.75),na.rm=TRUE), 1)
  plot_xlimsc <- c(min(c(constr_lgpdctc12.c,constr_lgpdctc21.c,constr_lgpdctc12.ch,constr_lgpdctc21.ch, 0.75),na.rm=TRUE), 1)
  pci <- ggplot(data.frame(Gamma12=constr_lgpdctc12.i, Gamma21=constr_lgpdctc21.i), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="Independent variables", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_segment(x=0.5, xend=0.5, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.5, xend=0.5, y=0, yend=-Inf, col=col_v21, size=1)
  pcih <- ggplot(data.frame(Gamma12=constr_lgpdctc12.ih, Gamma21=constr_lgpdctc21.ih), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="Common hidden confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsi) + geom_segment(x=0.75, xend=0.75, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.75, xend=0.75, y=0, yend=-Inf, col=col_v21, size=1)#+
  #geom_segment(x=0.5, xend=0.5, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.5, xend=0.5, y=0, yend=-Inf, col=col_v21, size=1)
  pcc <- ggplot(data.frame(Gamma12=constr_lgpdctc12.c, Gamma21=constr_lgpdctc21.c), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="X1 -> X2", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_segment(x=1, xend=1, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=0.75, xend=0.75, y=0, yend=-Inf, col=col_v21, size=1)
  pcch <- ggplot(data.frame(Gamma12=constr_lgpdctc12.ch, Gamma21=constr_lgpdctc21.ch), aes(x=x) ) +
    geom_histogram(aes(x=Gamma12, y=..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v12) +
    geom_histogram(aes(x=Gamma21, y=-..density..), binwidth=0.005,boundary=1,closed="right", fill=col_v21) + labs(title="X1 -> X2 with confounder", x=NULL, y=NULL) +
    expand_limits(x=plot_xlimsc) + geom_segment(x=1, xend=1, y=0, yend=Inf, col=col_v12, size=1) + geom_segment(x=7/8, xend=7/8, y=0, yend=-Inf, col=col_v21, size=1)#+
  #geom_segment(x=0.75, xend=0.75, y=0, yend=-Inf, col=col_v21, size=1)
  pvals_constrlgpdctc <- ggarrange(pci,pcih,pcc,pcch, labels=c("A", "B", "C", "D"), ncol=2, nrow=2) %>%
    annotate_figure(left=text_grob("Density", rot = 90), bottom="Constrained fit LGPD CTCs")
  #plot(pvals_constrlgpdctc)
  Sys.sleep(0.5)
  ggsave(paste(sims_code, "_vals_constr_lgpdctc.png", sep=""), plot=pvals_constrlgpdctc, device="png", path=save_path, width=Wvals, height=Hvals, units="mm", dpi=300)
  Sys.sleep(0.5)
  
}


end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)

