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
#library("plyr") library(reshape) library(ggplot2) library(tidyr) library(dplyr)
library("ggpubr")
#theme_set(theme_bw() + theme(legend.position = "top"))
library(qqplotr)

source("../R_functions/CTC_contribution_functions_OCP.R")

save_path <- "Results/Permutation_tests/power_plots_q90/"
sims_path <- "CTCs_simulations/Permutation_tests/q90/"
check_directory(save_path, recursive=TRUE)

distr_codes <- c("ln_H1p5_q90_10k2k_R1k", "ln_q90_10k2k_R1k", "Pa1-2_asym0p8_1_q90_10k2k_R1k", "Pa1-2_asym1_0p8_q90_10k2k_R1k",
                 "Pa1-2_H1-1_q90_10k2k_R1k", "Pa1-2_q90_10k2k_R1k", "Pa2p3_H1p8_q90_10k2k_R1k", "Pa3_H2_q90_10k2k_R1k",
                 "t4_H3_q90_10k2k_R1k", "t4_q90_10k2k_R1k")

start_time <- Sys.time()

for (distr_code in distr_codes) {
  
  source("../R_functions/CTC_contribution_functions_OCP.R")
  
  set.seed(1)
  
  beta_vals <- c(0,0.01,0.05,0.1,0.2)
  
  beta_codes <- paste0("_beta", sapply(beta_vals, function(x) {str_replace(toString(x),"([.])","p")}), "_1000")
  sims_codes <- paste0(distr_code, beta_codes)
  
  #Plots widths and heights for the 2 types of plots
  Wpval <- 200#200
  Hpval <- floor(1/3 * Wpval)
  Wqq <- 300#200
  Hqq <- floor(1/2 * Wqq)
  
  #Plot Colors
  col_dif <- rgb(83/255,81/255,84/255)# Greys: "grey30" rgb(128/255,133/255,133/255) rgb(83/255,81/255,84/255)
  col_v12 <- "#69b3a2"# Greens: "#69b3a2" rgb(132/255,186/255,91/255) rgb(62/255,150/255,81/255)
  col_v21 <- rgb(57/255,106/255,177/255)# Blues: "#404080" rgb(114/255,147/255,203/255) rgb(57/255,106/255,177/255)
  
  Pmc_ctc_t.c <- tibble(Pmc=double(), beta_21=double())
  Pmc_ctc_t.ch <- tibble(Pmc=double(), beta_21=double())
  Pmc_constr_ctc_t.c <- tibble(Pmc=double(), beta_21=double())
  Pmc_constr_ctc_t.ch <- tibble(Pmc=double(), beta_21=double())
  Pmc_pfcorr_ctc_t.c <- tibble(Pmc=double(), beta_21=double())
  Pmc_pfcorr_ctc_t.ch <- tibble(Pmc=double(), beta_21=double())
  
  for(i in 1:length(sims_codes)){
    
    sims_code <- sims_codes[i]
    beta21 <- beta_vals[i]
    
    load(paste0(sims_path, "CTCs_Permutation_test_power_", sims_code, ".RData"))
    
    cat("\n=================== Test powers: beta_21=",beta21," sim_code=",sims_code," ===================\n", sep="")
    
    cat("NAs (NP CTC):           A(",sum(is.na(Pmc_ctc.c)),") B(",sum(is.na(Pmc_ctc.ch)),")\n",
        "NAs (constr LGPDCTC):   A(",sum(is.na(Pmc_constr_lgpdctc.c)),") B(",sum(is.na(Pmc_constr_lgpdctc.ch)),")\n",
        "NAs (pf-corr. LGPDCTC): A(",sum(is.na(Pmc_pfcorr_lgpdctc.c)),") B(",sum(is.na(Pmc_pfcorr_lgpdctc.ch)),")\n", sep="")
    
    plot_xlims <- c(0,1)
    d.c <- ggplot(data.frame(Pmc=Pmc_ctc.c), aes(x=x) ) +
      geom_histogram(aes(x=Pmc, y=..density..), binwidth=0.05,boundary=1,closed="right", fill=col_dif) + labs(title="No confounder", x=NULL, y=NULL) +
      expand_limits(x=plot_xlims) + geom_vline(xintercept=0.05, col="black", size=1)
    d.ch <- ggplot(data.frame(Pmc=Pmc_ctc.ch), aes(x=x) ) +
      geom_histogram(aes(x=Pmc, y=..density..), binwidth=0.05,boundary=1,closed="right", fill=col_dif) + labs(title="With confounder", x=NULL, y=NULL) +
      expand_limits(x=plot_xlims) + geom_vline(xintercept=0.05, col="black", size=1)
    pPmc <- ggarrange(d.c, d.ch, labels=c("A", "B"), ncol=2, nrow=1) %>%
      annotate_figure(left=text_grob("Density", rot = 90), bottom="Non-parametric CTC test p-values")
    plot(pPmc)
    ggsave(paste(sims_code, "_Pmc_vals_npctc.pdf", sep=""), plot=pPmc, device="pdf", path=save_path, width=Wpval, height=Hpval, units="mm", dpi=300)
    
    
    d.c <- ggplot(data.frame(Pmc=Pmc_constr_lgpdctc.c), aes(x=x) ) +
      geom_histogram(aes(x=Pmc, y=..density..), binwidth=0.05,boundary=1,closed="right", fill=col_dif) + labs(title="No confounder", x=NULL, y=NULL) +
      expand_limits(x=plot_xlims) + geom_vline(xintercept=0.05, col="black", size=1)
    d.ch <- ggplot(data.frame(Pmc=Pmc_constr_lgpdctc.ch), aes(x=x) ) +
      geom_histogram(aes(x=Pmc, y=..density..), binwidth=0.05,boundary=1,closed="right", fill=col_dif) + labs(title="With confounder", x=NULL, y=NULL) +
      expand_limits(x=plot_xlims) + geom_vline(xintercept=0.05, col="black", size=1)
    pPmc <- ggarrange(d.c, d.ch, labels=c("A", "B"), ncol=2, nrow=1) %>%
      annotate_figure(left=text_grob("Density", rot = 90), bottom="Constrained fit LGPD CTC test p-values")
    plot(pPmc)
    ggsave(paste(sims_code, "_Pmc_vals_constr_lgpdctc.pdf", sep=""), plot=pPmc, device="pdf", path=save_path, width=Wpval, height=Hpval, units="mm", dpi=300)
    
    
    d.c <- ggplot(data.frame(Pmc=Pmc_pfcorr_lgpdctc.c), aes(x=x) ) +
      geom_histogram(aes(x=Pmc, y=..density..), binwidth=0.05,boundary=1,closed="right", fill=col_dif) + labs(title="No confounder", x=NULL, y=NULL) +
      expand_limits(x=plot_xlims) + geom_vline(xintercept=0.05, col="black", size=1)
    d.ch <- ggplot(data.frame(Pmc=Pmc_pfcorr_lgpdctc.ch), aes(x=x) ) +
      geom_histogram(aes(x=Pmc, y=..density..), binwidth=0.05,boundary=1,closed="right", fill=col_dif) + labs(title="With confounder", x=NULL, y=NULL) +
      expand_limits(x=plot_xlims) + geom_vline(xintercept=0.05, col="black", size=1)
    pPmc <- ggarrange(d.c, d.ch, labels=c("A", "B"), ncol=2, nrow=1) %>%
      annotate_figure(left=text_grob("Density", rot = 90), bottom="Post-fit corr. LGPD CTC test p-values")
    plot(pPmc)
    ggsave(paste(sims_code, "_Pmc_vals_pfcorr_lgpdctc.pdf", sep=""), plot=pPmc, device="pdf", path=save_path, width=Wpval, height=Hpval, units="mm", dpi=300)
    
    
    Pmc_ctc_t.c <- Pmc_ctc_t.c %>% add_row(Pmc=Pmc_ctc.c, beta_21=beta21)
    Pmc_ctc_t.ch <- Pmc_ctc_t.ch %>% add_row(Pmc=Pmc_ctc.ch, beta_21=beta21)
    Pmc_constr_ctc_t.c <- Pmc_constr_ctc_t.c %>% add_row(Pmc=Pmc_constr_lgpdctc.c, beta_21=beta21)
    Pmc_constr_ctc_t.ch <- Pmc_constr_ctc_t.ch %>% add_row(Pmc=Pmc_constr_lgpdctc.ch, beta_21=beta21)
    Pmc_pfcorr_ctc_t.c <- Pmc_pfcorr_ctc_t.c %>% add_row(Pmc=Pmc_pfcorr_lgpdctc.c, beta_21=beta21)
    Pmc_pfcorr_ctc_t.ch <- Pmc_pfcorr_ctc_t.ch %>% add_row(Pmc=Pmc_pfcorr_lgpdctc.ch, beta_21=beta21)
  }
  
  
  plot_qq_ctc_test <- function(ctc_tibble, confidence_method="ks", title_str=NULL, x_str=NULL, y_str=NULL){
    plot_xlims <- c(0,1)
    dpars <- list(min=0, max=1)
    detr <- FALSE
    bt <- confidence_method#"pointwise" "ks"
    qq_p <- ggplot(mapping = aes(sample=Pmc)) +#, color=beta_21
      stat_qq_band(data=filter(ctc_tibble,beta_21==0), distribution="unif", detrend=detr, identity=TRUE, bandType=bt, dparams=dpars, alpha=0.8) +
      stat_qq_line(data=filter(ctc_tibble,beta_21==0), distribution="unif", detrend=detr, identity=TRUE, dparams=dpars) +
      geom_hline(yintercept=0.05, col="black", size=1, linetype="solid") +
      stat_qq_point(data=filter(ctc_tibble,beta_21==0.2), distribution="unif", detrend=detr, identity=TRUE, dparams=dpars, mapping=aes(color="0.2"),size=1) +
      stat_qq_point(data=filter(ctc_tibble,beta_21==0.1), distribution="unif", detrend=detr, identity=TRUE, dparams=dpars, mapping=aes(color="0.1"),size=1) +
      stat_qq_point(data=filter(ctc_tibble,beta_21==0.05), distribution="unif", detrend=detr, identity=TRUE, dparams=dpars, mapping=aes(color="0.05"),size=1) +
      stat_qq_point(data=filter(ctc_tibble,beta_21==0.01), distribution="unif", detrend=detr, identity=TRUE, dparams=dpars, mapping=aes(color="0.01"),size=1) +
      stat_qq_point(data=filter(ctc_tibble,beta_21==0), distribution="unif", detrend=detr, identity=TRUE, dparams=dpars, mapping=aes(color="0"),size=1) +
      expand_limits(x=plot_xlims) + labs(title=title_str, x=x_str, y=y_str, color="Causal\nStrength") + scale_x_continuous(breaks=seq(0,10,0.2)) + 
      scale_y_continuous(breaks=seq(0,10,0.2)) + coord_fixed()#+ coord_flip()
    # + theme(legend.position=c(0.12, 0.74), legend.background=element_rect(fill="gray92",color="white"), legend.key=element_rect(fill="gray92", color=NA))
    return(qq_p)
  }
  
  p_qq_np.c <- plot_qq_ctc_test(Pmc_ctc_t.c, "ks", y_str="No confounder")
  p_qq_np.ch <- plot_qq_ctc_test(Pmc_ctc_t.ch, "ks", x_str="Non-parametric", y_str="With confounder")
  p_qq_constr.c <- plot_qq_ctc_test(Pmc_constr_ctc_t.c, "ks")
  p_qq_constr.ch <- plot_qq_ctc_test(Pmc_constr_ctc_t.ch, "ks", x_str="Constrained fit")
  p_qq_pfcorr.c <- plot_qq_ctc_test(Pmc_pfcorr_ctc_t.c, "ks")
  p_qq_pfcorr.ch <- plot_qq_ctc_test(Pmc_pfcorr_ctc_t.ch, "ks", x_str="Post-fit corrected")
  pqq <- ggarrange(p_qq_np.c, p_qq_constr.c, p_qq_pfcorr.c, p_qq_np.ch, p_qq_constr.ch, p_qq_pfcorr.ch, labels=NULL, ncol=3, nrow=2,
                   common.legend=TRUE, legend="right", align="hv") %>%
    annotate_figure(left=text_grob("CTC p-value Sample Quantiles", rot = 90), bottom="Standard Uniform Quantiles")
  plot(pqq)
  ggsave(paste(substr(sims_codes[1],1,str_length(sims_codes[1])-6), "s_Pmc_QQ.pdf", sep=""), plot=pqq, device="pdf",
         path=save_path, width=250, height=155, units="mm", dpi=300)
  ggsave(paste(substr(sims_codes[1],1,str_length(sims_codes[1])-6), "s_Pmc_QQ.png", sep=""), plot=pqq, device="png",
         path=save_path, width=250, height=155, units="mm", dpi=300)
  
  rm(list = setdiff(ls(), union(lsf.str(), c("distr_code", "distr_codes", "save_path", "sims_path", "start_time"))))
}


end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)

