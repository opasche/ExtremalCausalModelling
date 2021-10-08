library(tidyverse)
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

save_path <- "Results/"


## Simulations
cat("\n=================== SIMULATIONS: CTC distributions ===================\n")
nb_sim_gen <- 1000
n_sim_boot <- 1000
size_sim <- 1000000
param.k <- 2*floor(size_sim^0.4)


H <- matrix(rt(nb_sim_gen*size_sim, df=4), c(nb_sim_gen,size_sim))
X1.i <- matrix(rt(nb_sim_gen*size_sim, df=4), c(nb_sim_gen,size_sim))
X2.i <- matrix(rt(nb_sim_gen*size_sim, df=4), c(nb_sim_gen,size_sim))
#rlnorm(nb_sim_gen*size_sim, meanlog=0, sdlog=1)#rt(nb_sim_gen*size_sim, df=4)#rpareto(nb_sim_gen*size_sim, location=1, shape=2)

load("CTCs_simulations/CTCs_t4_1M2k1000.RData")

cat("Independent variables:\n")
for (i in 1:nb_sim_gen){
  ctc12.i[i] <- causal_tail_coeff(X1.i[i,], X2.i[i,], k=param.k, both_tails = FALSE)
  ctc21.i[i] <- causal_tail_coeff(X2.i[i,], X1.i[i,], k=param.k, both_tails = FALSE)
}
cat("X1 X2 with confounder:\n")
X1.ih <- H + X1.i
X2.ih <- H + X2.i
for (i in 1:nb_sim_gen){
  ctc12.ih[i] <- causal_tail_coeff(X1.ih[i,], X2.ih[i,], k=param.k, both_tails = FALSE)
  ctc21.ih[i] <- causal_tail_coeff(X2.ih[i,], X1.ih[i,], k=param.k, both_tails = FALSE)
}
cat("X1 -> X2:\n")
X1.c <- X1.i# no copy (same tracemem())
X2.c <- X1.c + X2.i
for (i in 1:nb_sim_gen){
  ctc12.c[i] <- causal_tail_coeff(X1.c[i,], X2.c[i,], k=param.k, both_tails = FALSE)
  ctc21.c[i] <- causal_tail_coeff(X2.c[i,], X1.c[i,], k=param.k, both_tails = FALSE)
}
cat("X1 -> X2 with confounder:\n")
X1.ch <- X1.ih# no copy (same tracemem())
X2.ch <- X1.ch + H + X2.i
for (i in 1:nb_sim_gen){
  ctc12.ch[i] <- causal_tail_coeff(X1.ch[i,], X2.ch[i,], k=param.k, both_tails = FALSE)
  ctc21.ch[i] <- causal_tail_coeff(X2.ch[i,], X1.ch[i,], k=param.k, both_tails = FALSE)
}
difftail.i <- ctc12.i - ctc21.i
difftail.ih <- ctc12.ih - ctc21.ih
difftail.c <- ctc12.c - ctc21.c
difftail.ch <- ctc12.ch - ctc21.ch


## Simulations 2
cat("\n=================== SIMULATIONS: bootstrap distributions ===================\n")
# Bootstrap 95% CI for the statistics
# function to obtain the statistics from the data
both_coeffs <- function(data, indices, k=param.k, both_tails=TRUE) {#k=floor(length(data[,1])^0.4) nrow(data)
  d <- data[indices,] # allows boot to select sample
  coeff12 <- causal_tail_coeff(d[[1]],d[[2]], k=k, both_tails=both_tails)
  coeff21 <- causal_tail_coeff(d[[2]],d[[1]], k=k, both_tails=both_tails)
  return(c(coeff12 = coeff12, coeff21 = coeff21))
}
# function to obtain the difference statistic from the data
diff_coeffs <- function(data, indices, k=param.k, both_tails=TRUE) {#k=floor(length(data[,1])^0.4) nrow(data)
  d <- data[indices,] # allows boot to select sample
  coeff12 <- causal_tail_coeff(d[[1]],d[[2]], k=k, both_tails=both_tails)
  coeff21 <- causal_tail_coeff(d[[2]],d[[1]], k=k, both_tails=both_tails)
  return(coeff12 - coeff21)
}

nb_rebootstraps <- 10#nb_sim_gen #among the n_sim
#save space
rm(H)
X1.i <- X1.i[1:nb_rebootstraps,]
X1.ih <- X1.ih[1:nb_rebootstraps,]
X1.c <- X1.c[1:nb_rebootstraps,]
X1.ch <- X1.ch[1:nb_rebootstraps,]
X2.i <- X2.i[1:nb_rebootstraps,]
X2.ih <- X2.ih[1:nb_rebootstraps,]
X2.c <- X2.c[1:nb_rebootstraps,]
X2.ch <- X2.ch[1:nb_rebootstraps,]
#declare objects
bs_i <- matrix(nrow = n_sim_boot, ncol = nb_rebootstraps+1)
bs_ih <- matrix(nrow = n_sim_boot, ncol = nb_rebootstraps+1)
bs_c <- matrix(nrow = n_sim_boot, ncol = nb_rebootstraps+1)
bs_ch <- matrix(nrow = n_sim_boot, ncol = nb_rebootstraps+1)
#First columns are the simulations distributions - true value
bs_i[,1] <- difftail.i - 0
bs_ih[,1] <- difftail.ih - 0
bs_c[,1] <- difftail.c - (1-0.75)
bs_ch[,1] <- difftail.ch - (1-0.875)
cat("Itteration (out of ",nb_rebootstraps,"): ")
for (i in 1:nb_rebootstraps){
  cat(i,", ")
  bs_i[,i+1] <- boot(data=data.frame(X1=X1.i[i,], X2=X2.i[i,]), statistic=diff_coeffs, R=n_sim_boot, k=param.k, both_tails=FALSE)$t - difftail.i[i]
  bs_ih[,i+1] <- boot(data=data.frame(X1=X1.ih[i,], X2=X2.ih[i,]), statistic=diff_coeffs, R=n_sim_boot, k=param.k, both_tails=FALSE)$t - difftail.ih[i]
  bs_c[,i+1] <- boot(data=data.frame(X1=X1.c[i,], X2=X2.c[i,]), statistic=diff_coeffs, R=n_sim_boot, k=param.k, both_tails=FALSE)$t - difftail.c[i]
  bs_ch[,i+1] <- boot(data=data.frame(X1=X1.ch[i,], X2=X2.ch[i,]), statistic=diff_coeffs, R=n_sim_boot, k=param.k, both_tails=FALSE)$t - difftail.ch[i]
}
df <- data.frame(bs_i)
names(df) <- c("Simuations",1:nb_rebootstraps)
pbi <- df %>% gather(key="Origin", value="Value", factor_key=TRUE) %>% ggplot( aes(x=Origin, y=Value, fill=Origin)) +
  geom_violin(width=0.8, color="grey", alpha=0.5) +
  stat_boxplot(geom = "errorbar", width=0.4)+geom_boxplot(width=0.4, color="black", alpha=1) +
  labs(title = "Independent variables") +
  theme(legend.position="none", axis.title.x=element_blank(), panel.grid.major.x = element_blank())
df <- data.frame(bs_ih)
names(df) <- c("Simuations",1:nb_rebootstraps)
pbih <- df %>% gather(key="Origin", value="Value", factor_key=TRUE) %>% ggplot( aes(x=Origin, y=Value, fill=Origin)) +
  geom_violin(width=0.8, color="grey", alpha=0.5) +
  stat_boxplot(geom = "errorbar", width=0.4)+geom_boxplot(width=0.4, color="black", alpha=1) +
  labs(title = "Common hidden confounder") +
  theme(legend.position="none", axis.title.x=element_blank(), panel.grid.major.x = element_blank())
df <- data.frame(bs_c)
names(df) <- c("Simuations",1:nb_rebootstraps)
pbc <- df %>% gather(key="Origin", value="Value", factor_key=TRUE) %>% ggplot( aes(x=Origin, y=Value, fill=Origin)) +
  geom_violin(width=0.8, color="grey", alpha=0.5) +
  stat_boxplot(geom = "errorbar", width=0.4)+geom_boxplot(width=0.4, color="black", alpha=1) +
  labs(title = "X1 -> X2") +
  theme(legend.position="none", axis.title.x=element_blank(), panel.grid.major.x = element_blank())
df <- data.frame(bs_ch)
names(df) <- c("Simuations",1:nb_rebootstraps)
pbch <- df %>% gather(key="Origin", value="Value", factor_key=TRUE) %>% ggplot( aes(x=Origin, y=Value, fill=Origin)) +
  geom_violin(width=0.8, color="grey", alpha=0.5) +
  stat_boxplot(geom = "errorbar", width=0.4)+geom_boxplot(width=0.4, color="black", alpha=1) +
  labs(title = "X1 -> X2 with confounder") +
  theme(legend.position="none", axis.title.x=element_blank(), panel.grid.major.x = element_blank())
fboot <- ggarrange(pbi,pbih,pbc,pbch, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2) %>%
  annotate_figure(top = text_grob("Causal tail coefficient differences: simulations vs. bootstrap", face="bold", size=14))
plot(fboot)
ggsave(paste("ctc1", "_t4_10k2k_boot_last.png", sep=""), plot=fboot, device="png", path=save_path, width=220, height=floor(0.6*220), units="mm", dpi=300)

save(bs_i,bs_ih,bs_c,bs_ch, file="Results/ctc1_t4_10k2k_boot_last.RData")


end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)

#Sys.sleep(300)
#system('shutdown -t 30 -s')

