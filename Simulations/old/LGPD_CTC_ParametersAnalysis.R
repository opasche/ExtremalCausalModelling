library(tidyverse)
library("evd")
library("evdbayes")
library("ismev")
library(latex2exp)
library(causalXtreme)
library(EnvStats)
library(boot)
library("ggpubr")
library(maxLik)

source("../R_functions/CTC_contribution_functions_OCP.R")

set.seed(1)

start_time <- Sys.time()

## Simulations
cat("\n=================== SIMULATIONS: Linear scale GPD CTC distributions ===================\n")
t_df <- 4
n_sim <- 5
size_sim <- 1000000
param.k <- 2*floor(size_sim^0.4) #floor(size_sim^0.4)

lgpdctc12.i <- list()
lgpdctc12.ih <- list()
lgpdctc12.c <- list()
lgpdctc21.c <- list()
lgpdctc12.ch <- list()
lgpdctc21.ch <- list()
t_lgpdctc12.i <- list()
t_lgpdctc12.ih <- list()
t_lgpdctc12.c <- list()
t_lgpdctc21.c <- list()
t_lgpdctc12.ch <- list()
t_lgpdctc21.ch <- list()
l_lgpdctc12.i <- list()
l_lgpdctc12.ih <- list()
l_lgpdctc12.c <- list()
l_lgpdctc21.c <- list()
l_lgpdctc12.ch <- list()
l_lgpdctc21.ch <- list()

H <- matrix(rt(n_sim*size_sim, df=t_df), c(n_sim,size_sim))
X1.i <- matrix(rt(n_sim*size_sim, df=t_df), c(n_sim,size_sim))
X2.i <- matrix(rt(n_sim*size_sim, df=t_df), c(n_sim,size_sim))
#rlnorm(n_sim*size_sim, meanlog=0, sdlog=1)#rt(n_sim*size_sim, df=t_df)#rpareto(n_sim*size_sim, location=1, shape=2)

cat("\n================ X1 indep. X2: ================\n")
for (i in 1:n_sim){
  lgpdctc12.i[[i]] <- LGPD_causal_tail_coeff_diagnostic(X1.i[i,], X2.i[i,], H[i,], k=param.k, parametric_F1=TRUE)
  t_lgpdctc12.i[[i]] <- student_constrained_LGPD_CTC(X1.i[i,], X2.i[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE, return_MLEs=TRUE)
  l_lgpdctc12.i[[i]] <- constrained_LGPD_CTC(X1.i[i,], X2.i[i,], H[i,], k=param.k, parametric_F1=TRUE, return_MLEs=TRUE)
}
cat("\n==== Gamma_12 ====\n")
cat("Post-fit correction:\n")
print(unnest_wider(tibble(lgpdctc12.i),"lgpdctc12.i"))
cat("Student box-constraints:\n")
print(unnest_wider(tibble(t_lgpdctc12.i),"t_lgpdctc12.i"))
cat("Linear constraints:\n")
print(unnest_wider(tibble(l_lgpdctc12.i),"l_lgpdctc12.i"))

cat("\n================ X1 X2 with confounder: ================\n")
X1.ih <- H + X1.i
X2.ih <- H + X2.i
for (i in 1:n_sim){
  lgpdctc12.ih[[i]] <- LGPD_causal_tail_coeff_diagnostic(X1.ih[i,], X2.ih[i,], H[i,], k=param.k, parametric_F1=TRUE)
  t_lgpdctc12.ih[[i]] <- student_constrained_LGPD_CTC(X1.ih[i,], X2.ih[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE, return_MLEs=TRUE)
  l_lgpdctc12.ih[[i]] <- constrained_LGPD_CTC(X1.ih[i,], X2.ih[i,], H[i,], k=param.k, parametric_F1=TRUE, return_MLEs=TRUE)
}
cat("\n==== Gamma_12 ====\n")
cat("Post-fit correction:\n")
print(unnest_wider(tibble(lgpdctc12.ih),"lgpdctc12.ih"))
cat("Student box-constraints:\n")
print(unnest_wider(tibble(t_lgpdctc12.ih),"t_lgpdctc12.ih"))
cat("Linear constraints:\n")
print(unnest_wider(tibble(l_lgpdctc12.ih),"l_lgpdctc12.ih"))

cat("\n================ X1 -> X2: ================\n")
X1.c <- X1.i# no copy (same tracemem())
X2.c <- X1.c + X2.i
for (i in 1:n_sim){
  lgpdctc12.c[[i]] <- LGPD_causal_tail_coeff_diagnostic(X1.c[i,], X2.c[i,], H[i,], k=param.k, parametric_F1=TRUE)
  lgpdctc21.c[[i]] <- LGPD_causal_tail_coeff_diagnostic(X2.c[i,], X1.c[i,], H[i,], k=param.k, parametric_F1=TRUE)
  t_lgpdctc12.c[[i]] <- student_constrained_LGPD_CTC(X1.c[i,], X2.c[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE, return_MLEs=TRUE)
  t_lgpdctc21.c[[i]] <- student_constrained_LGPD_CTC(X2.c[i,], X1.c[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE, return_MLEs=TRUE)
  l_lgpdctc12.c[[i]] <- constrained_LGPD_CTC(X1.c[i,], X2.c[i,], H[i,], k=param.k, parametric_F1=TRUE, return_MLEs=TRUE)
  l_lgpdctc21.c[[i]] <- constrained_LGPD_CTC(X2.c[i,], X1.c[i,], H[i,], k=param.k, parametric_F1=TRUE, return_MLEs=TRUE)
}
cat("\n==== Gamma_12 ====\n")
cat("Post-fit correction:\n")
print(unnest_wider(tibble(lgpdctc12.c),"lgpdctc12.c"))
cat("Student box-constraints:\n")
print(unnest_wider(tibble(t_lgpdctc12.c),"t_lgpdctc12.c"))
cat("Linear constraints:\n")
print(unnest_wider(tibble(l_lgpdctc12.c),"l_lgpdctc12.c"))
cat("\n==== Gamma_21 ====\n")
cat("Post-fit correction:\n")
print(unnest_wider(tibble(lgpdctc21.c),"lgpdctc21.c"))
cat("Student box-constraints:\n")
print(unnest_wider(tibble(t_lgpdctc21.c),"t_lgpdctc21.c"))
cat("Linear constraints:\n")
print(unnest_wider(tibble(l_lgpdctc21.c),"l_lgpdctc21.c"))

cat("\n================ X1 -> X2 with confounder: ================\n")
X1.ch <- X1.ih# no copy (same tracemem())
X2.ch <- X1.ch + H + X2.i
for (i in 1:n_sim){
  lgpdctc12.ch[[i]] <- LGPD_causal_tail_coeff_diagnostic(X1.ch[i,], X2.ch[i,], H[i,], k=param.k, parametric_F1=TRUE)
  lgpdctc21.ch[[i]] <- LGPD_causal_tail_coeff_diagnostic(X2.ch[i,], X1.ch[i,], H[i,], k=param.k, parametric_F1=TRUE)
  t_lgpdctc12.ch[[i]] <- student_constrained_LGPD_CTC(X1.ch[i,], X2.ch[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE, return_MLEs=TRUE)
  t_lgpdctc21.ch[[i]] <- student_constrained_LGPD_CTC(X2.ch[i,], X1.ch[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE, return_MLEs=TRUE)
  l_lgpdctc12.ch[[i]] <- constrained_LGPD_CTC(X1.ch[i,], X2.ch[i,], H[i,], k=param.k, parametric_F1=TRUE, return_MLEs=TRUE)
  l_lgpdctc21.ch[[i]] <- constrained_LGPD_CTC(X2.ch[i,], X1.ch[i,], H[i,], k=param.k, parametric_F1=TRUE, return_MLEs=TRUE)
}
cat("\n==== Gamma_12 ====\n")
cat("Post-fit correction:\n")
print(unnest_wider(tibble(lgpdctc12.ch),"lgpdctc12.ch"))
cat("Student box-constraints:\n")
print(unnest_wider(tibble(t_lgpdctc12.ch),"t_lgpdctc12.ch"))
cat("Linear constraints:\n")
print(unnest_wider(tibble(l_lgpdctc12.ch),"l_lgpdctc12.ch"))
cat("\n==== Gamma_21 ====\n")
cat("Post-fit correction:\n")
print(unnest_wider(tibble(lgpdctc21.ch),"lgpdctc21.ch"))
cat("Student box-constraints:\n")
print(unnest_wider(tibble(t_lgpdctc21.ch),"t_lgpdctc21.ch"))
cat("Linear constraints:\n")
print(unnest_wider(tibble(l_lgpdctc21.ch),"l_lgpdctc21.ch"))



end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)

#Sys.sleep(300)
#system('shutdown -t 30 -s')

