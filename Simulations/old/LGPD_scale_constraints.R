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


#upper=c(Inf,u1/(t_df*max(Hcs)),Inf)
#upper=c(Inf,u2/(t_df*max(Hcs)),Inf)
#lower=c(-Inf,-u1/(t_df*max(Hcs)),-Inf), upper=c(Inf,-u1/(t_df*min(Hcs)),Inf)


## Simulations
cat("\n=================== SIMULATIONS: Linear scale GPD unbounded and bounded fits ===================\n")
t_df <- 4
n_sim <- 200
size_sim <- 1000000
param.k <- 2*floor(size_sim^0.4) #floor(size_sim^0.4)

save_path <- "Results/constrained_lgpdctc/constrained_lgpdctc_1M2k200.RData"

H <- matrix(rt(n_sim*size_sim, df=t_df), c(n_sim,size_sim))
X1.i <- matrix(rt(n_sim*size_sim, df=t_df), c(n_sim,size_sim))
X2.i <- matrix(rt(n_sim*size_sim, df=t_df), c(n_sim,size_sim))
#rlnorm(n_sim*size_sim, meanlog=0, sdlog=1)#rt(n_sim*size_sim, df=t_df)#rpareto(n_sim*size_sim, location=1, shape=2)

cat("\n==== X1 indep. X2: ====\n")
minimums.i <- tibble(min_scale_F1 = double(), min_scale_F2 = double(), min_scale_F1_tconstr = double(), min_scale_F2_tconstr = double(),
                     min_scale_F1_lconstr = double(), min_scale_F2_lconstr = double())
for (i in 1:n_sim){
  Hcs <- scale(H[i,])
  u1 <- quantile(X1.i[i,],0.95)
  u2 <- quantile(X2.i[i,],0.95)
  fit1 <- gpd.fit(X1.i[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE)#,siglink=exp)
  min1 <- min(fit1$mle[1]+fit1$mle[2]*Hcs)
  fit2 <- gpd.fit(X2.i[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE)#,siglink=exp)
  min2 <- min(fit2$mle[1]+fit2$mle[2]*Hcs)
  cfit1 <- gpd_constraints_fit_optim(X1.i[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE, method = "L-BFGS-B",
                                     lower=c(-Inf,-u1/(t_df*max(Hcs)),-Inf), upper=c(Inf,-u1/(t_df*min(Hcs)),Inf))#,siglink=exp)
  cmin1 <- min(cfit1$mle[1]+cfit1$mle[2]*Hcs)
  cfit2 <- gpd_constraints_fit_optim(X2.i[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE, method = "L-BFGS-B",
                                     lower=c(-Inf,-u2/(t_df*max(Hcs)),-Inf), upper=c(Inf,-u2/(t_df*min(Hcs)),Inf))
  cmin2 <- min(cfit2$mle[1]+cfit2$mle[2]*Hcs)
  lfit1 <- gpd_constraints_fit_maxLik(X1.i[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE, return_se=FALSE,
                                      constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
  lmin1 <- min(cfit1$mle[1]+cfit1$mle[2]*Hcs)
  lfit2 <- gpd_constraints_fit_maxLik(X2.i[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE, return_se=FALSE,
                                      constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
  lmin2 <- min(cfit2$mle[1]+cfit2$mle[2]*Hcs)
  minimums.i <- minimums.i %>% add_row(min_scale_F1=min1, min_scale_F2=min2, min_scale_F1_tconstr=cmin1, min_scale_F2_tconstr=cmin2, min_scale_F1_lconstr=lmin1, min_scale_F2_lconstr=lmin2)
}
cat("First 15 values:\n")
print(head(minimums.i,15))
cat("Number of neative min. scales among the",n_sim,"simulations:\n")
print(apply(X = minimums.i<0, 2, sum))

cat("\n==== X1 X2 with confounder: ====\n")
X1.ih <- H + X1.i
X2.ih <- H + X2.i
minimums.ih <- tibble(min_scale_F1 = double(), min_scale_F2 = double(), min_scale_F1_tconstr = double(), min_scale_F2_tconstr = double(),
                      min_scale_F1_lconstr = double(), min_scale_F2_lconstr = double())
for (i in 1:n_sim){
  Hcs <- scale(H[i,])
  u1 <- quantile(X1.ih[i,],0.95)
  u2 <- quantile(X2.ih[i,],0.95)
  fit1 <- gpd.fit(X1.ih[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE)#,siglink=exp)
  min1 <- min(fit1$mle[1]+fit1$mle[2]*Hcs)
  fit2 <- gpd.fit(X2.ih[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE)#,siglink=exp)
  min2 <- min(fit2$mle[1]+fit2$mle[2]*Hcs)
  cfit1 <- gpd_constraints_fit_optim(X1.ih[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE, method = "L-BFGS-B",
                               lower=c(-Inf,-u1/(t_df*max(Hcs)),-Inf), upper=c(Inf,-u1/(t_df*min(Hcs)),Inf))
  cmin1 <- min(cfit1$mle[1]+cfit1$mle[2]*Hcs)
  cfit2 <- gpd_constraints_fit_optim(X2.ih[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE, method = "L-BFGS-B",
                               lower=c(-Inf,-u2/(t_df*max(Hcs)),-Inf), upper=c(Inf,-u2/(t_df*min(Hcs)),Inf))
  cmin2 <- min(cfit2$mle[1]+cfit2$mle[2]*Hcs)
  lfit1 <- gpd_constraints_fit_maxLik(X1.ih[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE, return_se=FALSE,
                                      constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
  lmin1 <- min(cfit1$mle[1]+cfit1$mle[2]*Hcs)
  lfit2 <- gpd_constraints_fit_maxLik(X2.ih[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE, return_se=FALSE,
                                      constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
  lmin2 <- min(cfit2$mle[1]+cfit2$mle[2]*Hcs)
  minimums.ih <- minimums.ih %>% add_row(min_scale_F1=min1, min_scale_F2=min2, min_scale_F1_tconstr=cmin1, min_scale_F2_tconstr=cmin2, min_scale_F1_lconstr=lmin1, min_scale_F2_lconstr=lmin2)
}
cat("First 15 values:\n")
print(head(minimums.ih,15))
cat("Number of neative min. scales among the",n_sim,"simulations:\n")
print(apply(X = minimums.ih<0, 2, sum))

cat("\n==== X1 -> X2: ====\n")
X1.c <- X1.i# no copy (same tracemem())
X2.c <- X1.c + X2.i
minimums.c <- tibble(min_scale_F1 = double(), min_scale_F2 = double(), min_scale_F1_tconstr = double(), min_scale_F2_tconstr = double(),
                     min_scale_F1_lconstr = double(), min_scale_F2_lconstr = double())
for (i in 1:n_sim){
  Hcs <- scale(H[i,])
  u1 <- quantile(X1.c[i,],0.95)
  u2 <- quantile(X2.c[i,],0.95)
  fit1 <- gpd.fit(X1.c[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE)#,siglink=exp)
  min1 <- min(fit1$mle[1]+fit1$mle[2]*Hcs)
  fit2 <- gpd.fit(X2.c[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE)#,siglink=exp)
  min2 <- min(fit2$mle[1]+fit2$mle[2]*Hcs)
  cfit1 <- gpd_constraints_fit_optim(X1.c[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE, method = "L-BFGS-B",
                               lower=c(-Inf,-u1/(t_df*max(Hcs)),-Inf), upper=c(Inf,-u1/(t_df*min(Hcs)),Inf))
  cmin1 <- min(cfit1$mle[1]+cfit1$mle[2]*Hcs)
  cfit2 <- gpd_constraints_fit_optim(X2.c[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE, method = "L-BFGS-B",
                               lower=c(-Inf,-u2/(t_df*max(Hcs)),-Inf), upper=c(Inf,-u2/(t_df*min(Hcs)),Inf))
  cmin2 <- min(cfit2$mle[1]+cfit2$mle[2]*Hcs)
  lfit1 <- gpd_constraints_fit_maxLik(X1.c[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE, return_se=FALSE,
                                      constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
  lmin1 <- min(cfit1$mle[1]+cfit1$mle[2]*Hcs)
  lfit2 <- gpd_constraints_fit_maxLik(X2.c[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE, return_se=FALSE,
                                      constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
  lmin2 <- min(cfit2$mle[1]+cfit2$mle[2]*Hcs)
  minimums.c <- minimums.c %>% add_row(min_scale_F1=min1, min_scale_F2=min2, min_scale_F1_tconstr=cmin1, min_scale_F2_tconstr=cmin2, min_scale_F1_lconstr=lmin1, min_scale_F2_lconstr=lmin2)
}
cat("First 15 values:\n")
print(head(minimums.c,15))
cat("Number of neative min. scales among the",n_sim,"simulations:\n")
print(apply(X = minimums.c<0, 2, sum))

cat("\n==== X1 -> X2 with confounder: ====\n")
X1.ch <- X1.ih# no copy (same tracemem())
X2.ch <- X1.ch + H + X2.i
minimums.ch <- tibble(min_scale_F1 = double(), min_scale_F2 = double(), min_scale_F1_tconstr = double(), min_scale_F2_tconstr = double(),
                      min_scale_F1_lconstr = double(), min_scale_F2_lconstr = double())
for (i in 1:n_sim){
  Hcs <- scale(H[i,])
  u1 <- quantile(X1.ch[i,],0.95)
  u2 <- quantile(X2.ch[i,],0.95)
  fit1 <- gpd.fit(X1.ch[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE)#,siglink=exp)
  min1 <- min(fit1$mle[1]+fit1$mle[2]*Hcs)
  fit2 <- gpd.fit(X2.ch[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE)#,siglink=exp)
  min2 <- min(fit2$mle[1]+fit2$mle[2]*Hcs)
  cfit1 <- gpd_constraints_fit_optim(X1.ch[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE, method = "L-BFGS-B",
                               lower=c(-Inf,-u1/(t_df*max(Hcs)),-Inf), upper=c(Inf,-u1/(t_df*min(Hcs)),Inf))
  cmin1 <- min(cfit1$mle[1]+cfit1$mle[2]*Hcs)
  cfit2 <- gpd_constraints_fit_optim(X2.ch[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE, method = "L-BFGS-B",
                               lower=c(-Inf,-u2/(t_df*max(Hcs)),-Inf), upper=c(Inf,-u2/(t_df*min(Hcs)),Inf))
  cmin2 <- min(cfit2$mle[1]+cfit2$mle[2]*Hcs)
  lfit1 <- gpd_constraints_fit_maxLik(X1.ch[i,],threshold=u1,ydat=as.matrix(Hcs),sigl=1, show=FALSE, return_se=FALSE,
                                      constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
  lmin1 <- min(cfit1$mle[1]+cfit1$mle[2]*Hcs)
  lfit2 <- gpd_constraints_fit_maxLik(X2.ch[i,],threshold=u2,ydat=as.matrix(Hcs),sigl=1, show=FALSE, return_se=FALSE,
                                      constraints=list(ineqA=matrix(c(1,min(Hcs),0,1,max(Hcs),0), nrow=2, ncol=3, byrow=TRUE), ineqB=matrix(0, nrow=2, ncol=1, byrow=TRUE)))
  lmin2 <- min(cfit2$mle[1]+cfit2$mle[2]*Hcs)
  minimums.ch <- minimums.ch %>% add_row(min_scale_F1=min1, min_scale_F2=min2, min_scale_F1_tconstr=cmin1, min_scale_F2_tconstr=cmin2, min_scale_F1_lconstr=lmin1, min_scale_F2_lconstr=lmin2)
}
cat("First 15 values:\n")
print(head(minimums.ch,15))
cat("Number of neative min. scales among the",n_sim,"simulations:\n")
print(apply(X = minimums.ch<0, 2, sum))



## fit bounded GPD CTC

t_lgpdctc12.i <- rep(as.double(NA), n_sim)
t_lgpdctc21.i <- rep(as.double(NA), n_sim)
t_lgpdctc12.ih <- rep(as.double(NA), n_sim)
t_lgpdctc21.ih <- rep(as.double(NA), n_sim)
t_lgpdctc12.c <- rep(as.double(NA), n_sim)
t_lgpdctc21.c <- rep(as.double(NA), n_sim)
t_lgpdctc12.ch <- rep(as.double(NA), n_sim)
t_lgpdctc21.ch <- rep(as.double(NA), n_sim)
l_lgpdctc12.i <- rep(as.double(NA), n_sim)
l_lgpdctc21.i <- rep(as.double(NA), n_sim)
l_lgpdctc12.ih <- rep(as.double(NA), n_sim)
l_lgpdctc21.ih <- rep(as.double(NA), n_sim)
l_lgpdctc12.c <- rep(as.double(NA), n_sim)
l_lgpdctc21.c <- rep(as.double(NA), n_sim)
l_lgpdctc12.ch <- rep(as.double(NA), n_sim)
l_lgpdctc21.ch <- rep(as.double(NA), n_sim)
cat("\nConstrained CTCs Calculation Progress out of",n_sim,":\n")
for (i in 1:n_sim){
  if(i %% 25 == 0){cat("",i)}
  t_lgpdctc12.i[i] <- student_constrained_LGPD_CTC(X1.i[i,], X2.i[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE)
  t_lgpdctc21.i[i] <- student_constrained_LGPD_CTC(X2.i[i,], X1.i[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE)
  t_lgpdctc12.ih[i] <- student_constrained_LGPD_CTC(X1.ih[i,], X2.ih[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE)
  t_lgpdctc21.ih[i] <- student_constrained_LGPD_CTC(X2.ih[i,], X1.ih[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE)
  t_lgpdctc12.c[i] <- student_constrained_LGPD_CTC(X1.c[i,], X2.c[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE)
  t_lgpdctc21.c[i] <- student_constrained_LGPD_CTC(X2.c[i,], X1.c[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE)
  t_lgpdctc12.ch[i] <- student_constrained_LGPD_CTC(X1.ch[i,], X2.ch[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE)
  t_lgpdctc21.ch[i] <- student_constrained_LGPD_CTC(X2.ch[i,], X1.ch[i,], H[i,], k=param.k, t_df=t_df, parametric_F1=TRUE)
  l_lgpdctc12.i[i] <- constrained_LGPD_CTC(X1.i[i,], X2.i[i,], H[i,], k=param.k, parametric_F1=TRUE)
  l_lgpdctc21.i[i] <- constrained_LGPD_CTC(X2.i[i,], X1.i[i,], H[i,], k=param.k, parametric_F1=TRUE)
  l_lgpdctc12.ih[i] <- constrained_LGPD_CTC(X1.ih[i,], X2.ih[i,], H[i,], k=param.k, parametric_F1=TRUE)
  l_lgpdctc21.ih[i] <- constrained_LGPD_CTC(X2.ih[i,], X1.ih[i,], H[i,], k=param.k, parametric_F1=TRUE)
  l_lgpdctc12.c[i] <- constrained_LGPD_CTC(X1.c[i,], X2.c[i,], H[i,], k=param.k, parametric_F1=TRUE)
  l_lgpdctc21.c[i] <- constrained_LGPD_CTC(X2.c[i,], X1.c[i,], H[i,], k=param.k, parametric_F1=TRUE)
  l_lgpdctc12.ch[i] <- constrained_LGPD_CTC(X1.ch[i,], X2.ch[i,], H[i,], k=param.k, parametric_F1=TRUE)
  l_lgpdctc21.ch[i] <- constrained_LGPD_CTC(X2.ch[i,], X1.ch[i,], H[i,], k=param.k, parametric_F1=TRUE)
}
cat("\n")
t_lgpddifftail.i <- t_lgpdctc12.i - t_lgpdctc21.i
t_lgpddifftail.ih <- t_lgpdctc12.ih - t_lgpdctc21.ih
t_lgpddifftail.c <- t_lgpdctc12.c - t_lgpdctc21.c
t_lgpddifftail.ch <- t_lgpdctc12.ch - t_lgpdctc21.ch
l_lgpddifftail.i <- l_lgpdctc12.i - l_lgpdctc21.i
l_lgpddifftail.ih <- l_lgpdctc12.ih - l_lgpdctc21.ih
l_lgpddifftail.c <- l_lgpdctc12.c - l_lgpdctc21.c
l_lgpddifftail.ch <- l_lgpdctc12.ch - l_lgpdctc21.ch

if(!is.null(save_path)){
  save(t_lgpdctc12.i,t_lgpdctc21.i,t_lgpdctc12.ih,t_lgpdctc21.ih,t_lgpdctc12.c,t_lgpdctc21.c,t_lgpdctc12.ch,t_lgpdctc21.ch,
       l_lgpdctc12.i,l_lgpdctc21.i,l_lgpdctc12.ih,l_lgpdctc21.ih,l_lgpdctc12.c,l_lgpdctc21.c,l_lgpdctc12.ch,l_lgpdctc21.ch,
       file=save_path)
}

par(mfrow=c(2,2))
hist(t_lgpddifftail.i, breaks = 20, freq = FALSE, main = "Independent variables", xlab = "Student constrained GPD CTC diff.")
abline(v=0, col="blue", lwd=1)#points(0, 0, lwd=4, col="blue")
hist(t_lgpddifftail.ih, breaks = 20, freq = FALSE, main = "Common hidden confounder", xlab = "Student constrained GPD CTC diff.")
abline(v=0, col="blue", lwd=1)#points(0, 0, lwd=4, col="blue")
hist(t_lgpddifftail.c, breaks = 20, freq = FALSE, main = "X1 -> X2", xlab = "Student constrained GPD CTC diff.")
abline(v=0.25, col="blue", lwd=1)#points(0, 0, lwd=4, col="blue")
hist(t_lgpddifftail.ch, breaks = 20, freq = FALSE, main = "X1 -> X2 with confounder", xlab = "Student constrained GPD CTC diff.")
abline(v=1/8, col="blue", lwd=1)#points(0, 0, lwd=4, col="blue")
par(mfrow=c(1,1))

pci <- ggplot(data.frame(Gamma12=t_lgpdctc12.i, Gamma21=t_lgpdctc21.i), aes(x=x) ) +
  geom_histogram( aes(x = Gamma12, y = ..density..), bins = 50, fill="#69b3a2" ) +
  geom_histogram( aes(x = Gamma21, y = -..density..), bins = 50, fill= "#404080") + labs(title="Independent variables", x="Student constrained GPD CTCs", y="Density")
pcih <- ggplot(data.frame(Gamma12=t_lgpdctc12.ih, Gamma21=t_lgpdctc21.ih), aes(x=x) ) +
  geom_histogram( aes(x = Gamma12, y = ..density..), bins = 50, fill="#69b3a2" ) +
  geom_histogram( aes(x = Gamma21, y = -..density..), bins = 50, fill= "#404080") + labs(title="Common hidden confounder", x="Student constrained GPD CTCs", y="Density")
pcc <- ggplot(data.frame(Gamma12=t_lgpdctc12.c, Gamma21=t_lgpdctc21.c), aes(x=x) ) +
  geom_histogram( aes(x = Gamma12, y = ..density..), bins = 50, fill="#69b3a2" ) +
  geom_histogram( aes(x = Gamma21, y = -..density..), bins = 50, fill= "#404080") + labs(title="X1 -> X2", x="Student constrained GPD CTCs", y="Density")
pcch <- ggplot(data.frame(Gamma12=t_lgpdctc12.ch, Gamma21=t_lgpdctc21.ch), aes(x=x) ) +
  geom_histogram( aes(x = Gamma12, y = ..density..), bins = 50, fill="#69b3a2" ) +
  geom_histogram( aes(x = Gamma21, y = -..density..), bins = 50, fill= "#404080") + labs(title="X1 -> X2 with confounder", x="Student constrained GPD CTCs", y="Density")
plot(ggarrange(pci,pcih,pcc,pcch, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2))

par(mfrow=c(2,2))
hist(l_lgpddifftail.i, breaks = 20, freq = FALSE, main = "Independent variables", xlab = "Linearly constrained GPD CTC diff.")
abline(v=0, col="blue", lwd=1)#points(0, 0, lwd=4, col="blue")
hist(l_lgpddifftail.ih, breaks = 20, freq = FALSE, main = "Common hidden confounder", xlab = "Linearly constrained GPD CTC diff.")
abline(v=0, col="blue", lwd=1)#points(0, 0, lwd=4, col="blue")
hist(l_lgpddifftail.c, breaks = 20, freq = FALSE, main = "X1 -> X2", xlab = "Linearly constrained GPD CTC diff.")
abline(v=0.25, col="blue", lwd=1)#points(0, 0, lwd=4, col="blue")
hist(l_lgpddifftail.ch, breaks = 20, freq = FALSE, main = "X1 -> X2 with confounder", xlab = "Linearly constrained GPD CTC diff.")
abline(v=1/8, col="blue", lwd=1)#points(0, 0, lwd=4, col="blue")
par(mfrow=c(1,1))

pci <- ggplot(data.frame(Gamma12=l_lgpdctc12.i, Gamma21=l_lgpdctc21.i), aes(x=x) ) +
  geom_histogram( aes(x = Gamma12, y = ..density..), bins = 50, fill="#69b3a2" ) +
  geom_histogram( aes(x = Gamma21, y = -..density..), bins = 50, fill= "#404080") + labs(title="Independent variables", x="Linearly constrained GPD CTCs", y="Density")
pcih <- ggplot(data.frame(Gamma12=l_lgpdctc12.ih, Gamma21=l_lgpdctc21.ih), aes(x=x) ) +
  geom_histogram( aes(x = Gamma12, y = ..density..), bins = 50, fill="#69b3a2" ) +
  geom_histogram( aes(x = Gamma21, y = -..density..), bins = 50, fill= "#404080") + labs(title="Common hidden confounder", x="Linearly constrained GPD CTCs", y="Density")
pcc <- ggplot(data.frame(Gamma12=l_lgpdctc12.c, Gamma21=l_lgpdctc21.c), aes(x=x) ) +
  geom_histogram( aes(x = Gamma12, y = ..density..), bins = 50, fill="#69b3a2" ) +
  geom_histogram( aes(x = Gamma21, y = -..density..), bins = 50, fill= "#404080") + labs(title="X1 -> X2", x="Linearly constrained GPD CTCs", y="Density")
pcch <- ggplot(data.frame(Gamma12=l_lgpdctc12.ch, Gamma21=l_lgpdctc21.ch), aes(x=x) ) +
  geom_histogram( aes(x = Gamma12, y = ..density..), bins = 50, fill="#69b3a2" ) +
  geom_histogram( aes(x = Gamma21, y = -..density..), bins = 50, fill= "#404080") + labs(title="X1 -> X2 with confounder", x="Linearly constrained GPD CTCs", y="Density")
plot(ggarrange(pci,pcih,pcc,pcch, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2))



end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)

#Sys.sleep(300)
#system('shutdown -t 30 -s')

