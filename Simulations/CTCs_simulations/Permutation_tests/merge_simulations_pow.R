

files = c("CTCs_Permutation_test_power_t4_10k2k_R1k_beta0p1_500.RData",
          "CTCs_Permutation_test_power_t4_10k2k_R1k_beta0p1_500_2.RData")#0 0p01 0p1

export_to = "CTCs_Permutation_test_power_t4_10k2k_R1k_beta0p1_1000.RData"

m_res_ctc.c <- list()
m_Pmc_ctc.c <- NULL
m_res_ctc.ch <- list()
m_Pmc_ctc.ch <- NULL
m_res_pfcorr_lgpdctc.c <- list()
m_Pmc_pfcorr_lgpdctc.c <- NULL
m_res_pfcorr_lgpdctc.ch <- list()
m_Pmc_pfcorr_lgpdctc.ch <- NULL
m_res_constr_lgpdctc.c <- list()
m_Pmc_constr_lgpdctc.c <- NULL
m_res_constr_lgpdctc.ch <- list()
m_Pmc_constr_lgpdctc.ch <- NULL

for (file in files){
  load(file)
  Sys.sleep(2)
  m_res_ctc.c <- c(m_res_ctc.c, res_ctc.c)
  m_Pmc_ctc.c <- c(m_Pmc_ctc.c, Pmc_ctc.c)
  m_res_ctc.ch <- c(m_res_ctc.ch, res_ctc.ch)
  m_Pmc_ctc.ch <- c(m_Pmc_ctc.ch, Pmc_ctc.ch)
  m_res_pfcorr_lgpdctc.c <- c(m_res_pfcorr_lgpdctc.c, res_pfcorr_lgpdctc.c)
  m_Pmc_pfcorr_lgpdctc.c <- c(m_Pmc_pfcorr_lgpdctc.c, Pmc_pfcorr_lgpdctc.c)
  m_res_pfcorr_lgpdctc.ch <- c(m_res_pfcorr_lgpdctc.ch, res_pfcorr_lgpdctc.ch)
  m_Pmc_pfcorr_lgpdctc.ch <- c(m_Pmc_pfcorr_lgpdctc.ch, Pmc_pfcorr_lgpdctc.ch)
  m_res_constr_lgpdctc.c <- c(m_res_constr_lgpdctc.c, res_constr_lgpdctc.c)
  m_Pmc_constr_lgpdctc.c <- c(m_Pmc_constr_lgpdctc.c, Pmc_constr_lgpdctc.c)
  m_res_constr_lgpdctc.ch <- c(m_res_constr_lgpdctc.ch, res_constr_lgpdctc.ch)
  m_Pmc_constr_lgpdctc.ch <- c(m_Pmc_constr_lgpdctc.ch, Pmc_constr_lgpdctc.ch)
  Sys.sleep(2)
}


res_ctc.c <- m_res_ctc.c
Pmc_ctc.c <- m_Pmc_ctc.c
res_ctc.ch <- m_res_ctc.ch
Pmc_ctc.ch <- m_Pmc_ctc.ch
res_pfcorr_lgpdctc.c <- m_res_pfcorr_lgpdctc.c
Pmc_pfcorr_lgpdctc.c <- m_Pmc_pfcorr_lgpdctc.c
res_pfcorr_lgpdctc.ch <- m_res_pfcorr_lgpdctc.ch
Pmc_pfcorr_lgpdctc.ch <- m_Pmc_pfcorr_lgpdctc.ch
res_constr_lgpdctc.c <- m_res_constr_lgpdctc.c
Pmc_constr_lgpdctc.c <- m_Pmc_constr_lgpdctc.c
res_constr_lgpdctc.ch <- m_res_constr_lgpdctc.ch
Pmc_constr_lgpdctc.ch <- m_Pmc_constr_lgpdctc.ch

if(!is.null(export_to)){
  save(res_ctc.c,Pmc_ctc.c,res_ctc.ch,Pmc_ctc.ch,
       res_pfcorr_lgpdctc.c,Pmc_pfcorr_lgpdctc.c,res_pfcorr_lgpdctc.ch,Pmc_pfcorr_lgpdctc.ch,
       res_constr_lgpdctc.c,Pmc_constr_lgpdctc.c,res_constr_lgpdctc.ch,Pmc_constr_lgpdctc.ch,
       file=export_to)
}

print(length(unique(res_ctc.c)))
print(length(unique(Pmc_ctc.c)))
print(length(unique(res_ctc.ch)))
print(length(unique(Pmc_ctc.ch)))
print(length(unique(res_pfcorr_lgpdctc.c)))
print(length(unique(Pmc_pfcorr_lgpdctc.c)))
print(length(unique(res_pfcorr_lgpdctc.ch)))
print(length(unique(Pmc_pfcorr_lgpdctc.ch)))
print(length(unique(res_constr_lgpdctc.c)))
print(length(unique(Pmc_constr_lgpdctc.c)))
print(length(unique(res_constr_lgpdctc.ch)))
print(length(unique(Pmc_constr_lgpdctc.ch)))
