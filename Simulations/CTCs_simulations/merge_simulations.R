

files = c("CTCs_Pa3_H1p5_q90_1M2k1000_1.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_2.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_3.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_4.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_5.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_6.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_7.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_8.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_9.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_10.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_11.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_12.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_13.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_14.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_15.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_16.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_17.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_18.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_19.RData",
          "CTCs_Pa3_H1p5_q90_1M2k1000_20.RData")

export_to = "CTCs_Pa3_H1p5_q90_1M2k1000.RData"

m_ctc12.i <- NULL
m_ctc21.i <- NULL
m_ctc12.ih <- NULL
m_ctc21.ih <- NULL
m_ctc12.c <- NULL
m_ctc21.c <- NULL
m_ctc12.ch <- NULL
m_ctc21.ch <- NULL
m_gpdctc12.i <- NULL
m_gpdctc21.i <- NULL
m_gpdctc12.ih <- NULL
m_gpdctc21.ih <- NULL
m_gpdctc12.c <- NULL
m_gpdctc21.c <- NULL
m_gpdctc12.ch <- NULL
m_gpdctc21.ch <- NULL
m_pfcorr_lgpdctc12.i <- NULL
m_pfcorr_lgpdctc21.i <- NULL
m_pfcorr_lgpdctc12.ih <- NULL
m_pfcorr_lgpdctc21.ih <- NULL
m_pfcorr_lgpdctc12.c <- NULL
m_pfcorr_lgpdctc21.c <- NULL
m_pfcorr_lgpdctc12.ch <- NULL
m_pfcorr_lgpdctc21.ch <- NULL
m_constr_lgpdctc12.i <- NULL
m_constr_lgpdctc21.i <- NULL
m_constr_lgpdctc12.ih <- NULL
m_constr_lgpdctc21.ih <- NULL
m_constr_lgpdctc12.c <- NULL
m_constr_lgpdctc21.c <- NULL
m_constr_lgpdctc12.ch <- NULL
m_constr_lgpdctc21.ch <- NULL

for (file in files){
  load(file)
  Sys.sleep(2)
  m_ctc12.i <- c(m_ctc12.i, ctc12.i)
  m_ctc21.i <- c(m_ctc21.i, ctc21.i)
  m_ctc12.ih <- c(m_ctc12.ih, ctc12.ih)
  m_ctc21.ih <- c(m_ctc21.ih, ctc21.ih)
  m_ctc12.c <- c(m_ctc12.c, ctc12.c)
  m_ctc21.c <- c(m_ctc21.c, ctc21.c)
  m_ctc12.ch <- c(m_ctc12.ch, ctc12.ch)
  m_ctc21.ch <- c(m_ctc21.ch, ctc21.ch)
  m_gpdctc12.i <- c(m_gpdctc12.i, gpdctc12.i)
  m_gpdctc21.i <- c(m_gpdctc21.i, gpdctc21.i)
  m_gpdctc12.ih <- c(m_gpdctc12.ih, gpdctc12.ih)
  m_gpdctc21.ih <- c(m_gpdctc21.ih, gpdctc21.ih)
  m_gpdctc12.c <- c(m_gpdctc12.c, gpdctc12.c)
  m_gpdctc21.c <- c(m_gpdctc21.c, gpdctc21.c)
  m_gpdctc12.ch <- c(m_gpdctc12.ch, gpdctc12.ch)
  m_gpdctc21.ch <- c(m_gpdctc21.ch, gpdctc21.ch)
  m_pfcorr_lgpdctc12.i <- c(m_pfcorr_lgpdctc12.i, pfcorr_lgpdctc12.i)
  m_pfcorr_lgpdctc21.i <- c(m_pfcorr_lgpdctc21.i, pfcorr_lgpdctc21.i)
  m_pfcorr_lgpdctc12.ih <- c(m_pfcorr_lgpdctc12.ih, pfcorr_lgpdctc12.ih)
  m_pfcorr_lgpdctc21.ih <- c(m_pfcorr_lgpdctc21.ih, pfcorr_lgpdctc21.ih)
  m_pfcorr_lgpdctc12.c <- c(m_pfcorr_lgpdctc12.c, pfcorr_lgpdctc12.c)
  m_pfcorr_lgpdctc21.c <- c(m_pfcorr_lgpdctc21.c, pfcorr_lgpdctc21.c)
  m_pfcorr_lgpdctc12.ch <- c(m_pfcorr_lgpdctc12.ch, pfcorr_lgpdctc12.ch)
  m_pfcorr_lgpdctc21.ch <- c(m_pfcorr_lgpdctc21.ch, pfcorr_lgpdctc21.ch)
  m_constr_lgpdctc12.i <- c(m_constr_lgpdctc12.i, constr_lgpdctc12.i)
  m_constr_lgpdctc21.i <- c(m_constr_lgpdctc21.i, constr_lgpdctc21.i)
  m_constr_lgpdctc12.ih <- c(m_constr_lgpdctc12.ih, constr_lgpdctc12.ih)
  m_constr_lgpdctc21.ih <- c(m_constr_lgpdctc21.ih, constr_lgpdctc21.ih)
  m_constr_lgpdctc12.c <- c(m_constr_lgpdctc12.c, constr_lgpdctc12.c)
  m_constr_lgpdctc21.c <- c(m_constr_lgpdctc21.c, constr_lgpdctc21.c)
  m_constr_lgpdctc12.ch <- c(m_constr_lgpdctc12.ch, constr_lgpdctc12.ch)
  m_constr_lgpdctc21.ch <- c(m_constr_lgpdctc21.ch, constr_lgpdctc21.ch)
  Sys.sleep(2)
}

ctc12.i <- m_ctc12.i
ctc21.i <- m_ctc21.i
ctc12.ih <- m_ctc12.ih
ctc21.ih <- m_ctc21.ih
ctc12.c <- m_ctc12.c
ctc21.c <- m_ctc21.c
ctc12.ch <- m_ctc12.ch
ctc21.ch <- m_ctc21.ch
gpdctc12.i <- m_gpdctc12.i
gpdctc21.i <- m_gpdctc21.i
gpdctc12.ih <- m_gpdctc12.ih
gpdctc21.ih <- m_gpdctc21.ih
gpdctc12.c <- m_gpdctc12.c
gpdctc21.c <- m_gpdctc21.c
gpdctc12.ch <- m_gpdctc12.ch
gpdctc21.ch <- m_gpdctc21.ch
pfcorr_lgpdctc12.i <- m_pfcorr_lgpdctc12.i
pfcorr_lgpdctc21.i <- m_pfcorr_lgpdctc21.i
pfcorr_lgpdctc12.ih <- m_pfcorr_lgpdctc12.ih
pfcorr_lgpdctc21.ih <- m_pfcorr_lgpdctc21.ih
pfcorr_lgpdctc12.c <- m_pfcorr_lgpdctc12.c
pfcorr_lgpdctc21.c <- m_pfcorr_lgpdctc21.c
pfcorr_lgpdctc12.ch <- m_pfcorr_lgpdctc12.ch
pfcorr_lgpdctc21.ch <- m_pfcorr_lgpdctc21.ch
constr_lgpdctc12.i <- m_constr_lgpdctc12.i
constr_lgpdctc21.i <- m_constr_lgpdctc21.i
constr_lgpdctc12.ih <- m_constr_lgpdctc12.ih
constr_lgpdctc21.ih <- m_constr_lgpdctc21.ih
constr_lgpdctc12.c <- m_constr_lgpdctc12.c
constr_lgpdctc21.c <- m_constr_lgpdctc21.c
constr_lgpdctc12.ch <- m_constr_lgpdctc12.ch
constr_lgpdctc21.ch <- m_constr_lgpdctc21.ch

if(!is.null(export_to)){
  save(ctc12.i,ctc21.i,ctc12.ih,ctc21.ih,ctc12.c,ctc21.c,ctc12.ch,ctc21.ch,
       gpdctc12.i,gpdctc21.i,gpdctc12.ih,gpdctc21.ih,gpdctc12.c,gpdctc21.c,gpdctc12.ch,gpdctc21.ch,
       pfcorr_lgpdctc12.i,pfcorr_lgpdctc21.i,pfcorr_lgpdctc12.ih,pfcorr_lgpdctc21.ih,
       pfcorr_lgpdctc12.c,pfcorr_lgpdctc21.c,pfcorr_lgpdctc12.ch,pfcorr_lgpdctc21.ch,
       constr_lgpdctc12.i,constr_lgpdctc21.i,constr_lgpdctc12.ih,constr_lgpdctc21.ih,
       constr_lgpdctc12.c,constr_lgpdctc21.c,constr_lgpdctc12.ch,constr_lgpdctc21.ch,
       file=export_to)
}
#if(!is.null(export_to)){
#  save(constr_lgpdctc12.i,constr_lgpdctc21.i,constr_lgpdctc12.ih,constr_lgpdctc21.ih,
#       constr_lgpdctc12.c,constr_lgpdctc21.c,constr_lgpdctc12.ch,constr_lgpdctc21.ch,
#       file=export_to)
#}

print(length(unique(ctc12.i)))
print(length(unique(ctc21.i)))
print(length(unique(ctc12.ih)))
print(length(unique(ctc21.ih)))
print(length(unique(ctc12.c)))
print(length(unique(ctc21.c)))
print(length(unique(ctc12.ch)))
print(length(unique(ctc21.ch)))
print(length(unique(gpdctc12.i)))
print(length(unique(gpdctc21.i)))
print(length(unique(gpdctc12.ih)))
print(length(unique(gpdctc21.ih)))
print(length(unique(gpdctc12.c)))
print(length(unique(gpdctc21.c)))
print(length(unique(gpdctc12.ch)))
print(length(unique(gpdctc21.ch)))
print(length(unique(pfcorr_lgpdctc12.i)))
print(length(unique(pfcorr_lgpdctc21.i)))
print(length(unique(pfcorr_lgpdctc12.ih)))
print(length(unique(pfcorr_lgpdctc21.ih)))
print(length(unique(pfcorr_lgpdctc12.c)))
print(length(unique(pfcorr_lgpdctc21.c)))
print(length(unique(pfcorr_lgpdctc12.ch)))
print(length(unique(pfcorr_lgpdctc21.ch)))
print(length(unique(constr_lgpdctc12.i)))
print(length(unique(constr_lgpdctc21.i)))
print(length(unique(constr_lgpdctc12.ih)))
print(length(unique(constr_lgpdctc21.ih)))
print(length(unique(constr_lgpdctc12.c)))
print(length(unique(constr_lgpdctc21.c)))
print(length(unique(constr_lgpdctc12.ch)))
print(length(unique(constr_lgpdctc21.ch)))
