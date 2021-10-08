
file_constr = "constr_LGPDCTC_Pa1-2_H1-4_1M2k1000.RData"
file_unconstr = "unconstr_CTCs_Pa1-2_H1-4_1M2k1000.RData"

export_to = "CTCs_Pa1-2_H1-4_1M2k1000.RData"


load(file_constr)
Sys.sleep(2)

load(file_unconstr)
Sys.sleep(2)

## RENAME VARIABLES
#ctc12.i <- ctc12.i
#ctc21.i <- ctc21.i
#ctc12.ih <- ctc12.ih
#ctc21.ih <- ctc21.ih
#ctc12.c <- ctc12.c
#ctc21.c <- ctc21.c
#ctc12.ch <- ctc12.ch
#ctc21.ch <- ctc21.ch
#gpdctc12.i <- gpdctc12.i
#gpdctc21.i <- gpdctc21.i
#gpdctc12.ih <- gpdctc12.ih
#gpdctc21.ih <- gpdctc21.ih
#gpdctc12.c <- gpdctc12.c
#gpdctc21.c <- gpdctc21.c
#gpdctc12.ch <- gpdctc12.ch
#gpdctc21.ch <- gpdctc21.ch
pfcorr_lgpdctc12.i <- lgpdctc12.i
pfcorr_lgpdctc21.i <- lgpdctc21.i
pfcorr_lgpdctc12.ih <- lgpdctc12.ih
pfcorr_lgpdctc21.ih <- lgpdctc21.ih
pfcorr_lgpdctc12.c <- lgpdctc12.c
pfcorr_lgpdctc21.c <- lgpdctc21.c
pfcorr_lgpdctc12.ch <- lgpdctc12.ch
pfcorr_lgpdctc21.ch <- lgpdctc21.ch
#constr_lgpdctc12.i <- constr_lgpdctc12.i
#constr_lgpdctc21.i <- constr_lgpdctc21.i
#constr_lgpdctc12.ih <- constr_lgpdctc12.ih
#constr_lgpdctc21.ih <- constr_lgpdctc21.ih
#constr_lgpdctc12.c <- constr_lgpdctc12.c
#constr_lgpdctc21.c <- constr_lgpdctc21.c
#constr_lgpdctc12.ch <- constr_lgpdctc12.ch
#constr_lgpdctc21.ch <- constr_lgpdctc21.ch

Sys.sleep(2)

save(ctc12.i,ctc21.i,ctc12.ih,ctc21.ih,ctc12.c,ctc21.c,ctc12.ch,ctc21.ch,
     gpdctc12.i,gpdctc21.i,gpdctc12.ih,gpdctc21.ih,gpdctc12.c,gpdctc21.c,gpdctc12.ch,gpdctc21.ch,
     pfcorr_lgpdctc12.i,pfcorr_lgpdctc21.i,pfcorr_lgpdctc12.ih,pfcorr_lgpdctc21.ih,
     pfcorr_lgpdctc12.c,pfcorr_lgpdctc21.c,pfcorr_lgpdctc12.ch,pfcorr_lgpdctc21.ch,
     constr_lgpdctc12.i,constr_lgpdctc21.i,constr_lgpdctc12.ih,constr_lgpdctc21.ih,
     constr_lgpdctc12.c,constr_lgpdctc21.c,constr_lgpdctc12.ch,constr_lgpdctc21.ch,
     file=export_to)


Sys.sleep(2)
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
