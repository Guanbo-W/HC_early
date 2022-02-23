# survival outcome
source("~/ATE/Prep_Survival_ATE.R")
source("~/ATE/Methods_Survival_ATE.R")
table=NULL
betahat_T=NULL

# Frequentist
########## AFT
##### Separate Full Matching
start.time <- Sys.time()
tableALL.SFM.surv.aft=SFM.aft(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv)
rownames(tableALL.SFM.surv.aft)="SFM.aft"
tableALL.SFM.w.surv.aft=SFM.aft(X.w, nT.w, DA.w, A.w, N.w, Index.w_T_HC0, Index.w_T_HC1, ATT.surv.w)
rownames(tableALL.SFM.w.surv.aft)="SFM.w.aft"
end.time <- Sys.time()
end.time - start.time
######## SIPTW
start.time <- Sys.time()
tableALL.SIPTW.surv.aft=SIPTW.aft(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv)
rownames(tableALL.SIPTW.surv.aft)="SIPTW.aft"
tableALL.SIPTW.surv.aft
tableALL.SIPTW.w.surv.aft=SIPTW.aft(X.w, nT.w, DA.w, A.w, N.w, Index.w_T_HC0, Index.w_T_HC1, ATT.surv.w)
tableALL.SIPTW.w.surv.aft
rownames(tableALL.SIPTW.w.surv.aft)="SIPTW.w.aft"
end.time <- Sys.time()
end.time - start.time
##### JFM2
tableALL.JFM2.surv.aft=JFM2.aft(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv, D, DD, Index_noC)
rownames(tableALL.JFM2.surv.aft)="JFM2.aft"
tableALL.JFM2.w.surv.aft=JFM2.aft(X.w, nT.w, DA.w, A.w, N.w, Index.w_T_HC0, Index.w_T_HC1, ATT.surv.w, D.w, DD.w, Index.w_noC)
rownames(tableALL.JFM2.w.surv.aft)="JFM2.w.aft"
####### JIPTW2
tableALL.JIPTW2.surv.aft=JIPTW2.aft(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv, D, DD, Index_noC)
rownames(tableALL.JIPTW2.surv.aft)="JIPTW2.aft"
tableALL.JIPTW2.w.surv.aft=JIPTW2.aft(X.w, nT.w, DA.w, A.w, N.w, Index.w_T_HC0, Index.w_T_HC1, ATT.surv.w, D.w, DD.w, Index.w_noC)
rownames(tableALL.JIPTW2.w.surv.aft)="JIPTW2.w.aft"


########## CoxPH
##### Separate Full Matching
tableALL.SFM.surv=SFM(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv)
rownames(tableALL.SFM.surv)="SFM"
tableALL.SFM.w.surv=SFM(X.w, nT.w, DA.w, A.w, N.w, Index.w_T_HC0, Index.w_T_HC1, ATT.surv.w)
rownames(tableALL.SFM.w.surv)="SFM.w"
######## SIPTW
tableALL.SIPTW.surv=SIPTW(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv)
tableALL.SIPTW.surv
rownames(tableALL.SIPTW.surv)="SIPTW"
tableALL.SIPTW.w.surv=SIPTW(X.w, nT.w, DA.w, A.w, N.w, Index.w_T_HC0, Index.w_T_HC1, ATT.surv.w)
rownames(tableALL.SIPTW.w.surv)="SIPTW.w"
##### JFM2
tableALL.JFM2.surv=JFM2(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv, D, DD, Index_noC)
tableALL.JFM2.surv
rownames(tableALL.JFM2.surv)="JFM2"
tableALL.JFM2.w.surv=JFM2(X.w, nT.w, DA.w, A.w, N.w, Index.w_T_HC0, Index.w_T_HC1, ATT.surv.w, D.w, DD.w, Index.w_noC)
rownames(tableALL.JFM2.w.surv)="JFM2.w"
####### JIPTW2
tableALL.JIPTW2.surv=JIPTW2(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv, D, DD, Index_noC)
rownames(tableALL.JIPTW2.surv)="JIPTW2"
tableALL.JIPTW2.w.surv=JIPTW2(X.w, nT.w, DA.w, A.w, N.w, Index.w_T_HC0, Index.w_T_HC1, ATT.surv.w, D.w, DD.w, Index.w_noC)
rownames(tableALL.JIPTW2.w.surv)="JIPTW2.w"



######## Bayesian
# NPD
start.time <- Sys.time()
tableALL.NPD=NPD(nHC0, nHC1,nT, DA, A, N, ATT.surv, "NPD.surv", trt, extD)
rownames(tableALL.NPD)="NPD"
# tableALL.NPD.w=NPD(nHC0, nHC1,nT.w, DA.w, A.w, N.w, ATT.surv.w, "NPD.surv", trt.w, extD.w)
# rownames(tableALL.NPD.w)="NPD.w"
end.time <- Sys.time()
end.time - start.time
# NPS
start.time <- Sys.time()
tableALL.NPS=NPS(nT, DA, A, N, ATT.surv, "NPS.surv", trt, extS)
rownames(tableALL.NPS)="NPS"
# tableALL.NPS.w=NPS(nT.w, DA.w, A.w, N.w, ATT.surv.w, "NPS.surv", trt.w, extS.w)
# rownames(tableALL.NPS.w)="NPS.w"
end.time <- Sys.time()
end.time - start.time
# IPD
start.time <- Sys.time()
# tableALL.IPD=IPD(nHC0, nHC1,nT, DA, A, N, ATT.surv, "IPD.surv", trt, extD)
# rownames(tableALL.IPD)="IPD"
tableALL.IPD.w=IPD(nHC0, nHC1,nT.w, DA.w, A.w, N.w, ATT.surv.w, "IPD.surv", trt.w, extD.w)
rownames(tableALL.IPD.w)="IPD.w"
end.time <- Sys.time()
end.time - start.time
# IPS
start.time <- Sys.time()
# tableALL.IPS=IPS(nT, DA, A, N, ATT.surv, "IPS.surv", trt, extS)
# rownames(tableALL.IPS)="IPS"
tableALL.IPS.w=IPS(nT.w, DA.w, A.w, N.w, ATT.surv.w, "IPS.surv", trt.w, extS.w)
rownames(tableALL.IPS.w)="IPS.w"
end.time <- Sys.time()
end.time - start.time
# WPD
start.time <- Sys.time()
# tableALL.WPD=WPD(nT, DA, A, N, ATT.surv, "WPD.surv", trt, extD)
# rownames(tableALL.WPD)="WPD"
tableALL.WPD.w=WPD(nHC0, nHC1,nT.w, DA.w, A.w, N.w, ATT.surv.w, "WPD.surv", trt.w, extD.w)
rownames(tableALL.WPD.w)="WPD.w"
end.time <- Sys.time()
end.time - start.time
# WPS
start.time <- Sys.time()
# tableALL.WPS=WPS(nT, DA, A, N, ATT.surv, "WPS.surv", trt, extS)
# rownames(tableALL.WPS)="WPS"
tableALL.WPS.w=WPS(nT.w, DA.w, A.w, N.w, ATT.surv.w, "WPS.surv", trt.w, extS.w)
rownames(tableALL.WPS.w)="WPS.w"
end.time <- Sys.time()
end.time - start.time
# NB
start.time <- Sys.time()
tableALL.NB=NB(nT, DA, A, nT+nC, ATT.surv, "NB.surv", trt)
rownames(tableALL.NB)="NB"
tableALL.NB.w=NB(nT.w, DA.w, A.w, nT.w+nC, ATT.surv.w, "NB.surv", trt.w)
rownames(tableALL.NB.w)="NB.w"
end.time <- Sys.time()
end.time - start.time
# FB
start.time <- Sys.time()
tableALL.FB=FB(nT, DA, A, N, ATT.surv, "FB.surv", trt, extD)
rownames(tableALL.FB)="FB"
tableALL.FB.w=FB(nT.w, DA.w, A.w, N.w, ATT.surv.w, "FB.surv", trt.w, extD.w)
rownames(tableALL.FB.w)="FB.w"
end.time <- Sys.time()
end.time - start.time


