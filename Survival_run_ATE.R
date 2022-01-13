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



################### Results
index_bias=c(1,4,7,10)
index_variance=c(2,5,8,11)
index_MSE=c(3,6,9,12)
# BIAS
Methods_grp=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM.l","SIPTW.l","JFM1.l","JIPTW1.l","JFM2.l","JIPTW2.l"),4)
Methods=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2"),4)
MT=rep(c(rep("16",6),rep("30",6)),4)
Scenario=c(rep("S1",12),rep("S2",12),rep("S3",12),rep("S4",12))

data_bias_freq.aft=abs(rbind(tableALL.SFM.surv.aft,tableALL.SIPTW.surv.aft,tableALL.JFM1.surv.aft,tableALL.JIPTW1.surv.aft,tableALL.JFM2.surv.aft,tableALL.JIPTW2.surv.aft,tableALL.SFM.w.surv.aft,tableALL.SIPTW.w.surv.aft,tableALL.JFM1.w.surv.aft,tableALL.JIPTW1.w.surv.aft,tableALL.JFM2.w.surv.aft,tableALL.JIPTW2.w.surv.aft)[,index_bias])
Bias=c(data_bias_freq.aft[,1],data_bias_freq.aft[,2],data_bias_freq.aft[,3],data_bias_freq.aft[,4])
data_bias_freq.aft=data.frame(Bias,Scenario,MT,Methods,Methods_grp)

gg11=ggplot(data_bias_freq.aft, aes(x=Scenario, y=Bias, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  geom_point(aes(color=Methods))+
  ggtitle("Biases of Frequentist methods (AFT)") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1.1)

Methods_grp=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM.l","SIPTW.l","JFM1.l","JIPTW1.l","JFM2.l","JIPTW2.l"),4)
Methods=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2"),4)
MT=rep(c(rep("16",6),rep("30",6)),4)
Scenario=c(rep("S1",12),rep("S2",12),rep("S3",12),rep("S4",12))

data_bias_freq=abs(rbind(tableALL.SFM.surv,tableALL.SIPTW.surv,tableALL.JFM1.surv,tableALL.JIPTW1.surv,tableALL.JFM2.surv,tableALL.JIPTW2.surv,tableALL.SFM.w.surv,tableALL.SIPTW.w.surv,tableALL.JFM1.w.surv,tableALL.JIPTW1.w.surv,tableALL.JFM2.w.surv,tableALL.JIPTW2.w.surv)[,index_bias])
Bias=c(data_bias_freq[,1],data_bias_freq[,2],data_bias_freq[,3],data_bias_freq[,4])
data_bias_freq=data.frame(Bias,Scenario,MT,Methods,Methods_grp)

gg1=ggplot(data_bias_freq, aes(x=Scenario, y=Bias, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  geom_point(aes(color=Methods))+
  ggtitle("Biases of Frequentist methods (Cox)") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1.1)

Methods_grp=rep(c("IPD","IPS","NPD","NPS","WPD","WPS", "NB", "FB", "IPD.l","IPS.l","NPD.l","NPS.l","WPD.l","WPS.l", "NB.l", "FB.l"),4)
Methods=rep(c("IPD","IPS","NPD","NPS","WPD","WPS","NB", "FB", "IPD","IPS","NPD","NPS","WPD","WPS", "NB", "FB"),4)
MT=rep(c(rep("16",8),rep("30",8)),4)
Scenario=c(rep("S1",16),rep("S2",16),rep("S3",16),rep("S4",16))

data_bias_Bayes=abs(rbind(tableALL.IPD,tableALL.IPS,tableALL.NPD,tableALL.NPS,tableALL.WPD,tableALL.WPS,tableALL.NB,tableALL.FB,tableALL.IPD.w,tableALL.IPS.w,tableALL.NPD.w,tableALL.NPS.w,tableALL.WPD,tableALL.WPS,tableALL.NB.w,tableALL.FB.w)[,index_bias])
Bias=c(data_bias_Bayes[,1],data_bias_Bayes[,2],data_bias_Bayes[,3],data_bias_Bayes[,4])
data_bias_Bayes=data.frame(Bias,Scenario,MT,Methods,Methods_grp)

gg2=ggplot(data_bias_Bayes, aes(x=Scenario, y=Bias, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  scale_colour_manual(values=c(IPD="#000066",IPS="#663399",NPD="#339999",NPS="#CC0033",WPD="#00FF33",WPS="#FFFF33",NB="#FF6600",FB="#FF9933"))+
  geom_point(aes(color=Methods))+
  ggtitle("Biases of Bayesian methods (Survival)") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,1.1)

grid.arrange(gg11, gg1, gg2, ncol=3)

# VARIANCE
Methods_grp=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM.l","SIPTW.l","JFM1.l","JIPTW1.l","JFM2.l","JIPTW2.l"),4)
Methods=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2"),4)
MT=rep(c(rep("16",6),rep("30",6)),4)
Scenario=c(rep("S1",12),rep("S2",12),rep("S3",12),rep("S4",12))

data_v_freq.aft=abs(rbind(tableALL.SFM.surv.aft,tableALL.SIPTW.surv.aft,tableALL.JFM1.surv.aft,tableALL.JIPTW1.surv.aft,tableALL.JFM2.surv.aft,tableALL.JIPTW2.surv.aft,tableALL.SFM.w.surv.aft,tableALL.SIPTW.w.surv.aft,tableALL.JFM1.w.surv.aft,tableALL.JIPTW1.w.surv.aft,tableALL.JFM2.w.surv.aft,tableALL.JIPTW2.w.surv.aft)[,index_variance])
Variance=c(data_v_freq.aft[,1],data_v_freq.aft[,2],data_v_freq.aft[,3],data_v_freq.aft[,4])
data_v_freq.aft=data.frame(Variance,Scenario,MT,Methods,Methods_grp)

gg33=ggplot(data_v_freq.aft, aes(x=Scenario, y=Variance, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  geom_point(aes(color=Methods))+
  ggtitle("Variances of Frequentist methods (AFT)") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,0.16)

Methods_grp=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM.l","SIPTW.l","JFM1.l","JIPTW1.l","JFM2.l","JIPTW2.l"),4)
Methods=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2"),4)
MT=rep(c(rep("16",6),rep("30",6)),4)
Scenario=c(rep("S1",12),rep("S2",12),rep("S3",12),rep("S4",12))

data_v_freq=abs(rbind(tableALL.SFM.surv,tableALL.SIPTW.surv,tableALL.JFM1.surv,tableALL.JIPTW1.surv,tableALL.JFM2.surv,tableALL.JIPTW2.surv,tableALL.SFM.w.surv,tableALL.SIPTW.w.surv,tableALL.JFM1.w.surv,tableALL.JIPTW1.w.surv,tableALL.JFM2.w.surv,tableALL.JIPTW2.w.surv)[,index_variance])
Variance=c(data_v_freq[,1],data_v_freq[,2],data_v_freq[,3],data_v_freq[,4])
data_v_freq=data.frame(Variance,Scenario,MT,Methods,Methods_grp)

gg3=ggplot(data_v_freq, aes(x=Scenario, y=Variance, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  geom_point(aes(color=Methods))+
  ggtitle("Variances of Frequentist methods (Cox)") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,0.16)

Methods_grp=rep(c("IPD","IPS","NPD","NPS","WPD","WPS", "NB", "FB", "IPD.l","IPS.l","NPD.l","NPS.l","WPD.l","WPS.l", "NB.l", "FB.l"),4)
Methods=rep(c("IPD","IPS","NPD","NPS","WPD","WPS","NB", "FB", "IPD","IPS","NPD","NPS","WPD","WPS", "NB", "FB"),4)
MT=rep(c(rep("16",8),rep("30",8)),4)
Scenario=c(rep("S1",16),rep("S2",16),rep("S3",16),rep("S4",16))

data_v_Bayes=abs(rbind(tableALL.IPD,tableALL.IPS,tableALL.NPD,tableALL.NPS,tableALL.WPD,tableALL.WPS,tableALL.NB,tableALL.FB,tableALL.IPD.w,tableALL.IPS.w,tableALL.NPD.w,tableALL.NPS.w,tableALL.WPD,tableALL.WPS,tableALL.NB.w,tableALL.FB.w)[,index_variance])
Variance=c(data_v_Bayes[,1],data_v_Bayes[,2],data_v_Bayes[,3],data_v_Bayes[,4])
data_v_Bayes=data.frame(Variance,Scenario,MT,Methods,Methods_grp)

gg4=ggplot(data_v_Bayes, aes(x=Scenario, y=Variance, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  geom_point(aes(color=Methods))+
  scale_colour_manual(values=c(IPD="#000066",IPS="#663399",NPD="#339999",NPS="#CC0033",WPD="#00FF33",WPS="#FFFF33",NB="#FF6600",FB="#FF9933"))+
  ggtitle("Variances of Bayesian methods (Survival)") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,0.16)

grid.arrange(gg33, gg3, gg4, ncol=3)


# MSE
Methods_grp=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM.l","SIPTW.l","JFM1.l","JIPTW1.l","JFM2.l","JIPTW2.l"),4)
Methods=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2"),4)
MT=rep(c(rep("16",6),rep("30",6)),4)
Scenario=c(rep("S1",12),rep("S2",12),rep("S3",12),rep("S4",12))

data_mse_freq.aft=abs(rbind(tableALL.SFM.surv.aft,tableALL.SIPTW.surv.aft,tableALL.JFM1.surv.aft,tableALL.JIPTW1.surv.aft,tableALL.JFM2.surv.aft,tableALL.JIPTW2.surv.aft,tableALL.SFM.w.surv.aft,tableALL.SIPTW.w.surv.aft,tableALL.JFM1.w.surv.aft,tableALL.JIPTW1.w.surv.aft,tableALL.JFM2.w.surv.aft,tableALL.JIPTW2.w.surv.aft)[,index_MSE])
MSE=c(data_mse_freq.aft[,1],data_mse_freq.aft[,2],data_mse_freq.aft[,3],data_mse_freq.aft[,4])
data_mse_freq.aft=data.frame(MSE,Scenario,MT,Methods,Methods_grp)

gg55=ggplot(data_mse_freq.aft, aes(x=Scenario, y=MSE, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  geom_point(aes(color=Methods))+
  ggtitle("MSEs of Frequentist methods (AFT)") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,0.375)

Methods_grp=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM.l","SIPTW.l","JFM1.l","JIPTW1.l","JFM2.l","JIPTW2.l"),4)
Methods=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2"),4)
MT=rep(c(rep("16",6),rep("30",6)),4)
Scenario=c(rep("S1",12),rep("S2",12),rep("S3",12),rep("S4",12))

data_mse_freq=abs(rbind(tableALL.SFM.surv,tableALL.SIPTW.surv,tableALL.JFM1.surv,tableALL.JIPTW1.surv,tableALL.JFM2.surv,tableALL.JIPTW2.surv,tableALL.SFM.w.surv,tableALL.SIPTW.w.surv,tableALL.JFM1.w.surv,tableALL.JIPTW1.w.surv,tableALL.JFM2.w.surv,tableALL.JIPTW2.w.surv)[,index_MSE])
MSE=c(data_mse_freq[,1],data_mse_freq[,2],data_mse_freq[,3],data_mse_freq[,4])
data_mse_freq=data.frame(MSE,Scenario,MT,Methods,Methods_grp)

gg5=ggplot(data_mse_freq, aes(x=Scenario, y=MSE, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  geom_point(aes(color=Methods))+
  ggtitle("MSEs of Frequentist methods (Cox)") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,0.375)

Methods_grp=rep(c("IPD","IPS","NPD","NPS","WPD","WPS", "NB", "FB", "IPD.l","IPS.l","NPD.l","NPS.l","WPD.l","WPS.l", "NB.l", "FB.l"),4)
Methods=rep(c("IPD","IPS","NPD","NPS","WPD","WPS","NB", "FB", "IPD","IPS","NPD","NPS","WPD","WPS", "NB", "FB"),4)
MT=rep(c(rep("16",8),rep("30",8)),4)
Scenario=c(rep("S1",16),rep("S2",16),rep("S3",16),rep("S4",16))

data_mse_Bayes=abs(rbind(tableALL.IPD,tableALL.IPS,tableALL.NPD,tableALL.NPS,tableALL.WPD,tableALL.WPS,tableALL.NB,tableALL.FB,tableALL.IPD.w,tableALL.IPS.w,tableALL.NPD.w,tableALL.NPS.w,tableALL.WPD,tableALL.WPS,tableALL.NB.w,tableALL.FB.w)[,index_MSE])
MSE=c(data_mse_Bayes[,1],data_mse_Bayes[,2],data_mse_Bayes[,3],data_mse_Bayes[,4])
data_mse_Bayes=data.frame(MSE,Scenario,MT,Methods,Methods_grp)

gg6=ggplot(data_mse_Bayes, aes(x=Scenario, y=MSE, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  geom_point(aes(color=Methods))+
  scale_colour_manual(values=c(IPD="#000066",IPS="#663399",NPD="#339999",NPS="#CC0033",WPD="#00FF33",WPS="#FFFF33",NB="#FF6600",FB="#FF9933"))+
  ggtitle("MSEs of Bayesian methods (Survival)") +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(0,0.375)


grid.arrange(gg55, gg5, gg6, ncol=3) 

