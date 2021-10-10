# Binary outcome
source("~/Prep_Binary.R")
source("~/Methods_Binary.R")
table=NULL
betahat_T=NULL

## Frequentist

#### SFM
tableALL.SFM.b=SFM.b(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT)
rownames(tableALL.SFM.b)="SFM.b"
tableALL.SFM.b.w=SFM.b(X.w, nT.w, DA.w, A.w, N.w, Index.w_T_HC0, Index.w_T_HC1, ATT.w)
rownames(tableALL.SFM.b.w)="SFM.b.w"
#### SIPTW
tableALL.SIPTW.b=SIPTW.b(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT)
rownames(tableALL.SIPTW.b)="SIPTW.b"
tableALL.SIPTW.b.w=SIPTW.b(X.w, nT.w, DA.w, A.w, N.w, Index.w_T_HC0, Index.w_T_HC1, ATT.w)
rownames(tableALL.SIPTW.b.w)="SIPTW.b.w"
#### JFM1
tableALL.SFM1.b=SFM1.b(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC)
rownames(tableALL.SFM1.b)="JFM1.b"
tableALL.SFM1.b.w=SFM1.b(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC)
rownames(tableALL.SFM1.b.w)="JFM1.b.w"
#### JIPTW1
tableALL.SIPTW1.b=SIPTW1.b(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC)
rownames(tableALL.SIPTW1.b)="SIPTW1.b"
tableALL.SIPTW1.b.w=SIPTW1.b(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC)
rownames(tableALL.SIPTW1.b.w)="SIPTW1.b.w"
#### JFM2
tableALL.SFM2.b=SFM2.b(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC)
rownames(tableALL.SFM2.b)="JFM2.b"
tableALL.SFM2.b.w=SFM2.b(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC)
rownames(tableALL.SFM2.b.w)="JFM2.b.w"
#### JIPTW2
tableALL.SIPTW2.b=SIPTW2.b(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC)
rownames(tableALL.SIPTW2.b)="SIPTW2.b"
tableALL.SIPTW2.b.w=SIPTW2.b(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC)
rownames(tableALL.SIPTW2.b.w)="SIPTW2.b.w"
############ Bayesian
#### IPD
Sys.time()
tableALL.IPD.b=IPD.b(nT, DA, A, N, ATT, "NPD", trt, extD)
rownames(tableALL.IPD.b)="IPD.b"
tableALL.IPD.b.w=IPD.b(nT.w, DA.w, A.w, N.w, ATT.w, "NPD", trt.w, extD.w)
rownames(tableALL.IPD.b.w)="IPD.b.w"
Sys.time()
#### IPS
Sys.time()
tableALL.IPS.b=IPS.b(nT, DA, A, N, ATT, "NPS", trt, extS)
rownames(tableALL.IPS.b)="IPS.b"
tableALL.IPS.b.w=IPS.b(nT.w, DA.w, A.w, N.w, ATT.w, "NPS", trt.w, extS.w)
rownames(tableALL.IPS.b.w)="IPS.b.w"
Sys.time()
#### NPD
Sys.time()
tableALL.NPD.b=NPD.b(nT, DA, A, N, ATT, "NNPD", trt, extD)
rownames(tableALL.NPD.b)="NPD.b"
tableALL.NPD.b.w=NPD.b(nT.w, DA.w, A.w, N.w, ATT.w, "NNPD", trt.w, extD.w)
rownames(tableALL.NPD.b.w)="NPD.b.w"
Sys.time()
#### NPS
Sys.time()
tableALL.NPS.b=NPS.b(nT, DA, A, N, ATT, "NNPS", trt, extS)
rownames(tableALL.NPS.b)="NPS.b"
tableALL.NPS.b.w=NPS.b(nT.w, DA.w, A.w, N.w, ATT.w, "NNPS", trt.w, extS.w)
rownames(tableALL.NPS.b.w)="NPS.b.w"
Sys.time()
#### WPD
Sys.time()
tableALL.WPD.b=WPD.b(nT, DA, A, N, ATT, "WPD", trt, extD)
rownames(tableALL.WPD.b)="WPD.b"
tableALL.WPD.b.w=WPD.b(nT.w, DA.w, A.w, N.w, ATT.w, "WPD", trt.w, extD.w)
rownames(tableALL.WPD.b.w)="WPD.b.w"
Sys.time()
#### WPS
Sys.time()
tableALL.WPS.b=WPS.b(nT, DA, A, N, ATT, "WPS", trt, extS)
rownames(tableALL.WPS.b)="WPS.b"
tableALL.WPS.b.w=NPS.b(nT.w, DA.w, A.w, N.w, ATT.w, "WPS", trt.w, extS.w)
rownames(tableALL.WPS.b.w)="WPS.b.w"
Sys.time()
#### NB
Sys.time()
tableALL.NB.binary=NB(nT, DA, A, nT+nC, ATT, "NB", trt)
rownames(tableALL.NB.b)="NB.b"
tableALL.NB.b.w=NB(nT.w, DA.w, A.w, nT.w+nC, ATT.w, "NB", trt.w)
rownames(tableALL.NB.b.w)="NB.b.w"
Sys.time()
#### FB
Sys.time()
tableALL.FB.binary=FB(nT, DA, A, N, ATT, "FB", trt, extD)
rownames(tableALL.FB.b)="FB.b"
tableALL.FB.w.binary=FB(nT.w, DA.w, A.w, N.w, ATT.w, "FB", trt.w, extD.w)
rownames(tableALL.FB.b.w)="FB.b.w"
Sys.time()
################### Results
index_bias=c(1,4,7,10)
index_variance=c(2,5,8,11)
index_MSE=c(3,6,9,12)
# BIAS
Methods_grp=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM.l","SIPTW.l","JFM1.l","JIPTW1.l","JFM2.l","JIPTW2.l"),4)
Methods=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2"),4)
MT=rep(c(rep("16",6),rep("30",6)),4)
Scenario=c(rep("S1",12),rep("S2",12),rep("S3",12),rep("S4",12))

data_bias_freq=abs(rbind(tableALL.SFM.b,tableALL.SIPTW.b,tableALL.JFM1.b,tableALL.JIPTW1.b,tableALL.JFM2.b,tableALL.JIPTW2.b,tableALL.SFM.b.w,tableALL.SIPTW.b.w,tableALL.JFM1.b.w,tableALL.JIPTW1.b.w,tableALL.JFM2.b.w,tableALL.JIPTW2.b.w)[,index_bias])
Bias=c(data_bias_freq[,1],data_bias_freq[,2],data_bias_freq[,3],data_bias_freq[,4])
data_bias_freq=data.frame(Bias,Scenario,MT,Methods,Methods_grp)

gg1=ggplot(data_bias_freq, aes(x=Scenario, y=Bias, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  geom_point(aes(color=Methods))+
  ggtitle("Biases of Frequentist methods (Binary)") +
  theme(plot.title = element_text(hjust = 0.5))  + ylim(0,6.5)

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
  ggtitle("Biases of Bayesian methods (Binary)") +
  theme(plot.title = element_text(hjust = 0.5))   +ylim(0,6.5)

grid.arrange(gg1, gg2, ncol=2)

# VARIANCE
Methods_grp=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM.l","SIPTW.l","JFM1.l","JIPTW1.l","JFM2.l","JIPTW2.l"),4)
Methods=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2"),4)
MT=rep(c(rep("16",6),rep("30",6)),4)
Scenario=c(rep("S1",12),rep("S2",12),rep("S3",12),rep("S4",12))

data_v_freq=abs(rbind(tableALL.SFM.b,tableALL.SIPTW.b,tableALL.JFM1.b,tableALL.JIPTW1.b,tableALL.JFM2.b,tableALL.JIPTW2.b,tableALL.SFM.b.w,tableALL.SIPTW.b.w,tableALL.JFM1.b.w,tableALL.JIPTW1.b.w,tableALL.JFM2.b.w,tableALL.JIPTW2.b.w)[,index_variance])
Variance=c(data_v_freq[,1],data_v_freq[,2],data_v_freq[,3],data_v_freq[,4])
data_v_freq=data.frame(Variance,Scenario,MT,Methods,Methods_grp)

gg3=ggplot(data_v_freq, aes(x=Scenario, y=Variance, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  geom_point(aes(color=Methods))+
  ggtitle("Variances of Frequentist methods (Binary)") +
  theme(plot.title = element_text(hjust = 0.5))   +ylim(0,130)

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
  ggtitle("Variances of Bayesian methods (Binary)") +
  theme(plot.title = element_text(hjust = 0.5))   + ylim(0,130)

grid.arrange(gg3, gg4, ncol=2)


# MSE
Methods_grp=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM.l","SIPTW.l","JFM1.l","JIPTW1.l","JFM2.l","JIPTW2.l"),4)
Methods=rep(c("SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2","SFM","SIPTW","JFM1","JIPTW1","JFM2","JIPTW2"),4)
MT=rep(c(rep("16",6),rep("30",6)),4)
Scenario=c(rep("S1",12),rep("S2",12),rep("S3",12),rep("S4",12))

data_mse_freq=abs(rbind(tableALL.SFM.b,tableALL.SIPTW.b,tableALL.JFM1.b,tableALL.JIPTW1.b,tableALL.JFM2.b,tableALL.JIPTW2.b,tableALL.SFM.b.w,tableALL.SIPTW.b.w,tableALL.JFM1.b.w,tableALL.JIPTW1.b.w,tableALL.JFM2.b.w,tableALL.JIPTW2.b.w)[,index_MSE])
MSE=c(data_mse_freq[,1],data_mse_freq[,2],data_mse_freq[,3],data_mse_freq[,4])
data_mse_freq=data.frame(MSE,Scenario,MT,Methods,Methods_grp)

gg5=ggplot(data_mse_freq, aes(x=Scenario, y=MSE, group=Methods_grp)) +
  geom_line(aes(linetype=MT,color=Methods))+
  geom_point(aes(color=Methods))+
  ggtitle("MSEs of Frequentist methods (Binary)") +
  theme(plot.title = element_text(hjust = 0.5))  + ylim(0,155)

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
  ggtitle("MSEs of Bayesian methods (Binary)") +
  theme(plot.title = element_text(hjust = 0.5))  + ylim(0,155)

grid.arrange(gg5, gg6, ncol=2)

