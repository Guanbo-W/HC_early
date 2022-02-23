# Evaluating hybrid controls methodology in the metastatic urothelial carcinoma Morpheus trial
# Phase 1b/II uMC Morpheous Trail. Borrowing two historical control
# Survival Outcome

######
############ load library
library(dummies)
library(MatchIt)
library(optmatch)
library(doParallel)
library(foreach)
library(fastDummies)
library(tinytex)
library(coda)
library(rjags)
library(ggplot2)
require(gridExtra)
library(survival)
set.seed(2021)

Nsim=500

###### 
############ Sample size
nT=16; nT.w=40; nC=22; nHC0=145; nHC1=467 
N=nT+nC+nHC0+nHC1
N.w=nT.w+nC+nHC0+nHC1
############ Arm
A_MT=rep("T",nT)
A_MT.w=rep("T.w",nT.w)
A_MC=rep("C",nC)
A_HC0=rep("HC0",nHC0)
A_HC1=rep("HC1",nHC1)
A=relevel(factor(c(A_MT, A_MC, A_HC0, A_HC1)),ref="C")
DA=dummy_cols(A)
A.w=relevel(factor(c(A_MT.w, A_MC, A_HC0, A_HC1)),ref="C")
DA.w=dummy_cols(A.w)

############ Conditional Arm effect 
# beta (log HR)=c(HC0, HC1, T)
# Scenario 1 (All equivalent)
beta1=c(log(1), log(1), log(1))

# New Scenario 2 (HC0=HC1=C<T)
beta2=c(log(1), log(1), log(1/2))

# New Scenario 3 (HC0=HC1<C<T)
beta3=c(log(3), log(3), log(1/2))

# New Scenario 4  (HC0<HC1<C<T)
beta4=c(log(12), log(3), log(1/2))

beta=cbind(beta1,beta2,beta3,beta4)

############ 13 Potential baseline confounders: 2 continuous, 11 binary
# Age
X_gen=function(nT,nC,nHC0,nHC1){
  X1_T=round(rnorm(nT, mean=65, sd=6))
  X1_C=round(rnorm(nC, mean=70, sd=8))
  X1_HC0=round(rnorm(nHC0, mean=65, sd=10))
  X1_HC1=round(rnorm(nHC1, mean=67, sd=10))
  X1=c(X1_T, X1_C, X1_HC0, X1_HC1)
  
  # BMI
  X2_T=rnorm(nT, mean=24.91, sd=6)
  X2_C=rnorm(nC, mean=25.36, sd=6)
  X2_HC0=rnorm(nHC0, mean=26.87, sd=9)
  X2_HC1=rnorm(nHC1, mean=25.84, sd=8)
  X2=c(X2_T, X2_C, X2_HC0, X2_HC1)
  
  # Sex (Female)
  X3_T=rbinom(nT, 1, 0.25)
  X3_C=rbinom(nC, 1, 0.14)
  X3_HC0=rbinom(nHC0, 1, 0.19)
  X3_HC1=rbinom(nHC1, 1, 0.24)
  X3=c(X3_T, X3_C, X3_HC0, X3_HC1)
  
  # ECOG (1)
  X4_T=rbinom(nT, 1, 0.69)
  X4_C=rbinom(nC, 1, 0.60)
  X4_HC0=rbinom(nHC0, 1, 0.57)
  X4_HC1=rbinom(nHC1, 1, 0.53)
  X4=c(X4_T, X4_C, X4_HC0, X4_HC1)
  
  # Race (White)
  X5_T=rbinom(nT, 1, 0.56)
  X5_C=rbinom(nC, 1, 0.64)
  X5_HC0=rbinom(nHC0, 1, 0.95)
  X5_HC1=rbinom(nHC1, 1, 0.71)
  X5=c(X5_T, X5_C, X5_HC0, X5_HC1)
  
  # Smoking status (Yes)
  X6_T=rbinom(nT, 1, 0.75)
  X6_C=rbinom(nC, 1, 0.82)
  X6_HC0=rbinom(nHC0, 1, 0.61)
  X6_HC1=rbinom(nHC1, 1, 0.70)
  X6=c(X6_T, X6_C, X6_HC0, X6_HC1)
  
  # Surgery (Yes)
  X7_T=rbinom(nT, 1, 0.94)
  X7_C=rbinom(nC, 1, 0.90)
  X7_HC0=rbinom(nHC0, 1, 0.87)
  X7_HC1=rbinom(nHC1, 1, 0.94)
  X7=c(X7_T, X7_C, X7_HC0, X7_HC1)
  
  # PDL1 (Positive)
  X8_T=rbinom(nT, 1, 0.38)
  X8_C=rbinom(nC, 1, 0.55)
  X8_HC0=rbinom(nHC0, 1, 0.31)
  X8_HC1=rbinom(nHC1, 1, 0.25)
  X8=c(X8_T, X8_C, X8_HC0, X8_HC1)
  
  # Alkaline Phosphatase (High)
  X9_T=rbinom(nT, 1, 0.25)
  X9_C=rbinom(nC, 1, 0.18)
  X9_HC0=rbinom(nHC0, 1, 0.27)
  X9_HC1=rbinom(nHC1, 1, 0.24)
  X9=c(X9_T, X9_C, X9_HC0, X9_HC1)
  
  # CRP (High)
  X10_T=rbinom(nT, 1, 0.75)
  X10_C=rbinom(nC, 1, 0.68)
  X10_HC0=rbinom(nHC0, 1, 0.81)
  X10_HC1=rbinom(nHC1, 1, 0.54)
  X10=c(X10_T, X10_C, X10_HC0, X10_HC1)
  
  # Lymph node only site of metastasis (Yes)
  X11_T=rbinom(nT, 1, 0.19)
  X11_C=rbinom(nC, 1, 0.14)
  X11_HC0=rbinom(nHC0, 1, 0.18)
  X11_HC1=rbinom(nHC1, 1, 0.12)
  X11=c(X11_T, X11_C, X11_HC0, X11_HC1)
  
  # Visceral site of metastasis (Yes)
  X12_T=rbinom(nT, 1, 0.69)
  X12_C=rbinom(nC, 1, 0.84)
  X12_HC0=rbinom(nHC0, 1, 0.72)
  X12_HC1=rbinom(nHC1, 1, 0.77)
  X12=c(X12_T, X12_C, X12_HC0, X12_HC1)
  
  # Liver site of metastasis (Yes)
  X13_T=rbinom(nT, 1, 0.06)
  X13_C=rbinom(nC, 1, 0.14)
  X13_HC0=rbinom(nHC0, 1, 0.23)
  X13_HC1=rbinom(nHC1, 1, 0.30)
  X13=c(X13_T, X13_C, X13_HC0, X13_HC1)
  
  X=cbind(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13)
}
X=X_gen(nT,nC,nHC0,nHC1)
X.w=X_gen(nT.w,nC,nHC0,nHC1)


############ Confounder coefficient
alpha=seq(0.01,0.07,0.005)
############ EM coefficient
gamma=0.01
############ beta0
beta0=0.2

######
############ Estimate true ATT: marginal OR(T/C)
Gen_Y=function(N, arm, con_effect, cov, nT){
  v=2
  A=-log(runif(N,0,1))
  B= exp(as.matrix(arm[,3:5]) %*% con_effect + cov %*% alpha + cov[, 8]*arm[,5]*gamma)
  Y=(A/B)^(1/v)
  Y1=Y[1:nT];Y2=Y[(nT+1):(nT+nC)];Y3=Y[(nT+nC+1):(nT+nC+nHC0)];Y4=Y[(nT+nC+nHC0+1):(nT+nC+nHC0+nHC1)]
  C1=rnorm(nT,mean(1.5*Y1),0.1);c1=ifelse(C1>Y1,1,0)
  C2=rnorm(nC,mean(1.5*Y2),0.1);c2=ifelse(C2>Y2,1,0)
  C3=rnorm(nHC0,mean(1.5*Y3),0.1);c3=ifelse(C3>Y3,1,0)
  C4=rnorm(nHC1,mean(1.5*Y4),0.1);c4=ifelse(C4>Y4,1,0)
  return(Surv(Y,c(c1,c2,c3,c4)))
}
# ### Generating survival outcomes
# Gen_Y=function(N, arm, con_effect, cov, nT){
#   v=2
#   A=-log(runif(N,0,1))
#   B= exp(as.matrix(arm[,3:5]) %*% con_effect + cov %*% alpha + cov[, 8]*arm[,5]*gamma)
#   Y=(A/B)^(1/v)
#   Y1=Y[1:nT];Y2=Y[(nT+1):(nT+nC)];Y3=Y[(nT+nC+1):(nT+nC+nHC0)];Y4=Y[(nT+nC+nHC0+1):(nT+nC+nHC0+nHC1)]
#   C1=rnorm(nT,mean(Y1),0.1);c1=ifelse(C1>Y1,1,0)
#   C2=rnorm(nC,mean(Y2),0.1);c2=ifelse(C2>Y2,1,0)
#   C3=rnorm(nHC0,mean(Y3),0.1);c3=ifelse(C3>Y3,1,0)
#   C4=rnorm(nHC1,mean(Y4),0.1);c4=ifelse(C4>Y4,1,0)
#   return(Surv(Y,c(c1,c2,c3,c4)))
# }

Gen_Y_Bayes=function(N, arm, con_effect, cov, nT){
  v=2
  A=-log(runif(N,0,1))
  B= exp(as.matrix(arm[,3:5]) %*% con_effect + cov %*% alpha + cov[, 8]*arm[,5]*gamma)
  Y=(A/B)^(1/v)
  Y1=Y[1:nT];Y2=Y[(nT+1):(nT+nC)];Y3=Y[(nT+nC+1):(nT+nC+nHC0)];Y4=Y[(nT+nC+nHC0+1):(nT+nC+nHC0+nHC1)]
  C1=rnorm(nT,mean(Y1),0.1);c1=ifelse(C1>Y1,1,0)
  C2=rnorm(nC,mean(Y2),0.1);c2=ifelse(C2>Y2,1,0)
  C3=rnorm(nHC0,mean(Y3),0.1);c3=ifelse(C3>Y3,1,0)
  C4=rnorm(nHC1,mean(Y4),0.1);c4=ifelse(C4>Y4,1,0)
  return(cbind(Y,c(c1,c2,c3,c4)))
}

# #hist(Gen_Y(N, DA, beta[,2],X))
#
#
# # survival small MT ATT

# Gen_Y_ATT=function(AA, X){
#   v=2
#   A=-log(runif(length(AA),0,1))
#   B=exp(as.matrix(AA*log(1/2) + X %*% alpha + X[, 8]*AA*gamma))
#   Time=(A/B)^(1/v)
#   return(Time)}
# 
# samplenumber=200
# result=NULL
# registerDoParallel(16)
#  result=foreach (i=1:100000, .combine=rbind) %dopar% {
#   X=X_gen(nT,nC,nHC0,nHC1)
#   index=sample(seq(1,(nT+nC),1),samplenumber,replace=TRUE)
#   Y0=Gen_Y_ATT(rep(0,samplenumber),X[index,])
#   Y1=Gen_Y_ATT(rep(1,samplenumber),X[index,])
#   new.Y=c(Y0,Y1);new.trt=c(rep(0,samplenumber),rep(1,samplenumber))
#   result[i]=summary(coxph(Surv(new.Y)~new.trt))$coef[1]
# }
# result=mean(result)
# # #
result=-0.6810789
ATT.surv=c(0,result,result,result)
#
# survival big MT
# samplenumber=200
# result=NULL
# registerDoParallel(16)
# result=foreach (i=1:100000, .combine=rbind) %dopar% {
#   X.w=X_gen(nT.w,nC,nHC0,nHC1)
#   index=sample(seq(1,(nT.w+nC),1),samplenumber,replace=TRUE)
#   Y0=Gen_Y_ATT(rep(0,samplenumber),X.w[index,])
#   Y1=Gen_Y_ATT(rep(1,samplenumber),X.w[index,])
#   new.Y=c(Y0,Y1);new.trt=c(rep(0,samplenumber),rep(1,samplenumber))
#   result[i]=summary(coxph(Surv(new.Y)~new.trt))$coef[1]
# }
# result=mean(result)
result=-0.6807652
ATT.surv.w=c(0,result,result,result)


#######
Index_T_HC0=c(1:nT,(nT+nC+1):(nT+nC+nHC0))
Index_T_HC1=c(1:nT,(nT+nC+nHC0+1):N)
Index.w_T_HC0=c(1:nT.w,(nT.w+nC+1):(nT.w+nC+nHC0))
Index.w_T_HC1=c(1:nT.w,(nT.w+nC+nHC0+1):N.w)
D=c(rep("TrialM",nT+nC),rep("TrialHC0",nHC0),rep("TrialHC1",nHC1))
D.w=c(rep("TrialM",nT.w+nC),rep("TrialHC0",nHC0),rep("TrialHC1",nHC1))
DD=as.matrix(dummy_cols(D))
DD.w=as.matrix(dummy_cols(D.w))
Index_noC=c(1:nT,(nT+nC+1):N)
Index.w_noC=c(1:nT.w,(nT.w+nC+1):N.w)

extD=c(rep(1,nT+nC),rep(2,nHC0),rep(3,nHC1))
extD.w=c(rep(1,nT.w+nC),rep(2,nHC0),rep(3,nHC1))
extS=c(rep(1,nT+nC),rep(2,nHC0+nHC1))
extS.w=c(rep(1,nT.w+nC),rep(2,nHC0+nHC1))
trt=c(rep(1,nT),rep(0,N-nT))
trt.w=c(rep(1,nT.w),rep(0,N.w-nT.w))

