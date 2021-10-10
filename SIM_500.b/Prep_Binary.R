library(dummies)
library(MatchIt)
library(optmatch)
library(doParallel)
library(foreach)
library(fastDummies)
library(coda)
library(rjags)
library(ggplot2)
require(gridExtra)
library(survival)
set.seed(2021)

Nsim=500

############ Sample size
nT=16; nT.w=30; nC=22; nHC0=145; nHC1=467 
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
# beta (log OR/log HR)=c(HC0, HC1, T)
# Scenario 1 (All equivalent)
beta1=c(log(1), log(1), log(1))

# Scenario 2 (HC0=HC1=C<T)
beta2=c(log(1), log(1), log(2))

# Scenario 3 (HC0=HC1<C<T)
beta3=c(log(0.7), log(0.7), log(2))

# Scenario 4  (HC0<HC1<C<T)
beta4=c(log(0.5), log(0.7), log(2))

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
gamma=0.05

############ beta0
beta0=0.2

############ Estimate true ATT: marginal OR(T/C)

### ORR samll MT
# samplenumber=200
# result=NULL
# for(i in 1:10000){
#   X=X_gen(nT,nC,nHC0,nHC1)
#   index=sample(seq(1,nT,1),samplenumber,replace=TRUE)
#   pi0=plogis(beta0 + X[index,]%*% alpha + X[index, 8] *DA[index, 5] * gamma)
#   pi1=plogis(beta0 + log(2)*rep(1,samplenumber) + X[index,] %*% alpha + X[index, 8] * gamma)
#   Y0=rbinom(samplenumber, size=1, prob=pi0)
#   Y1=rbinom(samplenumber, size=1, prob=pi1)
#   new.Y=c(Y0,Y1);new.trt=c(rep(0,samplenumber),rep(1,samplenumber))
#   result[i]=glm(new.Y~new.trt,family = binomial)$coef[2]
# }
# result=mean(result)
result=0.7031033
ATT=c(0,result,result,result)
# 
# ### ORR large MT
# samplenumber=200
# result=NULL
# for(i in 1:10000){
#   X.w=X_gen(nT.w,nC,nHC0,nHC1)
#   index=sample(seq(1,nT.w,1),samplenumber,replace=TRUE)
#   pi0=plogis(beta0 + X.w[index,]%*% alpha + X.w[index, 8] *DA[index, 5]* gamma)
#   pi1=plogis(beta0 + log(2)*rep(1,samplenumber) + X.w[index,] %*% alpha + X.w[index, 8]*DA[index, 5] * gamma)
#   Y0=rbinom(samplenumber, size=1, prob=pi0)
#   Y1=rbinom(samplenumber, size=1, prob=pi1)
#   new.Y=c(Y0,Y1);new.trt=c(rep(0,samplenumber),rep(1,samplenumber))
#   result[i]=glm(new.Y~new.trt,family = binomial)$coef[2]
# }
# result=mean(result)
result=0.7038467 
ATT.w=c(0,result,result,result)

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