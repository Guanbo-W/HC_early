# Frequentist
########## AFT
##### Separate Full Matching
SFM.aft=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      #for (i in 1:Nsim){
      X=X_gen(nT,nC,nHC0,nHC1)
      Y=Gen_Y(N, DA, beta[,j], X, nT)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
      m.out_HC0=matchit(DA[Index_T_HC0,5] ~ X[Index_T_HC0,], method= "full",estimand="ATT")
      weights_HC0=m.out_HC0$weights[-c(1:nT)]
      m.out_HC1=matchit(DA[Index_T_HC1,5] ~ X[Index_T_HC1,], method = "full",estimand="ATT")
      weights_HC1=m.out_HC1$weights[-c(1:nT)]
      weights=as.vector(c(rep(1,nT+nC),weights_HC0,weights_HC1))
      M=survreg(formula = Y ~ A, data=data.frame(Y,A),weights = weights, dist="weibull")
      result1=sum(weights[(nT+nC):N])^2/sum((weights[(nT+nC):N])^2)
      result2=-M$coef[4]/M$scale
      betahat_T=c(result1,result2)}
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

######## SIPTW
SIPTW.aft=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv){
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      Y=Gen_Y(N, DA, beta[,j], X, nT)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
      weights_HC0=predict(glm(DA[Index_T_HC0,5] ~ X[Index_T_HC0,],family=binomial(link = "logit")),type="response")[-c(1:nT)]
      weights_HC1=predict(glm(DA[Index_T_HC1,5] ~ X[Index_T_HC1,],family=binomial(link = "logit")),type="response")[-c(1:nT)]
      weights=as.vector(c(rep(1,nT+nC),weights_HC0/(1-weights_HC0),weights_HC1/(1-weights_HC1)))
      M=survreg(formula = Y ~ A, data=data.frame(Y,A),weights = weights, dist="weibull")
      result1=sum(weights[(nT+nC):N])^2/sum((weights[(nT+nC):N])^2)
      result2=-M$coef[4]/M$scale
      betahat_T=c(result1,result2)}
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

##### JFM2
JFM2.aft=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv, D, DD, Index_noC){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      Y=Gen_Y(N, DA, beta[,j], X, nT)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
      m.out=matchit(c(rep(1,nT),rep(0,nHC0+nHC1)) ~ X[Index_noC,], method= "full",estimand="ATT")
      weights_HC01=m.out$weights[-c(1:nT)]
      weights=as.vector(c(rep(1,nT+nC),weights_HC01))
      M=survreg(formula = Y ~ A, data=data.frame(Y,A),weights = weights, dist="weibull")
      result1=sum(weights[(nT+nC):N])^2/sum((weights[(nT+nC):N])^2)
      result2=-M$coef[4]/M$scale
      betahat_T=c(result1,result2)}
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

####### JIPTW2
JIPTW2.aft=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv, D, DD, Index_noC){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      Y=Gen_Y(N, DA, beta[,j], X, nT)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
      weights_HC01=predict(glm(c(rep(1,nT),rep(0,nHC0+nHC1)) ~ X[Index_noC,],family=binomial(link = "logit")),type="response")[-c(1:nT)]
      weights=as.vector(c(rep(1,nT+nC),weights_HC01/(1-weights_HC01))) 
      M=survreg(formula = Y ~ A, data=data.frame(Y,A),weights = weights, dist="weibull")
      result1=sum(weights[(nT+nC):N])^2/sum((weights[(nT+nC):N])^2)
      result2=-M$coef[4]/M$scale
      betahat_T=c(result1,result2)}
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

########## CoxPH
##### Separate Full Matching
SFM=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      Y=Gen_Y(N, DA, beta[,j], X, nT)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
      m.out_HC0=matchit(DA[Index_T_HC0,5] ~ X[Index_T_HC0,], method= "full",estimand="ATT")
      weights_HC0=m.out_HC0$weights[-c(1:nT)]
      m.out_HC1=matchit(DA[Index_T_HC1,5] ~ X[Index_T_HC1,], method = "full",estimand="ATT")
      weights_HC1=m.out_HC1$weights[-c(1:nT)]
      weights=as.vector(c(rep(1,nT+nC),weights_HC0,weights_HC1))
      result1=sum(weights[(nT+nC):N])^2/sum((weights[(nT+nC):N])^2)
      result2=summary(coxph(formula = Y ~ A, weights = weights, robust = TRUE))$coef[3]
      betahat_T=c(result1,result2)}
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

######## SIPTW
SIPTW=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      Y=Gen_Y(N, DA, beta[,j], X, nT)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
      weights_HC0=predict(glm(DA[Index_T_HC0,5] ~ X[Index_T_HC0,],family=binomial(link = "logit")),type="response")[-c(1:nT)]
      weights_HC1=predict(glm(DA[Index_T_HC1,5] ~ X[Index_T_HC1,],family=binomial(link = "logit")),type="response")[-c(1:nT)]
      weights=as.vector(c(rep(1,nT+nC),weights_HC0/(1-weights_HC0),weights_HC1/(1-weights_HC1)))
      result1=sum(weights[(nT+nC):N])^2/sum((weights[(nT+nC):N])^2)
      result2=summary(coxph(formula = Y ~ A, weights = weights, robust = TRUE))$coef[3]
      betahat_T=c(result1,result2)}
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}


####### JIPTW2
JIPTW2=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT.surv, D, DD, Index_noC){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      Y=Gen_Y(N, DA, beta[,j], X, nT)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
      weights_HC01=predict(glm(c(rep(1,nT),rep(0,nHC0+nHC1)) ~ X[Index_noC,],family=binomial(link = "logit")),type="response")[-c(1:nT)]
      weights=as.vector(c(rep(1,nT+nC),weights_HC01/(1-weights_HC01))) 
      result1=sum(weights[(nT+nC):N])^2/sum((weights[(nT+nC):N])^2)
      result2=summary(coxph(formula = Y ~ A, weights = weights, robust = TRUE))$coef[3]
      betahat_T=c(result1,result2)}
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

######## Bayesian
# NPD
NPD=function(nT, DA, A, N, ATT.surv, NPD.surv, trt, extD){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      outcome=Gen_Y_Bayes(N, DA, beta[,j], X, nT)
      Y=outcome[,1]; event=outcome[,2]
      name=paste ("~/Library/Mobile Documents/com~apple~CloudDocs/ROCHE/SIM/Bayesian jags/", NPD.surv, sep = "", collapse = NULL)
      jags.out <- jags.model(name,data=list(N=N,Y=Y,event=event,trt = trt,extD = extD, nHC0=nHC0, nHC1=nHC1),
                             inits=list(alpha=c(0,0,0), beta=0, tau=0.03),
                             n.chains=n.chains,n.adapt=n.adapt, quiet=TRUE)
      results=summary(coda.samples(jags.out, "logHR_TC", n.iter=n.iter))
      result1=(results$statistics[2])^2*N
      result2=results$statistics[1]
      betahat_T=c(result1,result2)
      }
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

# NPS
NPS=function(nT, DA, A, N, ATT.surv, NPS.surv, trt, extS){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      outcome=Gen_Y_Bayes(N, DA, beta[,j], X, nT)
      Y=outcome[,1]; event=outcome[,2]
      name=paste ("~/Library/Mobile Documents/com~apple~CloudDocs/ROCHE/SIM/Bayesian jags/", NPS.surv, sep = "", collapse = NULL)
      jags.out <- jags.model(name,data=list(N=N,Y=Y,event=event,trt = trt,extS = extS, nHC0=nHC0, nHC1=nHC1),
                             inits=list(alpha=c(0,0), beta=0, tau1=0.03),
                             n.chains=n.chains,n.adapt=n.adapt, quiet=TRUE)
      results=summary(coda.samples(jags.out, "logHR_TC", n.iter=n.iter))
      result1=(results$statistics[2])^2*N
      result2=results$statistics[1]
      betahat_T=c(result1,result2)
    }
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

# IPD
IPD=function(nT, DA, A, N, ATT.surv, IPD.surv, trt, extD){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      outcome=Gen_Y_Bayes(N, DA, beta[,j], X, nT)
      Y=outcome[,1]; event=outcome[,2]
      name=paste ("~/Library/Mobile Documents/com~apple~CloudDocs/ROCHE/SIM/Bayesian jags/", IPD.surv, sep = "", collapse = NULL)
      jags.out <- jags.model(name,data=list(N=N,Y=Y,event=event,trt = trt,extD = extD, nHC0=nHC0, nHC1=nHC1),
                             inits=list(alpha=c(0,0,0), beta=0, tau1=0.03),
                             n.chains=n.chains,n.adapt=n.adapt, quiet=TRUE)
      results=summary(coda.samples(jags.out, "logHR_TC", n.iter=n.iter))
      result1=(results$statistics[2])^2*N
      result2=results$statistics[1]
      betahat_T=c(result1,result2)
    }
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

# IPS
IPS=function(nT, DA, A, N, ATT.surv, IPS.surv, trt, extS){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      outcome=Gen_Y_Bayes(N, DA, beta[,j], X, nT)
      Y=outcome[,1]; event=outcome[,2]
      name=paste ("~/Library/Mobile Documents/com~apple~CloudDocs/ROCHE/SIM/Bayesian jags/", IPS.surv, sep = "", collapse = NULL)
      jags.out <- jags.model(name,data=list(N=N,Y=Y,event=event,trt = trt,extS = extS, nHC0=nHC0, nHC1=nHC1),
                             inits=list(alpha=c(0,0), beta=0, tau1=0.03),
                             n.chains=n.chains,n.adapt=n.adapt, quiet=TRUE)
      results=summary(coda.samples(jags.out, "logHR_TC", n.iter=n.iter))
      result1=(results$statistics[2])^2*N
      result2=results$statistics[1]
      betahat_T=c(result1,result2)
    }
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

# WPD
WPD=function(nT, DA, A, N, ATT.surv, WPD.surv, trt, extD){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      outcome=Gen_Y_Bayes(N, DA, beta[,j], X, nT)
      Y=outcome[,1]; event=outcome[,2]
      name=paste ("~/Library/Mobile Documents/com~apple~CloudDocs/ROCHE/SIM/Bayesian jags/", WPD.surv, sep = "", collapse = NULL)
      jags.out <- jags.model(name,data=list(N=N,Y=Y,event=event,trt = trt,extD = extD,nHC0=nHC0,nHC1=nHC1),
                             inits=list(alpha=c(0,0,0), beta=0),
                             n.chains=n.chains,n.adapt=n.adapt, quiet=TRUE)
      results=summary(coda.samples(jags.out, "logHR_TC", n.iter=n.iter))
      result1=(results$statistics[2])^2*N
      result2=results$statistics[1]
      betahat_T=c(result1,result2)
    }
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

# WPS
WPS=function(nT, DA, A, N, ATT.surv, WPS.surv, trt, extS){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      outcome=Gen_Y_Bayes(N, DA, beta[,j], X, nT)
      Y=outcome[,1]; event=outcome[,2]
      name=paste ("~/Library/Mobile Documents/com~apple~CloudDocs/ROCHE/SIM/Bayesian jags/", WPS.surv, sep = "", collapse = NULL)
      jags.out <- jags.model(name,data=list(N=N,Y=Y,event=event,trt = trt,extS = extS, nHC0=nHC0, nHC1=nHC1),
                             inits=list(alpha=c(0,0), beta=0),
                             n.chains=n.chains,n.adapt=n.adapt, quiet=TRUE)
      results=summary(coda.samples(jags.out, "logHR_TC", n.iter=n.iter))
      result1=(results$statistics[2])^2*N
      result2=results$statistics[1]
      betahat_T=c(result1,result2)
    }
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

NB=function(nT, DA, A, N, ATT.surv, NB.surv, trt){
  for(j in 1:4){
    for (i in 1:Nsim) {
      X=X_gen(nT,nC,nHC0,nHC1)
      outcome=Gen_Y_Bayes(N, DA, beta[,j], X, nT)
      Y=outcome[,1]; event=outcome[,2]
      #name=paste ("~/Library/Mobile Documents/com~apple~CloudDocs/ROCHE/SIM/Bayesian jags/", NB.surv, sep = "", collapse = NULL)
      jags.out <- jags.model(NB.surv,data=list(N=N,Y=Y,event=event,trt = trt, nHC0=nHC0, nHC1=nHC1),
                             inits=list(alpha1=0, beta=0, r0=0.03),
                             n.chains=n.chains,n.adapt=n.adapt, quiet=TRUE)
      results=summary(coda.samples(jags.out, "logHR_TC", n.iter=n.iter))
      result1=(results$statistics[2])^2*N
      result2=results$statistics[1]
      betahat_T=c(result1,result2)
    }
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}


# FB
FB=function(nT, DA, A, N, ATT.surv, FB.surv, trt, extD){
  registerDoParallel(Ncore)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      outcome=Gen_Y_Bayes(N, DA, beta[,j], X, nT)
      Y=outcome[,1]; event=outcome[,2]
      name=paste ("~/Library/Mobile Documents/com~apple~CloudDocs/ROCHE/SIM/Bayesian jags/", FB.surv, sep = "", collapse = NULL)
      jags.out <- jags.model(name,data=list(N=N,Y=Y,event=event,trt = trt, extD=extD, nHC0=nHC0, nHC1=nHC1),
                             inits=list(r0=1, alpha1=0,beta=0),
                             n.chains=n.chains,n.adapt=n.adapt, quiet=TRUE)
      results=summary(coda.samples(jags.out, "logHR_TC", n.iter=n.iter))
      result1=(results$statistics[2])^2*N
      result2=results$statistics[1]
      betahat_T=c(result1,result2)
    }
    Bias=mean(betahat_T[,2])-ATT.surv[j]
    Variance=sum((betahat_T[,2]-mean(betahat_T[,2]))^2)/(Nsim-1)
    MSE=mean((betahat_T[,2]-ATT.surv[j])^2)
    ESS=mean(betahat_T[,1])
    table[[j]]=round(cbind(Bias,Variance,MSE,ESS),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}
