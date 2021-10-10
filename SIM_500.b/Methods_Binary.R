## Frequentist

#### SFM
SFM.b=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha +X[,8]*DA[,5] * gamma)
      Y=rbinom(N, size=1, prob=pi)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
      m.out_HC0=matchit(DA[Index_T_HC0,5] ~ X[Index_T_HC0,], method= "full",estimand="ATT")
      weights_HC0=m.out_HC0$weights[-c(1:nT)]
      m.out_HC1=matchit(DA[Index_T_HC1,5] ~ X[Index_T_HC1,], method = "full",estimand="ATT")
      weights_HC1=m.out_HC1$weights[-c(1:nT)]
      weights=as.vector(c(rep(1,nT+nC),weights_HC0,weights_HC1))
      model=glm(Y~A,family = quasibinomial(link = "logit"), weights=weights)
      betahat_T[i]=model$coef[4] }
  Bias=mean(betahat_T)-ATT[j]
  Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
  MSE=mean((betahat_T-ATT[j])^2)
  table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

#### SIPTW
SIPTW.b=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha +X[,8]*DA[,5] * gamma)
      Y=rbinom(N, size=1, prob=pi)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
    weights_HC0=predict(glm(DA[Index_T_HC0,5] ~ X[Index_T_HC0,],family=binomial(link = "logit")),type="response")[-c(1:nT)]
    weights_HC1=predict(glm(DA[Index_T_HC1,5] ~ X[Index_T_HC1,],family=binomial(link = "logit")),type="response")[-c(1:nT)]
    weights=as.vector(c(rep(1,nT+nC),weights_HC0/(1-weights_HC0),weights_HC1/(1-weights_HC1)))
    model=glm(Y~A,family = quasibinomial(link = "logit"), weights=weights)
    betahat_T[i]=model$coef[4] }
    Bias=mean(betahat_T)-ATT[j]
    Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
    MSE=mean((betahat_T-ATT[j])^2)
    table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

#### JFM1
JFM1.b=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha +X[,8]*DA[,5] * gamma)
      Y=rbinom(N, size=1, prob=pi)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
    m.out=matchit(c(rep(1,nT),rep(0,nHC0+nHC1)) ~ X[Index_noC,] + DD[Index_noC,2] + DD[Index_noC,3], method= "full",estimand="ATT")
    weights_HC01=m.out$weights[-c(1:nT)]
    weights=as.vector(c(rep(1,nT+nC),weights_HC01))
    model=glm(Y~A,family = quasibinomial(link = "logit"), weights=weights)
    betahat_T[i]=model$coef[4] }
    Bias=mean(betahat_T)-ATT[j]
    Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
    MSE=mean((betahat_T-ATT[j])^2)
    table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

#### JIPTW 1 
JIPTW1.b=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha +X[,8]*DA[,5] * gamma)
      Y=rbinom(N, size=1, prob=pi)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
    weights_HC01=predict(glm(c(rep(1,nT),rep(0,nHC0+nHC1)) ~ X[Index_noC,]  + DD[Index_noC,2] + DD[Index_noC,3],family=binomial(link = "logit")),type="response")[-c(1:nT)]
    weights=as.vector(c(rep(1,nT+nC),weights_HC01/(1-weights_HC01)))
    model=glm(Y~A,family = quasibinomial(link = "logit"), weights=weights)
    betahat_T[i]=model$coef[4] }
    Bias=mean(betahat_T)-ATT[j]
    Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
    MSE=mean((betahat_T-ATT[j])^2)
    table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

#### JFM2 
JFM2.b=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha +X[,8]*DA[,5] * gamma)
      Y=rbinom(N, size=1, prob=pi)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
    m.out=matchit(c(rep(1,nT),rep(0,nHC0+nHC1)) ~ X[Index_noC,], method= "full",estimand="ATT")
    weights_HC01=m.out$weights[-c(1:nT)]
    weights=as.vector(c(rep(1,nT+nC),weights_HC01))
    model=glm(Y~A + DD[,2] + DD[,3],family = quasibinomial(link = "logit"), weights=weights)
    betahat_T[i]=model$coef[4] }
    Bias=mean(betahat_T)-ATT[j]
    Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
    MSE=mean((betahat_T-ATT[j])^2)
    table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

#### JIPTW 2 
JIPTW2.b=function(X, nT, DA, A, N, Index_T_HC0, Index_T_HC1, ATT, DD, Index_noC){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha +X[,8]*DA[,5] * gamma)
      Y=rbinom(N, size=1, prob=pi)
      X=X_gen(nT,nC,nHC0,nHC1)[,-c(11:13)]
    weights_HC01=predict(glm(c(rep(1,nT),rep(0,nHC0+nHC1)) ~ X[Index_noC,],family=binomial(link = "logit")),type="response")[-c(1:nT)]
    weights=as.vector(c(rep(1,nT+nC),weights_HC01/(1-weights_HC01)))
    model=glm(Y~A+DD[,2]+DD[,3],family = quasibinomial(link = "logit"), weights=weights)
    betahat_T[i]=model$coef[4] }
  Bias=mean(betahat_T)-ATT[j]
  Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
  MSE=mean((betahat_T-ATT[j])^2)
  table[[j]]=round(cbind(Bias,Variance,MSE),3)}
table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
return(table)}

#### Bayesian
#### IPD
IPD.b=function(nT, DA, A, N, ATT, NPD, trt, extD){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
    X=X_gen(nT,nC,nHC0,nHC1)
    pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha+X[,8]*DA[,5] * gamma)
    Y=rbinom(N, size=1, prob=pi)
    name=paste ("~/jags_Binary/", NPD, sep = "", collapse = NULL)
    jags.out <- jags.model(name,
                           data=list(N=N,Y=Y,trt = trt,extD = extD),
                           inits=list(alpha=c(0,0,0), beta=0, tau1=0.03),
                           n.chains=3,n.adapt=1000, quiet=TRUE)
    betahat_T[i]=summary(coda.samples(jags.out, "logOR_TC", n.iter=2000))$statistics[1]}
    Bias=mean(betahat_T)-ATT[j]
    Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
    MSE=mean((betahat_T-ATT[j])^2)
    table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

#### IPS 
IPS.b=function(nT, DA, A, N, ATT, NPS, trt, extS){
registerDoParallel(16)
for(j in 1:4){
  betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
    X=X_gen(nT,nC,nHC0,nHC1)
    pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha+X[,8]*DA[,5] * gamma)
    Y=rbinom(N, size=1, prob=pi)
    name=paste ("~/jags_Binary/", NPS, sep = "", collapse = NULL)
    jags.out <- jags.model(name,
                           data=list(N=N,Y=Y,trt = trt,extS = extS),
                           inits=list(alpha=c(0,0), beta=0, tau1=0.03),
                           n.chains=3,n.adapt=1000, quiet=TRUE)
    betahat_T[i]=summary(coda.samples(jags.out, "logOR_TC", n.iter=2000))$statistics[1]}
  Bias=mean(betahat_T)-ATT[j]
  Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
  MSE=mean((betahat_T-ATT[j])^2)
  table[[j]]=round(cbind(Bias,Variance,MSE),3)}
table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
return(table)}

#### NPD
NPD.b=function(nT, DA, A, N, ATT.surv, NNPD, trt, extD){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha+X[,8]*DA[,5] * gamma)
      Y=rbinom(N, size=1, prob=pi)
      name=paste ("~/jags_Binary/", NNPD, sep = "", collapse = NULL)
      jags.out <- jags.model(name,
                           data=list(N=N,Y=Y,trt = trt,extD = extD),
                           inits=list(alpha=c(0,0,0), beta=0, tau1=0.03),
                           n.chains=3,n.adapt=1000, quiet=TRUE)
    betahat_T[i]=summary(coda.samples(jags.out, "logOR_TC", n.iter=2000))$statistics[1]}
  Bias=mean(betahat_T)-ATT[j]
  Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
  MSE=mean((betahat_T-ATT[j])^2)
  table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

#### NPS
NPS.b=function(nT, DA, A, N, ATT, NNPS, trt, extS){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
    X=X_gen(nT,nC,nHC0,nHC1)
    pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha+X[,8]*DA[,5] * gamma)
    Y=rbinom(N, size=1, prob=pi)
    name=paste ("~/jags_Binary/", NNPS, sep = "", collapse = NULL)
    jags.out <- jags.model(name,data=list(N=N,Y=Y,trt = trt,extS = extS),
                           inits=list(alpha=c(0,0), beta=0, tau1=0.03),
                           n.chains=3,n.adapt=1000, quiet=TRUE)
    betahat_T[i]=summary(coda.samples(jags.out, "logOR_TC", n.iter=2000))$statistics[1]}
  Bias=mean(betahat_T)-ATT[j]
  Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
  MSE=mean((betahat_T-ATT[j])^2)
  table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

#### WPD
WPD.b=function(nT, DA, A, N, ATT.surv, WPD, trt, extD){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha+X[,8]*DA[,5] * gamma)
      Y=rbinom(N, size=1, prob=pi)
      name=paste ("~/jags_Binary/", WPD, sep = "", collapse = NULL)
      jags.out <- jags.model(name,
                             data=list(N=N,Y=Y,trt = trt,extD = extD),
                             inits=list(alpha=c(0,0,0), beta=0, tau1=0.03),
                             n.chains=3,n.adapt=1000, quiet=TRUE)
      betahat_T[i]=summary(coda.samples(jags.out, "logOR_TC", n.iter=2000))$statistics[1]}
    Bias=mean(betahat_T)-ATT[j]
    Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
    MSE=mean((betahat_T-ATT[j])^2)
    table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

#### WPS
WPS.b=function(nT, DA, A, N, ATT, WPS, trt, extS){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha+X[,8]*DA[,5] * gamma)
      Y=rbinom(N, size=1, prob=pi)
      name=paste ("~/jags_Binary/", WPS, sep = "", collapse = NULL)
      jags.out <- jags.model(name,data=list(N=N,Y=Y,trt = trt,extS = extS),
                             inits=list(alpha=c(0,0), beta=0, tau1=0.03),
                             n.chains=3,n.adapt=1000, quiet=TRUE)
      betahat_T[i]=summary(coda.samples(jags.out, "logOR_TC", n.iter=2000))$statistics[1]}
    Bias=mean(betahat_T)-ATT[j]
    Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
    MSE=mean((betahat_T-ATT[j])^2)
    table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}

#### NB
NB.b=function(nT, DA, A, N, ATT, NB, trt){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha+X[,8]*DA[,5] * gamma)
      Y=rbinom(N, size=1, prob=pi)
      name=paste ("~/jags_Binary/", NB, sep = "", collapse = NULL)
      jags.out <- jags.model(name,data=list(N=N,Y=Y,trt = trt),
                             inits=list(alpha=c(0,0,0), beta=0, tau1=0.03),
                             n.chains=3,n.adapt=1000, quiet=TRUE)
      betahat_T[i]=summary(coda.samples(jags.out, "logOR_TC", n.iter=2000))$statistics[1]}
    Bias=mean(betahat_T)-ATT[j]
    Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
    MSE=mean((betahat_T-ATT[j])^2)
    table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}


##### FB
FB.b=function(nT, DA, A, N, ATT, FB, trt, extD){
  registerDoParallel(16)
  for(j in 1:4){
    betahat_T=foreach (i=1:Nsim, .combine=rbind) %dopar% {
      X=X_gen(nT,nC,nHC0,nHC1)
      pi <- plogis(beta0 + as.matrix(DA[,3:5]) %*% beta[,j] + X %*% alpha+X[,8]*DA[,5] * gamma)
      Y=rbinom(N, size=1, prob=pi)
      name=paste ("~/jags_Binary/", FB, sep = "", collapse = NULL)
      jags.out <- jags.model(name,data=list(N=N,Y=Y,trt = trt, extD=extD),
                             inits=list(alpha=c(0,0,0), beta=0, tau1=0.03),
                             n.chains=3,n.adapt=1000, quiet=TRUE)
      betahat_T[i]=summary(coda.samples(jags.out, "logOR_TC", n.iter=2000))$statistics[1]}
    Bias=mean(betahat_T)-ATT[j]
    Variance=sum((betahat_T-mean(betahat_T))^2)/(Nsim-1)
    MSE=mean((betahat_T-ATT[j])^2)
    table[[j]]=round(cbind(Bias,Variance,MSE),3)}
  table=cbind(table[[1]],table[[2]],table[[3]],table[[4]])
  return(table)}




