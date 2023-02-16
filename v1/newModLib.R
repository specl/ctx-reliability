library(R2jags)

# Model 1, One Task
# Model 2, Many Tasks
# Model 3, One Measure
# Model 4, Many Measures

write1=function(fileName){
  write('
 model{
  for (n in 1:N) {
   y[n]~dnorm(alpha[sub[n]]+x[n]*theta[sub[n]],pSigma)}
  for (i in 1:I){
    alpha[i]~dnorm(mu,1)
    theta[i]~dnorm(nu,pSigma*pGamma)}
  mu~dnorm(1,1)
  pSigma~dgamma(.5,.5)
  nu~dnorm(.07,1/.1^2)
  pGamma~dgamma(.5,.5*tune^2)
  }',fileName)}

write1("mod1.txt")

makeJagsDat1=function(dat,tune){
  dat$sub=as.integer(as.factor(dat$sub))
  y=dat$rt
  sub=dat$sub
  N=length(y)
  I=max(sub)
  x=dat$cond-1.5
  return(list("y"=y,"sub"=sub,"x"=x,"N"=N,"I"=I,"tune"=tune))}

makeInits1=function(dat){
  res=aov(y~as.factor(sub)+as.factor(x),data=dat)$residuals
  sig=sd(res)
  mrt=tapply(dat$y,list(dat$sub,dat$x),mean)
  inits=list(list(
    mu=mean(dat$y),
    pSigma=1/sig^2,
    alpha=tapply(dat$y,dat$sub,mean),
    theta=mrt[,2]-mrt[,1],
    nu=mean(mrt[,2]-mrt[,1]),
    pGamma=sig^2/var(mrt[,2]-mrt[,1])))
  return(inits)}

est1=function(dat,M,burn,tune=1){
  datJ=makeJagsDat1(dat,tune)
  inits=makeInits1(datJ) 
  parameters=c("alpha","theta","pSigma","pGamma")
  samples <- jags(datJ, 
                  inits=inits, 
                  parameters=parameters, 
                  model.file = "mod1.txt", 
                  n.chains=1,n.iter=M,n.burnin=burn,n.thin=1)
  alpha=samples$BUGSoutput$sims.list$alpha
  theta=samples$BUGSoutput$sims.list$theta
  pSigma=samples$BUGSoutput$sims.list$pSigma
  pGamma=samples$BUGSoutput$sims.list$pGamma
  return(list("alpha"=alpha,"theta"=theta,"pSigma"=pSigma,"pGamma"=pGamma))}


write2=function(fileName){
  write('
 model{
  for (n in 1:N) {
   y[n]~dnorm(alpha[sub[n],task[n]]+(x[n]*theta[sub[n],task[n]]),
          pSig[task[n]])}
  for (i in 1:I){
    for (j in 1:J){
      alpha[i,j]~dnorm(mu,1)}}
  for (i in 1:I){
    theta[i,(1:J)]~dmnorm(nu,B)}
  mu~dnorm(1,1)
  for (j in 1:J){
    pSig[j]~dgamma(.5,.5)
    sig2[j] <- 1/pSig[j]
    sigma[j] <- sqrt(sig2[j])
    nu[j]~dnorm(.07,1/.1^2)}
  B~dwish(W,J+1)
  }',fileName)}

write2("mod2.txt")

makeJagsDat2=function(dat,tune){
  dat$sub=as.integer(as.factor(dat$sub))
  y=dat$rt
  sub=dat$sub
  task=dat$task
  N=length(y)
  I=max(sub)
  J=max(task)
  x=dat$cond-1.5
  W=diag(J)*tune^2
  return(list("y"=y,"sub"=sub,"task"=task,"x"=x,"N"=N,"I"=I,"J"=J,"W"=W))}

makeInits2=function(datJ){
  res=aov(datJ$y~as.factor(datJ$sub)+as.factor(datJ$task)+as.factor(datJ$x))$residuals
  sig=tapply(res,datJ$task,sd)
  mrt=tapply(datJ$y,list(datJ$sub,datJ$task),mean)
  eff=(datJ$y-mrt[cbind(datJ$sub,datJ$task)])
  mrEff=tapply(eff,list(datJ$sub,datJ$task,datJ$x),mean)
    inits=list(list(
    mu=mean(datJ$y),
    pSig=1/sig^2,
    alpha=tapply(datJ$y,list(datJ$sub,datJ$task),mean),
    theta=mrEff[,,2]-mrEff[,,1],
    nu=apply(mrEff[,,2]-mrEff[,,1],2,mean),
    B=solve(cov(mrEff[,,2]-mrEff[,,1]))))
  return(inits)}

est2=function(dat,M,burn,tune){
  datJ=makeJagsDat2(dat,tune)
  inits=makeInits2(datJ) 
  parameters=c("alpha","theta","pSig","B")
  samples <- jags(datJ, 
                  inits=inits, 
                  parameters=parameters, 
                  model.file = "mod2.txt", 
                  n.chains=1,n.iter=M,n.burnin=burn,n.thin=1)
    alpha=samples$BUGSoutput$sims.list$alpha
    theta=samples$BUGSoutput$sims.list$theta
    pSig=samples$BUGSoutput$sims.list$pSig
    B=samples$BUGSoutput$sims.list$B
  return(list("alpha"=alpha,"theta"=theta,"pSig"=pSig,"B"=B))}


write3=function(fileName){
  write('
 model{
  for (n in 1:N) {
   y[n]~dnorm(theta[sub[n]],pSigma)}
  for (i in 1:I){
    theta[i]~dnorm(nu,pGamma*pSigma)}
  pSigma~dgamma(.5,.5)
  nu~dnorm(10,1/100)
  pGamma~dgamma(.5,.5*tune^2)
  }',fileName)}

write3("mod3.txt")

makeJagsDat3=function(dat,tune){
  dat$sub=as.integer(as.factor(dat$sub))
  y=dat$y
  sub=dat$sub
  N=length(y)
  I=max(sub)
  return(list("y"=y,"sub"=sub,"N"=N,"I"=I,"tune"=tune))}

makeInits3=function(dat){
  res=aov(y~as.factor(sub),data=dat)$residuals
  sig=sd(res)
  mrt=tapply(dat$y,dat$sub,mean)
  inits=list(list(
    nu=mean(dat$y),
    pSigma=1/sig^2,
    theta=mrt,
    pGamma=sig^2/var(mrt)))
  return(inits)}

est3=function(dat,M,burn,tune=1){
  datJ=makeJagsDat3(dat,tune)
  inits=makeInits3(datJ) 
  parameters=c("theta","pSigma","pGamma")
  samples <- jags(datJ, 
                  inits=inits, 
                  parameters=parameters, 
                  model.file = "mod3.txt", 
                  n.chains=1,n.iter=M,n.burnin=burn,n.thin=1)
  theta=samples$BUGSoutput$sims.list$theta
  pSigma=samples$BUGSoutput$sims.list$pSigma
  pGamma=samples$BUGSoutput$sims.list$pGamma
  return(list("theta"=theta,"pSigma"=pSigma,"pGamma"=pGamma))}



write4=function(fileName){
  write('
 model{
  for (n in 1:N) {
   y[n]~dnorm(theta[sub[n],task[n]],pSig[task[n]])}
  for (i in 1:I){
    theta[i,(1:J)]~dmnorm(nu,B)}
  for (j in 1:J){
    pSig[j]~dgamma(.5,.5)
    nu[j]~dnorm(10,1/100)}
  B~dwish(W,J+1)
  }',fileName)}

write4("mod4.txt")

makeJagsDat4=function(dat,tune){
  dat$sub=as.integer(as.factor(dat$sub))
  y=dat$y
  sub=dat$sub
  task=dat$task
  N=length(y)
  I=max(sub)
  J=max(task)
  W=diag(J)*tune^2
  return(list("y"=y,"sub"=sub,"task"=task,"N"=N,"I"=I,"J"=J,"W"=W))}

makeInits4=function(datJ){
  res=aov(datJ$y~as.factor(datJ$sub)+as.factor(datJ$task))$residuals
  sig=tapply(res,datJ$task,sd)
  m=tapply(datJ$y,list(datJ$sub,datJ$task),mean)
  inits=list(list(
    pSig=1/sig^2,
    theta=m,
    nu=apply(m,2,mean),
    B=solve(cov(m))))
  return(inits)}

est4=function(dat,M,burn,tune){
  datJ=makeJagsDat4(dat,tune)
  inits=makeInits4(datJ) 
  parameters=c("theta","pSig","B")
  samples <- jags(datJ, 
                  inits=inits, 
                  parameters=parameters, 
                  model.file = "mod4.txt", 
                  n.chains=1,n.iter=M,n.burnin=burn,n.thin=1)
  theta=samples$BUGSoutput$sims.list$theta
  pSig=samples$BUGSoutput$sims.list$pSig
  B=samples$BUGSoutput$sims.list$B
  return(list("theta"=theta,"pSig"=pSig,"B"=B))}



write5=function(fileName){
  write('
 model{
  for (n in 1:N) {
   y[n]~dnorm(alpha[sub[n]]+x[n]*theta[sub[n]],pSigma)}
  for (i in 1:I){
    alpha[i]~dnorm(nuAlpha,pSigma*pGammaAlpha)
    theta[i]~dnorm(nuTheta,pSigma*pGammaTheta)}
  nuAlpha~dnorm(1,1)
  pGammaAlpha~dgamma(.5,.5*tuneA^2)
  pSigma~dgamma(.5,.5)
  nuTheta~dnorm(.07,1/.1^2)
  pGammaTheta~dgamma(.5,.5*tuneT^2)
  }',fileName)}

write5("mod5.txt")

makeJagsDat5=function(dat,tuneA,tuneT){
  dat$sub=as.integer(as.factor(dat$sub))
  y=dat$rt
  sub=dat$sub
  N=length(y)
  I=max(sub)
  x=dat$cond-1.5
  return(list("y"=y,"sub"=sub,"x"=x,"N"=N,"I"=I,"tuneA"= tuneA, "tuneT" = tuneT))}



makeInits5=function(dat){
  res=aov(y~as.factor(sub)+as.factor(x),data=dat)$residuals
  sig=sd(res)
  mrt=tapply(dat$y,list(dat$sub,dat$x),mean)
  inits=list(list(
    nuAlpha=mean(dat$y),
    pSigma=1/sig^2,
    alpha=tapply(dat$y,dat$sub,mean),
    theta=mrt[,2]-mrt[,1],
    nuTheta=mean(mrt[,2]-mrt[,1]),
    pGammaTheta=sig^2/var(mrt[,2]-mrt[,1]),
    pGammaAlpha=sig^2/var((mrt[,2]+mrt[,1])/2)
    ))
  return(inits)
}

est5=function(dat,M,burn,tuneA=1, tuneT=1){
  datJ=makeJagsDat5(dat,tuneA, tuneT)
  inits=makeInits5(datJ) 
  parameters=c("alpha","theta","pSigma","pGammaTheta", "pGammaAlpha")
  samples <- jags(datJ, 
                  inits=inits, 
                  parameters=parameters, 
                  model.file = "mod5.txt", 
                  n.chains=1,n.iter=M,n.burnin=burn,n.thin=1)
  alpha=samples$BUGSoutput$sims.list$alpha
  theta=samples$BUGSoutput$sims.list$theta
  pSigma=samples$BUGSoutput$sims.list$pSigma
  pGammaAlpha=samples$BUGSoutput$sims.list$pGammaAlpha
  pGammaTheta=samples$BUGSoutput$sims.list$pGammaTheta
  return(list("alpha"=alpha,"theta"=theta,"pSigma"=pSigma,"pGammaTheta"=pGammaTheta, "pGammaAlpha"=pGammaAlpha))}
