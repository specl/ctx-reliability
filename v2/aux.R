rmCleanData=function(dat,
                     top=1.995,                 
                     hiTrlRT=1.995,
                     loTrlRT=.275,
                     loTop=.98,
                     loAcc=.90,
                     loRT=.99){
  excludeSub=c(132,164,368,402,429)  #origial author exclusions
  #exclusions
  bad0=dat$sub %in% excludeSub
  bad1=dat$block == 'practice'
  bad2=dat$trialType != 'exp'
  bad3=!(dat$acc %in% c(0,1)) #there are some strange things
  dat=dat[!(bad0 | bad1 | bad2 | bad3),]      
  sub=tapply(dat$sub,dat$sub,mean)
  topBySub=tapply(dat$rt<top,dat$sub,mean)
  accBySub=(tapply(dat$acc,dat$sub,mean))
  loBySub=tapply(dat$rt>loTrlRT,dat$sub,mean) 
  bad4=dat$sub %in% sub[
    accBySub<loAcc |
      topBySub<loTop |
      loBySub<loRT]
  bad5 = dat$rt<loTrlRT | dat$rt>min(hiTrlRT,top) | dat$acc==0 
  dat=dat[!(bad4 | bad5),]    
  return(dat)
}

readRMStroopI=function(){
  numStroop=read.table("https://raw.githubusercontent.com/PerceptionCognitionLab/data0/master/inhibitionTasks/ReyMermetJEPLMC2018/merged/numStroop.dat",head=T)
  a3=rmCleanData(numStroop)
  numStroop=a3[a3$acc==1,]
  numStroop$cond=as.integer(as.factor(numStroop$cond))
  dat=numStroop[numStroop$cond %in% 1:2,]
  dat$y = dat$rt
  return(dat)}

readRMFlankI=function(){
  letFlanker=read.table("https://raw.githubusercontent.com/PerceptionCognitionLab/data0/master/inhibitionTasks/ReyMermetJEPLMC2018/merged/letFlanker.dat",head=T)
  a3=rmCleanData(letFlanker)
  letFlanker=a3[a3$acc==1,]
  letFlanker$cond=as.integer(as.factor(letFlanker$cond))
  dat=letFlanker[letFlanker$cond %in% 1:2,]
  dat$y = dat$rt
  return(dat)}

readRMStroopII=function(){
  colStroop=read.table("https://raw.githubusercontent.com/PerceptionCognitionLab/data0/master/inhibitionTasks/ReyMermetJEPLMC2018/merged/colStroop.dat",head=T)
  a3=rmCleanData(colStroop)
  colStroop=a3[a3$acc==1,]
  colStroop$cond=as.integer(as.factor(colStroop$cond))
  dat=colStroop[colStroop$cond %in% 1:2,]
  dat$y = dat$rt
  return(dat)}



plot.eff=function(out,dat,runner=""){
  dat$sub=as.integer(as.factor(dat$sub))
  N=length(dat$sub)
  I=length(unique(dat$sub))
  L=round(N/(2*I),0)
  place=function(x,range) range[1]+x*(range[2]-range[1])
  lower=apply(out$theta,2,quantile,p=.025)
  upper=apply(out$theta,2,quantile,p=.975)
  theta=apply(out$theta,2,mean)
  m=tapply(dat$y,list(dat$sub,dat$cond),mean)
  eff=(m[,2]-m[,1])
  I=length(eff)
  o=order(theta)
  range=c(min(c(lower,eff)),max(c(upper,eff)))
  plot(1:I,eff[o],typ='n',xlab="Individuals",
       ylab=expression(paste("Effect, ",theta)),
       ylim=range,axes=F)
  axis(2)
  axis(1,at=c(1,I))
  polygon(c(1:I,I:1),c(lower[o],rev(upper[o])),
          col='lightblue',border=NA)
  lines(1:I,theta[o],lwd=2,col='darkblue')
  lines(1:I,eff[o],lty=2,lwd=2)
  modGamma=mean(sqrt(1/out$pGamma))
  gamma=round(modGamma,3)
  text(1,place(1,range),runner,adj=c(0,1))
  text(1,place(.8,range),labels=bquote(gamma==.(gamma)),adj=c(0,1))
  text(I,place(0,range),labels=bquote(L==.(L)),adj=c(1,0))
  #observed Gamma
  err2=(dat$rt-m[cbind(dat$sub,dat$cond)])^2
  S2=sum(err2)/(N-(2*I))
  V=var(eff)
  obsGamma=sqrt(V/S2-2/L)
  return(list("obsGamma"=obsGamma,"modGamma"=modGamma))
}


plot.meas=function(out,dat,runner=""){
  dat$sub=as.integer(as.factor(dat$sub))
  N=length(dat$sub)
  I=length(unique(dat$sub))
  L=round(N/I,0)
  place=function(x,range) range[1]+x*(range[2]-range[1])
  lower=apply(out$theta,2,quantile,p=.025)
  upper=apply(out$theta,2,quantile,p=.975)
  theta=apply(out$theta,2,mean)
  m=tapply(dat$y,dat$sub,mean)
  o=order(theta)
  range=c(min(c(lower,m)),max(c(upper,m)))
  plot(1:I,m[o],typ='n',xlab="Individuals",
       ylab=expression(paste("Effect, ",theta)),
       ylim=range,axes=F)
  axis(2)
  axis(1,at=c(1,I))
  polygon(c(1:I,I:1),c(lower[o],rev(upper[o])),
          col='lightblue',border=NA)
  lines(1:I,theta[o],lwd=2,col='darkblue')
  lines(1:I,m[o],lty=2,lwd=2)
  modGamma=mean(sqrt(1/out$pGamma))
  gamma=round(modGamma,3)
  text(1,place(1,range),runner,adj=c(0,1))
  text(1,place(.8,range),labels=bquote(gamma==.(gamma)),adj=c(0,1))
  text(I,place(0,range),labels=bquote(L==.(L)),adj=c(1,0))
  err2=(dat$y-m[dat$sub])^2
  S2=sum(err2)/(N-I)
  V=var(m)
  obsGamma=sqrt(V/S2-1/L)
  return(list("obsGamma"=obsGamma,"modGamma"=modGamma))
  
}



plot.cor=function(out,dat,a=1,b=2,runner="",yObs=3.7){
  if (length(dat$cond)==0){
    dat$cond=rep(1,length(dat$sub))}
  M=dim(out$B)[1]
  corr=1:M
  for (m in 1:M) corr[m]=cov2cor(solve(out$B[m,,]))[a,b]
  hist(corr,prob=T,main="",xlim=c(-1,1),ylim=c(0,4),
       xlab="Correlation Coefficient",col='lightblue',
       breaks=seq(-1,1,.1),border='grey')
  corVal=seq(-1.,1,.01)
  lines(c(-1,corVal,1),c(0,dunif(corVal,-1,1),0))
  modCor=c(quantile(corr,p=c(.025,.975)),mean(corr))
  names(modCor)=c("Q025","Q975","mean")
  m=tapply(dat$y,list(dat$sub,dat$task,dat$cond),mean)
  if (dim(m)[3]==1) eff=m[,,1]
  if (dim(m)[3]==2) eff=m[,,2]-m[,,1]
  conv=cor.test(eff[,1],eff[,2])
  obsCor=c(conv$conf.int,conv$estimate)
  text(-1,.5,adj=c(0,-.2),"Prior")
  text(-1,4,runner,adj=c(0,1))
  arrows(obsCor[1],yObs,obsCor[2],yObs,
         code=3,angle=90,length = .1)
  points(obsCor[3],yObs,cex=1.2,pch=21,bg='yellow')
  abline(lty=2,v=modCor[1:2])
  return(rbind(modCor,obsCor))
}

readIllusion = function(){
  rm(list=ls())
  link="https://raw.githubusercontent.com/PerceptionCognitionLab/data4/master/iBat1/merged.csv"
  alldat=read.csv(link)
  alldat$sub=as.numeric(as.factor(alldat$pid))
  badSub=c(3,33,12,51,57,10,11,45)
  badPz=c(62,69,38,97,19,99,49)
  bad=c(badSub)
  alldat=alldat[!alldat$sub%in%bad,]
  goodCol=c("sub","resp","ml_len")
  dat.ml=alldat[alldat$task=='ml',goodCol]
  goodCol=c("sub","targ","resp")
  dat.pog=alldat[alldat$task=='pog',goodCol]
  
  union(unique(dat.ml$sub),unique(dat.pog$sub)) == intersect(unique(dat.ml$sub),unique(dat.pog$sub))
  
  dat.ml$bias = dat.ml$resp/dat.ml$ml_len
  dat.pog$bias = dat.pog$resp-dat.pog$targ
  dat.ml$sub=as.integer(as.factor(dat.ml$sub))
  dat.pog$sub=as.integer(as.factor(dat.pog$sub))
  dat.pog$task = "pog"
  dat.ml$task = "ml"
  dat.pog$trial = rep(c(1:length(dat.ml[dat.ml$sub==1,]$sub)),length(unique(dat.ml$sub)))
  dat.ml$trial = rep(c(1:length(dat.ml[dat.ml$sub==1,]$sub)),length(unique(dat.ml$sub)))
  dat.pog = dat.pog[,c("sub","bias","task","trial")]
  dat.ml = dat.ml[,c("sub","bias","task","trial")]
  
  tdat = rbind(dat.pog,dat.ml)
  tdat$y = tdat$bias
  # print(x)
  # if (x == T){
  #   tdat = array(NA, dim = c(nrow(dat.pog)/nrow(dat.pog[dat.pog$sub==1,]),nrow(dat.pog[dat.pog$sub==1,]),2))
  #   tdat[,,1] = matrix(dat.ml$bias, dim(tdat)[1], dim(tdat)[2], byrow = T)
  #   tdat[,,2] = matrix(dat.pog$bias, dim(tdat)[1], dim(tdat)[2], byrow = T)
  # }
  return(tdat)
}


gMCMC = function(tdat, n.chains = 2000, prio){
  library(MCMCpack)
  I = ncol(tdat) # n trial
  J = nrow(tdat) # n partip
  y_bar = apply(tdat, 1, mean)
  M = n.chains
  mu = matrix(NA, M, J)
  mu[1,] = y_bar
  nu = 1:M
  nu[1] = mean(y_bar)
  s2 = 1:M
  s2[1] = mean(apply(tdat, 1, var))
  g2 = 1:M
  g2[1] = var(y_bar)
  
  
  h = prio[1]
  k2 = prio[2]
  a = prio[3]
  b = prio[4]
  e = prio[5]
  p = prio[6]
  
  
  for (m in 2:M){
    
    delta = s2[m-1]*g2[m-1]
    v = 1/(I/s2[m-1] + 1/delta)
    c = I*y_bar/s2[m-1] + nu[m-1]/delta
    x = rnorm(J, c*v, sqrt(v))
    mu[m,] = x
    # mu[m,] = samp.theta
    
    mu_bar = mean(x)
    v = 1/(J/delta + 1/k2)
    c = J*mu_bar/delta + h/k2
    nu[m] = rnorm(1, c*v, sqrt(v))
    # nu[m] = samp.nu
    
    SSE.1 = matrix(rep(mu[m,],I),J,I)
    SSE.1 = sum((tdat-SSE.1)^2)
    SSE.2 = sum((mu[m,] - nu[m])^2)
    
    ap = J/2*(1+I) + a
    bp = SSE.1/2 + SSE.2/(2*g2[m-1])  + b
    s2[m]=rinvgamma(1, shape=ap, scale=bp)
    # s2[m] = samp.var
    
    
    
    ap = J/2 + e
    bp = SSE.2/(2*s2[m-1])  + p
    g2[m]=rinvgamma(1, shape=ap, scale=bp)
    # delta[m] = samp.g2
  }
  return(list(theta = mu, nu = nu, s2 = s2, g2 = g2))
}

plotMCMC = function(para, t.para, n.chains = 2000, name.para = "" ){
  plot(1:length(para), para , typ='l', xlab = "Iteration", main = paste("Samples of", name.para, "from MCMC, M = ", n.chains))
  abline(h=t.para,col='red',lwd=2)
  hist(para, main = paste("Samples of", name.para, "from MCMC, M = ", n.chains))
  abline(v = t.para, col = "red")
}

makePiFig = function(){
  # png('percep-ill-task.png', width = 4, height = 4, units = 'in', res = 300)
  centerX=seq(.25,.75,length=2)
  centerY=seq(.5 ,.5,length=2)
  par(mar=c(0,0,0,0))
  plot(0:1,0:1,typ='n',axes=F)
  
  text(centerX[1],centerY[1]+.25,"Mueller-Lyer", cex = 1, adj = 1)
  rect(centerX-.2,centerY-.2,centerX+.2,centerY+.2,col='antiquewhite')
  segments(centerX[1]-.07,centerY[1],centerX[1]+.09,centerY[1])
  segments(centerX[1]-.09,centerY[1]+.03,centerX[1]-.07,centerY[1])
  segments(centerX[1]-.09,centerY[1]-.03,centerX[1]-.07,centerY[1])
  segments(centerX[1]+.07,centerY[1]+.03,centerX[1]+.09,centerY[1])
  segments(centerX[1]+.07,centerY[1]-.03,centerX[1]+.09,centerY[1])
  segments(centerX[1]+.03,centerY[1]+.03,centerX[1]+.008,centerY[1], col = "blue")
  segments(centerX[1]+.03,centerY[1]-.03,centerX[1]+.008,centerY[1], col = "blue")
  segments(centerX[1]-.07,centerY[1]-.25,centerX[1]+.09,centerY[1]-.25, lty = 2)
  segments(centerX[1]+.03,centerY[1]-.22,centerX[1]+.008,centerY[1]-.25, col = "blue")
  segments(centerX[1]+.03,centerY[1]-.28,centerX[1]+.008,centerY[1]-.25, col = "blue")
  text(centerX[1]+.13,centerY[1]-.25,"+", cex = 1, col = "darkred")
  text(centerX[1]-.10,centerY[1]-.25,"-", cex = 1, col = "darkred")
  # segments(centerX[1]+.02,centerY[1]+.03,centerX[1],centerY[1])
  # segments(centerX[1]+.02,centerY[1]-.03,centerX[1],centerY[1])
  
  
  text(centerX[2],centerY[2]+.25,"Poggendorf", cex = 1, adj = 1)
  segments(centerX[2]-.02,centerY[2]-.12,centerX[2]-.02,centerY[2]+.12)
  segments(centerX[2]+.02,centerY[2]-.12,centerX[2]+.02,centerY[2]+.12)
  segments(centerX[2]-.06,centerY[2]-.08,centerX[2]-.02,centerY[2]-.02)
  # segments(centerX[2]+.06,centerY[2]+.08,centerX[2]+.02,centerY[2]+.02)
  segments(centerX[2]+.06,centerY[2]+.1,centerX[2]+.02,centerY[2]+.04, col = "blue")
  segments(centerX[2]+.22,centerY[2]-.13,centerX[2]+.22,centerY[2]+.13, lty = 2)
  segments(centerX[2]+.26,centerY[2]+.1,centerX[2]+.22,centerY[2]+.04, col = "blue")
  text(centerX[2]+.22,centerY[2]+.17,"+", cex = 1, col = "darkred")
  text(centerX[2]+.22,centerY[2]-.17,"-", cex = 1, col = "darkred")
}

readRouderOtherI=function(){
  indat=read.table(url('https://raw.githubusercontent.com/PerceptionCognitionLab/data0/master/lexDec-dist5/ld5.all'))
  colnames(indat)=c('sub','block','trial','stim','resp','rt','error')
  bad1=indat$sub%in%c(34,43)
  bad2=indat$rt<250 | indat$rt>2000
  bad3=indat$err==1
  bad4=indat$block==0 & indat$trial<20
  bad5=indat$trial==0
  bad=bad1 | bad2 | bad3 |bad4 |bad5
  tmp=indat[!bad,]
  cond=rep(0,length(tmp$stim))
  cond[tmp$stim==0 | tmp$stim==5]=1
  cond[tmp$stim==2 | tmp$stim==3]=2
  dat=data.frame(tmp$sub[cond>0],cond[cond>0],tmp$rt[cond>0])
  colnames(dat)=c('sub','cond','rt')
  dat$rt=dat$rt/1000
  dat$y = dat$rt
  return(dat)}


sampEstG = function(tdat, contrast = F){
  if (contrast == F){
    
    nt = nrow(tdat)/length(unique(tdat$sub))
    res = aov(y~as.factor(sub),data=tdat)$residuals
    sig = sd(res)
    mrt = tapply(tdat$y,list(tdat$sub),mean)
    delta = var(mrt)
    g2_bar = delta/sig^2 - 1/nt
    g_bar = sqrt(g2_bar)
    
    return(g_bar)
  }
  if (contrast == T){
    
    nt = nrow(tdat)/length(unique(tdat$sub))
    res = aov(y~as.factor(sub)+as.factor(cond),data=tdat)$residuals
    sig = sd(res)
    mrt = tapply(tdat$y,list(tdat$sub, tdat$cond),mean)
    delta.a = var((mrt[,2]+mrt[,1])/2)
    delta.t = var(mrt[,2]-mrt[,1])
    g2_bar_theta = delta.t/sig^2 - 2/(nt/2)
    g_bar_theta = sqrt(g2_bar_theta)
    g2_bar_alpha = delta.a/sig^2 - 1/nt
    g_bar_alpha = sqrt(g2_bar_alpha)
    return(list(g_bar_alpha = g_bar_alpha, g_bar_theta = g_bar_theta))
  }
}

makeArmyWeight = function(){
  tdat1 = read.csv("ANSUR_II_FEMALE_Public.csv")
  tdat2 = read.csv("ANSUR_II_MALE_Public.csv")
  dat = rbind(tdat1,tdat2)
  weight=(dat$weightkg/10)*2.2
  return(weight)
}

lSizeForR = function(r, g2, cond = 1){return(cond/(g2/r-g2))}


makeTab = function(){
  contrast = c("N", "Y", "N", "Y", "N", "Y", "N", "N", "N")
  source('newModLib.R')
  source('aux.R')
  
  gammas_sample = 1:9
  gammas_model = 1:9
  l.r9 = 1:9
  l.r7 = 1:9
  
  ## We need the army gamma model_based estimates
  gammas_model[1] = NA
  datArmy = makeArmyWeight()
  gammas_sample[1] = sd(datArmy)/3
  
  stroop <- readRMStroopI()
  flank <- readRMFlankI()
  task=rep(1:2,c(length(stroop$sub),length(flank$sub)))
  dat <-rbind(stroop,flank)
  dat$task <-task
  ssub <- unique(stroop$sub)
  fsub <- unique(flank$sub)
  goodSub <- intersect(ssub,fsub)
  datRM <- dat[dat$sub %in% goodSub,]
  outFlank=est5(datRM[datRM$task==2,],M=3000,burn=200,tuneT=.16, tuneA = 1)
  outStroop=est5(datRM[datRM$task==1,],M=3000,burn=200,tuneT=.16, tuneA = 1)
  
  
  gammas_model[2] = mean(1/sqrt(outStroop$pGammaTheta))
  gammas_model[3] = mean(1/sqrt(outStroop$pGammaAlpha))
  gammas_model[4] = mean(1/sqrt(outFlank$pGammaTheta))
  gammas_model[5] = mean(1/sqrt(outFlank$pGammaAlpha))
  
  sampOutFlannk = sampEstG(datRM[datRM$task==2,], contrast = T)
  sampOutStroop = sampEstG(datRM[datRM$task==1,], contrast = T)  
  gammas_sample[2] = sampOutStroop$g_bar_theta
  gammas_sample[3] = sampOutStroop$g_bar_alpha
  gammas_sample[4] = sampOutFlannk$g_bar_theta
  gammas_sample[5] = sampOutFlannk$g_bar_alpha
  
  
  datFive = readRouderOtherI()
  outFive=est5(datFive,M=3000,burn=200,tuneT=.16, tuneA = 1)
  
  gammas_model[6] = mean(1/sqrt(outFive$pGammaTheta))
  gammas_model[7] = mean(1/sqrt(outFive$pGammaAlpha))
  
  
  sampOutFive = sampEstG(datFive, contrast = T)
  gammas_sample[6] = sampOutFive$g_bar_theta
  gammas_sample[7] = sampOutFive$g_bar_alpha
  
  datIll <- readIllusion()
  datIll$task = as.integer(as.factor(datIll$task)) 
  datIll$y= ifelse(datIll$task==1,-datIll$bias*100,datIll$bias)
  outML=est3(datIll[datIll$task==1,],M=3000,burn=200,tune=1)
  outPog=est3(datIll[datIll$task==2,],M=3000,burn=200,tune=1)
  
  gammas_model[8] = mean(1/sqrt(outML$pGamma))
  gammas_model[9] = mean(1/sqrt(outPog$pGamma))
  
  sampOutMl = sampEstG(datIll[datIll$task==1,])
  sampOutPog = sampEstG(datIll[datIll$task==2,])  
  
  gammas_sample[8] = sampOutMl
  gammas_sample[9] = sampOutPog

  
  l.r7[1] = lSizeForR(r=.7, g2=gammas_sample[1]^2)
  l.r9[1] = lSizeForR(r=.9, g2=gammas_sample[1]^2)
  l.r7[2:7] = lSizeForR(r=.7, g2=gammas_model[c(2:7)]^2, cond = T)
  l.r9[2:7] = lSizeForR(r=.9, g2=gammas_model[c(2:7)]^2, cond = T)
  l.r7[8:9] = lSizeForR(r=.7, g2=gammas_model[8:9]^2)
  l.r9[8:9] = lSizeForR(r=.9, g2=gammas_model[8:9]^2)
  
  tab = cbind(rep(NA,9), 
              gammas_sample, 
              gammas_model, 
              l.r7, 
              l.r9)
  tab = as.data.frame(tab)
  tab[,1] = contrast
  saveRDS(tab, "table.RDS")
  return(tab)
}

publishTab = function(tab){  
  tab[,2:3] = round(tab[,2:3], 2)
  tab[,4:5] = ceiling(tab[,4:5])
  tab[1,3] = "-" 
  library(papaja)
  studies=c("1. Weight",
                  "2. Effect",
                  "3. Speed",
                  "4. Effect",
                  "5. Speed",
                  "6. Effect",
                  "7. Speed",
                  "8. Mueller-Lyar",
                  "9. Poggendorf")
  tab = cbind(studies,tab)
  
  colnames(tab)=c(
    "Tasks",
    "Contrast",
    "Eq",
    "Model",
    "r = 0.7",
    "r = 0.9")
  
  tab.fig = apa_table(tab,
                      escape = F,
                      col_spanners = list(
                        'Design' = 2, 
                        'Signal-To-Noise SD
                        Estimates' = c(3, 4), 
                        'Needed Trial Size' = c(5,6)),
                      stub_indents = list(
                        'Body Measure ' = c(1,1),
                        "Stroop" = c(2,3),
                        "Flanker" = c(4,5),
                        "Lexical Distance " = c(6,7),
                        "Illusions" = c(8,9)),
                      note="The standard deviation of repeated weight measurements was assumed at 3 lbs.",
                      format.args = list(margin=1),
                      digits=c(0,0,0,2,0,0,0,0),
                      align = c("l",rep("c",5))
  )
  return(tab.fig)
}
