

TAGS <- function(donnee,inits){

  sortie.BG <- function(dat){
    BG <- array(c(dat$pinit[,1],as.vector(dat$ppinit)),dim=c(1,dat$pop+2*dat$test))
    BG[1,1:(dat$pop+dat$test)] <- 1-BG[1,1:(dat$pop+dat$test)]
    dimnames(BG) <- (list(c('Best Guess'),c(paste(c('Prevalence in population'),1:dat$pop,sep=''),paste(c('Specificity'),1:dat$test,sep=''),paste(c('Sensitivity'),1:dat$test,sep=''))))
    return(BG)
  }
  
  comb <- function (dat) {
    comb <- array(0,dim=c(dat$nombre,dat$test))
    for (p in 1:(dat$nombre-1)){
      q <- p
      for (t in 1:dat$test) {
      comb[p+1,t] <- q-2*floor(q/2)
      q <- floor(q/2)
      }
    }
  return(comb)
  }
  
  dat <- NULL
  dat$test <- donnee[1]
  dat$pop <- donnee[2]
  dat$nbref <- sum(donnee[3:4])
  dat$nombre <- 2**dat$test
  dat$par <- 2*dat$test+dat$pop
  dat$ddl <- (dat$nbref+dat$pop)*(2**dat$test-1)
  dat$a <- comb(dat)
  deb <- 5+dat$pop*dat$nombre
  dat$compte <- array(donnee[5:(deb-1)],dim=c(dat$nombre,dat$pop))
  dat$ref <- array(donnee[deb:length(donnee) ],dim=c(dat$nombre,2))
  dat$N <- apply(dat$compte,2,sum)
  dat$Nref <- apply(dat$ref,2,sum)
  dat$pinit  <- matrix(c(1-inits[1:dat$pop],inits[1:dat$pop]),ncol=2)
  dat$ppinit <- matrix(c(inits[(dat$pop+1):length(inits)]),ncol=2,byrow=T)
  dat$ppinit[,1] <- 1-dat$ppinit[,1]

  dat$rec <- cbind(dat$a,dat$compte,dat$ref)
  dimnames(dat$rec) <- list(c(1:dat$nombre),c(paste(c('Result test '),1:dat$test,sep=''),
  paste(c('Number in Population'),1:dat$pop,sep=''),'Reference: Disease Free','Reference: Infected'))
#  print(dat$rec)

  
   # FONCTION EVALUATE
  EM <- function(dat) {
    pp <- dat$ppinit
    p <-dat$pinit
    pitmoins <- rep(0,2*dat$pop+2*dat$test)
    for (i in 1:1000){
      y <- Vref(p,pp)
      Vind <- y %*% t(p)
      x <- matrix(y[,1],nrow=dat$nombre,ncol=dat$pop) * matrix(p[,1],ncol=dat$pop,nrow=dat$nombre,byrow=T) / Vind * dat$compte
      p[,1] <- apply(x,2,sum)/dat$N
      p[,2] <- 1 - p[,1]
      pp[,1] <- apply(t(dat$a) %*% (x+dat$ref[,1]),1,sum)/((sx<-sum(x))+dat$Nref[1])
      pp[,2] <- apply(t(dat$a) %*% (dat$compte-x+dat$ref[,2]),1,sum)/(sum(dat$N)-sx+dat$Nref[2])
      pit <- c(as.vector(p),as.vector(pp))
      if (all.equal(pit,pitmoins)==T) break else pitmoins <- pit
    }
    browser()
    if(sum(1-pp[,1],pp[,2]) < sum(pp[,1],1-pp[,2])){pp <- cbind(pp[,2],pp[,1]);p <- cbind(p[,2],p[,1])}
    ptot <- rbind(p,c(1,0),c(0,1))
    ll <- sum(cbind(dat$compte,dat$ref)*log(Vref(ptot,pp) %*% t(ptot)),na.rm=T)
    sol <- array(c(as.vector(p[,1]),as.vector(pp)),dim=c(1,dat$pop+2*dat$test))
    sol[1:(dat$pop+dat$test)] <- 1-sol[1:(dat$pop+dat$test)]
    dimnames(sol) <- list(c('Est'),
    c(paste(c('Prevalence'),1:dat$pop,sep=''),paste(c('Specifity test'),1:dat$test,sep=''),paste(c('Sensitivity test'),1:dat$test,sep='')))
    list(Iterations=i,LogLikelihood=ll,Estimations=round(sol,4),p=p,pp=pp)
  }
  
  vv <- function(i,pp) {t(pp[,i]*t(dat$a))+t((1-pp[,i])*t(1-dat$a))}
    Vref <- function(p,pp){
    z <- lapply(1:2,vv,pp)
    Vref <- lapply(z,apply,1,prod)
    matrix(c(Vref[[1]],Vref[[2]]),ncol=2)
  }
 
  # NEWTON RAPHSON
  ML <- function(dat,p,pp) {
    para <- logit(c(as.vector(p[,1]),as.vector(pp)))
    res.nlm <<- nlm(llik, para, hessian=TRUE, ndigit=10)
    x <- res.nlm$hessian
    sol <- res.nlm$estimate
    if (prod(diag(qr(x)$qr))*(-1)^(ncol(x)-1)==0) {
      remInv <- 'The Matrix is singular : no SE available';
      SE <- rep(NaN,2*dat$test+dat$pop)} 
    else { remInv <- "";
           SE <- sqrt(diag(solve(x)))}
    
    p <- inv.logit(sol[1:dat$pop])
    p <- cbind(p,(1-p))
    pp <- inv.logit(sol[(dat$pop+1):(dat$pop+2*dat$test)])		
    pp <- array(pp,dim=c(dat$test,2))		
    sol <- rbind(sol,sol - 1.96*SE,sol + 1.96*SE)
    sol <- inv.logit(sol)
    sol[,1:(dat$pop+dat$test)] <- 1-sol[,1:(dat$pop+dat$test)]
    sol[3:2,1:(dat$pop+dat$test)] <- sol[2:3,1:(dat$pop+dat$test)]
    dimnames(sol) <- list(c('Est','CIinf','CIsup'),
    c(paste(c('pre'),1:dat$pop,sep=''),paste(c('Sp'),1:dat$test,sep=''),paste(c('Se'),1:dat$test,sep='')))
    list(Iterations=res.nlm$iterations,LogLikelihood=-res.nlm$minimum,Estimations=round(sol,4),remInv=remInv,p=p,pp=pp)
  }
  
  # VRAISSEMBLANCE
  llik<-function(para) {
    p <- inv.logit(matrix(c(para[1:dat$pop],-para[1:dat$pop]),ncol=2,byrow=F))
    p <- rbind(p,c(1,0),c(0,1))
    pp <- inv.logit(array(para[(dat$pop+1):(dat$test*2+dat$pop)],dim=c(dat$test,2)))
    -sum(cbind(dat$compte,dat$ref)*log(Vref(p,pp) %*% t(p)))
  }
  
  logit <- function(p) {log(p/(1-p))}
  inv.logit <- function(x) {exp(x)/(1+exp(x))}
  
  # CHECK FOR CORRELATION
  check <- function(pest,ppest,dat){
  if (dat$par >= dat$ddl){return(list(ResCor=NULL, Commentary='The number of parameters = df: no meaningful residuals'))}
  if (dat$test == 1){return(list(ResCor=NULL, Commentary='The number of test=1: no residuals'))}
  res <- NULL
  for (p in 1:dat$pop){ 
    cor <- lab <- mue <- NULL
    muo <- apply(dat$compte[,p]*dat$a,2,sum)/sum(dat$compte[,p])
    for (i in 1:dat$test) {mue[i] <- pest[p,1]*(ppest[i,1])+pest[p,2]*ppest[i,2]}
    for (i in 1:(dat$test-1)) {
      for (j in (i+1):dat$test) {
      obs <-  sum(dat$a[,i]*dat$a[,j]*dat$compte[,p])/sum(dat$compte[,p]) 
      obs <- (obs-muo[i]*muo[j])/sqrt(muo[i]*(1-muo[i])*muo[j]*(1-muo[j]))
      att <- pest[p,1]*ppest[i,1]*ppest[j,1]+pest[p,2]*ppest[i,2]*ppest[j,2]
      att <- (att-mue[i]*mue[j])/sqrt(mue[i]*(1-mue[i])*mue[j]*(1-mue[j]))
      cor <- c(cor,obs-att)
      lab <- c(lab,paste('Correlation test ',i,' - test ',j,sep=''))
      }
    }
    res <- rbind(res,cor)
  }
  dimnames(res) <- list(paste('Population',1:dat$pop,':'),lab)
  list(ResCor=res, Commentary='The residuals should be randomly distributed around 0')
  }
  
  # LIKELIHOOD RATIO
  LR <- function(p,pp,L.ML,dat){
    p <- rbind(p,c(1,0),c(0,1))
    Vind <- Vref(p,pp) %*% t(p)
    Exp <- t(t(Vind[,1:dat$pop])*dat$N)
    Expref <- t(t(Vind[,(dat$pop+1):(dat$pop+2)])*dat$Nref)
    tot <- cbind(dat$rec,round(Exp,2),round(Expref,2))
    dimnames(tot) <- list(c(1:dat$nombre),c(paste(c('test'),1:dat$test,sep=''),paste(c('pop'),1:dat$pop,sep=''),'RefInd','RefInf',paste(c('ExpPop'),1:dat$pop,sep=''),'ExpRefInd','ExpRefInf'))
    df <- dat$ddl-dat$par
    if (dat$par<dat$ddl) {
      L.Tot <- sum(dat$compte*log(t(t(dat$compte)/dat$N)),na.rm=TRUE)+sum(dat$ref*log(t(t(dat$ref)/dat$Nref)),na.rm=TRUE)
      Dev <- 2*(L.Tot-L.ML)
      LRes <- matrix(c(L.Tot,L.ML,Dev,df,p <- 1-pchisq(Dev,df)),nrow=1)
      dimnames(LRes) <- list(c(''),c('Max LogLikelihood: Achievable','Obtained','Deviance','d.f.','p value'))
      if (p < 0.05) com <-'The model does not fit: Assumptions between test may be not justified'
      
      else com<-'The model could fit'
    }
    
    if (dat$par>=dat$ddl) { 
      tot <- NA
      LRes <- NA
      com <-'The number of parameters = ddl: no goodness-of-fit test possible'}
    
  list(Expected=tot,Test=LRes,Commentary=com)
  }
  
  # MODELISATION
#  cat('\\n EXPECTATION MAXIMISATION\\n')
  resEM <<- EM(dat)
#  print(resEM[1:3])
#  cat('\\n NEWTON-RAPHSON \\n')
  p <- array(pmax(0.0001,pmin(.9999,resEM$p)),dim=dim(resEM$p))
  pp <- array(pmax(0.0001,pmin(.9999,resEM$pp)),dim=dim(resEM$pp))
  resML <- ML(dat,p,pp)
#  print(resML[1:3])
#  cat('\\nExpected Results (NR) and Goodness-of-fit test\\n')
  resLR <- LR(resML$p,resML$pp,resML$LogLikelihood,dat)
#  print(resLR)
#  cat('\\nResiduals  correlations between test\\n')
  resCorr <- check(resML$p,resML$pp,dat)
#  print(resCorr)
#browser()

  sortie <- list(
    Com1 = paste('<h4>TAGS V.2.0 results</h4>',
             '<h5>DATA SUMMARY</h5>',
             dat$pop,'Population(s);',
             dat$test,'Tests; ',
             dat$nbref,'Reference Population(s)<br>',
             'degrees of freedom:',dat$ddl,
             '; parameters:',dat$par,'<br>'),
    BG = sortie.BG(dat),
    rec = dat$rec,
    EM = resEM,
    NR = resML,
    ER = resLR,
    Corr = resCorr,
    rem1 = '<strong>WARNING 1</strong>: test results are assumed to be independent conditionally on infection or disease status',
    rem2 = if ((dat$pop+dat$nbref)>1) '<strong>WARNING 2</strong>: tests are supposed to have constant sensitivity and specificity in all populations' else ""
  )
 

  return(sortie)
} #end TAGS

