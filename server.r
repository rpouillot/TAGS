
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
    if(any(is.na(pp))) return(NA)
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
  if(is.na(resEM)) return(NA)
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

shinyServer(function(input, output,session) {
  
  form1 <- reactive({
    test <- as.numeric(input$qtest)
    pop <- as.numeric(input$qpop)
    para <- test*2 + pop
    if(input$qref == "no reference population") {ref = 0; nbref <- 0}
    if(input$qref == "1 disease free population") {ref = 1; nbref <- 1}
    if(input$qref == "1 infected population") {ref = 2; nbref <- 1}
    if(input$qref =="1 infected and 1 disease free populations") {ref = 3; nbref <- 2}
    ddl <- (pop + nbref) * (2^test - 1)
    return(list(test=test, pop=pop,para=para,ref=ref,nbref=nbref,ddl=ddl,nombre=2^test))
  }) 
  
output$para <- renderText(form1()$para)  
output$ddl  <- renderText(form1()$ddl)  

Combinaison <- reactive({
  test <- form1()$test 
  comb <- comb2 <- NULL
  q <- 0
  for (p in 0:(2^test-1)) {
    comb[p+1] <- comb2[p+1] <- ""
    q <-  p
    for (t in 1:test){
      pipo <- q-2*floor(q/2)
      comb[p+1] <- paste(comb[p+1],"test #", t ," = ",ifelse(pipo,"Positif","Negatif"),sep="")
      if (t!=test) comb[p+1] <- paste(comb[p+1],"; ", sep="")
      comb2[p+1] <- paste(comb2[p+1]," ",pipo, sep="")
      q <- floor(q/2)
    }
  }
  return(list(comb=comb, comb2=comb2))
})

EntreRes <- function(nb, typePop, Pop){
  txt <- ""
  counter <- 1
  for (i in 1:nb){
    label <- paste("Enter the number of samples with a ",Combinaison()$comb[i]," (i.e. result: ",Combinaison()$comb2[i],"): ",sep="")
    txt <- paste(txt,numericInput(inputId=paste(typePop,"-",Pop,"-",counter,sep=""),label,value=0),"</p>")

#    cat(paste(typePop,"-",Pop,"-",counter,sep=""),"\n")    
    
    counter <- counter + 1
  }
  return(txt)
}

output$UI <- renderUI({

  if(input$qtest == 0 | input$qpop == 0) return(HTML("")) 
  if (form1()$para > form1()$ddl){
    return(HTML('The number of parameters > df'))}
  
  formu <- "<hr><font color='#000080'><b>4: Enter your experimental results </b></font><br><br>" 
  for (k in 1:(form1()$pop)){
    formu <- paste(formu,"<b> Population with an unknown infection status #",k,":</b><br>",sep="")
    formu <- paste(formu, EntreRes(form1()$nombre,"UN",k),sep="")
    formu <- paste(formu,"<hr>",sep="")
  }
  
  if (form1()$ref==1 || form1()$ref==3) {
    formu <- paste(formu,"<b> Disease Free Reference Population:</b><br>",sep="")
    formu <- paste(formu, EntreRes(form1()$nombre,"DF",1),sep="")
    formu <- paste(formu,"<hr>",sep="")	
  }
  
  if (form1()$ref==2 || form1()$ref==3) {
    formu <- paste(formu,"<b> Infected Reference Population:</b><br>",sep="")
    formu <- paste(formu,EntreRes(form1()$nombre,"IN",1),sep="")
    formu <- paste(formu,"<hr>",sep="")
  }
  
  formu <- paste(formu,"<font color='#000080'><b>5: (Facultative) Enter your best guesses (all values within ] 0-1 [)</b></font> <br>")
  
  for (i in 1:form1()$pop){
    formu <- paste(formu,"Disease prevalence in the population with unknown infection status #",i,": ")
    formu <- paste(formu,numericInput(inputId=paste("PR",i,sep=""),"",0.5),"<br>")
  }
  formu <- paste(formu,"<br>")
  
  for (i in 1:form1()$test){
    formu <- paste(formu,"test #",i,": Specificity: ",sep="")
    formu <- paste(formu,numericInput(inputId=paste("SP",i,sep=""),"",0.95), "<br>")
    formu <- paste(formu,"test #",i,": Sensitivity: ",sep="")
    formu <- paste(formu,numericInput(inputId=paste("SE",i,sep=""),"",0.8),"<br>")
  }

  return(HTML(formu))
  })

oTAGS <- reactive({
  if(input$button == 0) return(NULL)
  isolate({
  w <- c(form1()$test,
         form1()$pop,
         (form1()$ref==1 || form1()$ref==3)*1,
         (form1()$ref==2 || form1()$ref==3)*1)
  x <- NULL
  
  for(i in 1:form1()$pop){
    for(j in 1:form1()$nombre)
      x <- c(x,input[[paste("UN",i,j,sep="-")]])
  }

  for(j in 1:form1()$nombre){
    add <- if(form1()$ref %in% c(1,3)) input[[paste("DF",1,j,sep="-")]] else 0
    x <- c(x,add)
  }
  
  for(j in 1:form1()$nombre){
    add <- if(form1()$ref %in% c(2,3)) input[[paste("IN",1,j,sep="-")]] else 0
    x <- c(x,add)
  }
  
  y <- NULL
  for(j in 1:form1()$pop){
    y <- c(y, input[[paste("PR",j,sep="")]])
  }
  for(j in 1:form1()$test){
    y <- c(y, input[[paste("SP",j,sep="")]], input[[paste("SE",j,sep="")]])
  }
# cat("x\n",x,"y\n",y,"\n")
# debug(TAGS)
# print(res)

  if(sum(x) == 0) return("NoData")
  res <- TAGS(c(w,x),y) 
  if(is.na(res[1])) return("NA")
  return(res)})
})


output$sortie1 <- renderUI({
  if(input$button == 0) return(NULL)
  isolate({
    if (form1()$para > form1()$ddl){
      return(HTML('The number of parameters > df'))}
    if (oTAGS()[1] == "NoData"){
      return(HTML('No Data'))}
    if (oTAGS()[1] == "NA"){
      return(HTML('No Meaningful Evaluation'))}
    
    HTML(oTAGS()$Com1,
         '<br>',
         print(xtable(oTAGS()$rec,
            display = rep("d",form1()$test+form1()$pop+2+1)
            ),type="html",file="pipo"),
         '<br>', 
         print(xtable(oTAGS()$BG),type="html",file="pipo"),
         
         '<br><h5>EXPECTATION MAXIMISATION</h5>',
         'Iterations: ',
         oTAGS()$EM$Iterations,'<br>',
         'LogLikelihood: ',
         round(oTAGS()$EM$LogLikelihood,2),'<br>',
         'Estimations<br>',
         print(xtable(oTAGS()$EM$Estimations,
                                digits=rep(4,1 + form1()$pop + 2*form1()$test)
                                ),type="html",file="pipo"),
         
         '<br><h5>NEWTON-RAPHSON</h5>',
         'Iterations: ',
         oTAGS()$NR$Iterations,'<br>',
         'LogLikelihood: ',
         round(oTAGS()$NR$LogLikelihood,2),'<br>',
         'Estimates<br>',
         print(xtable(oTAGS()$NR$Estimations,
                                digits=rep(4,1 + form1()$pop+2*form1()$test)
                    ),type="html",file="pipo"),'<br>',
         
         oTAGS()$NR$remInv,
         '<br>',
         oTAGS()$rem1,'<br>',
         oTAGS()$rem2,'<br>',

         '<br><h5>Expected Results (NR) and Goodness-of-fit test</h5>',
         if(!is.na(oTAGS()$ER$Expected[1])) 
           print(xtable(oTAGS()$ER$Expected,
                                     digits = c(0,
                                                rep(0,form1()$test+form1()$pop+2),
                                                rep(2,form1()$pop+2))
                        ),type="html",file="pipo"),
         '<br><h5>Test</h5>',
         if(!is.na(oTAGS()$ER$Test[1])) 
           print(xtable(oTAGS()$ER$Test), type="html",file="pipo") else "",
         '<br>',
         'Commentary<br>',
         oTAGS()$ER$Com,'<br>',
         
         '<br><h5>Residuals  correlations between tests</h5><br>',
         if(!is.null(oTAGS()$Corr$ResCor)) 
           print(xtable(oTAGS()$Corr$ResCor),type="html",file="pipo") else "",'<br>',
         oTAGS()$Corr$Commentary,'<br>'
         
         
         
         
    ) #End HTML
  })
})




})
  