


rm(list=ls())

normal.sim.single <- function(sim,Mus,Sigmas,ns,data.id,trace)
{ 
  
  
  directional.p <- function(k,mu.tilde,ns,datai,plot=FALSE,npts=100,step.eps=10^-6)
  {  # directional p-value
    p <- length(mu.tilde)
    d <- p*(k-1)
    n <- sum(ns)
    
    
    lambdai <- sapply(datai, function(x) t(x$nu) %*% x$Deltahat %*% x$nu)
    
    tsup <- min(sqrt(1+1/lambdai))
    
    
    gg <- function(x)
    {
      
      detxi <-  sapply(datai,function(z){
        determinant(z$S + (1 - x^2) * z$nu %*% t(z$nu), logarithm = TRUE)$modulus 
      })
      C1i <- lapply(datai, function(z) z$n * z$Deltatilde %*% (  (1- x^2 * sum(diag(z$Deltatilde %*% z$A0  ))) * z$S.tilde - x^2 * z$A0 ) %*% z$Deltatilde     )
      C1 <- Reduce("+",C1i)
      detc1x <- as.numeric(determinant(C1,logarithm=TRUE)$modulus)
      
      gbar <- ( (d-1) * log(x) 
                + sum( 0.5 * (ns- p -2) *  detxi)
                + 0.5*detc1x)
      gbar
    }
    
    
    
    information <- function(x){
      
      .C1i <- lapply(datai, function(z) z$n * z$Deltatilde %*% (  (1- x^2 * sum(diag(z$Deltatilde %*% z$A0  ))) * z$S.tilde - x^2 * z$A0 ) %*% z$Deltatilde     )
      .C1 <- Reduce("+",.C1i)
      .invC1 <- solve(.C1)
      
      .deriv.C1i <- lapply(datai, function(z) z$n  * z$Deltatilde %*% ( sum(diag(z$Deltatilde %*% z$A0 )) * z$S.tilde + z$A0  ) %*% z$Deltatilde )
      .deriv.C1 <- Reduce("+",.deriv.C1i)
      
      tri <- sapply(datai, function(z){
        .Deltaxi <- solve( (z$S + (1 - x^2) * z$nu %*% t(z$nu)) )
        sum(diag( 4*x^2* .Deltaxi %*% z$A0 %*% .Deltaxi %*% z$A0 + 2 * .Deltaxi %*% z$A0 ) )
      })
      trc1 <- sum(diag( 4*x^2*.invC1 %*% .deriv.C1 %*% .invC1 %*% .deriv.C1 + 2* .invC1 %*% .deriv.C1  ))
      
      ( (d-1) / (x^2) + sum(0.5 * (ns-p-2) * tri) + 0.5*trc1)
    }
    
    tmin <- function(w) {
      if((t_hat-w*sig) < 0) 0
      else if((t_hat-w*sig) < 1 & (t_hat-w*sig) > 0  ) t_hat-w*sig
      else tmin(w+2)
    }
    
    tmax <- function(w) {
      if((t_hat+w*sig) >tsup ) tsup
      else if((t_hat+w*sig) > 1 & (t_hat+w*sig) < tsup  ) t_hat+w*sig
      else tmax(w+2)
    }
    
    
    t_hat <- optimize(gg, c(0, tsup),maximum=TRUE,tol = 500*.Machine$double.eps)$maximum
    gbarhat <- gg(t_hat)
    sig <- 1/sqrt(information(t_hat))
    t_min <- tmin(5)
    t_max <- tmax(5)
    
    
    ff <- function(x)
    {
      
      detxi <-  sapply(datai,function(z){
        determinant(z$S + (1 - x^2) * z$nu %*% t(z$nu), logarithm = TRUE)$modulus 
      })
      C1i <- lapply(datai, function(z) z$n * z$Deltatilde %*% (  (1- x^2 * sum(diag(z$Deltatilde %*% z$A0  ))) * z$S.tilde - x^2 * z$A0 ) %*% z$Deltatilde     )
      C1 <- Reduce("+",C1i)
      detc1x <- as.numeric(determinant(C1,logarithm=TRUE)$modulus)
      
      
      val <- exp ( ( (d-1) * log(x) 
                     + sum( 0.5 * (ns- p -2) *  detxi)
                     + 0.5*detc1x) - gbarhat) 
      val
    }
    
    ff.v <- Vectorize(ff)
    if (plot) plot(ff.v,t_min,t_max,n=npts)
    
    
    up <- integrate(ff.v,lower=1,upper=t_max)#;up
    down <- integrate(ff.v,lower=t_min,upper=1)#;down
    up$value/(down$value+up$value)
    
  }
  
  W.skovgaard <- function(k,mu.tilde,ns,datai)
  {  # W*
    p <- length(mu.tilde)
    d <- p * (k-1)
    
    C1i <- lapply(datai, function(x){
      (x$n * x$Deltatilde %*% (  (1-sum(diag(x$Deltatilde %*% x$A0  ))) * x$S.tilde - x$A0 ) %*% x$Deltatilde   )
    } )
    C1 <- Reduce("+",C1i)
    ilambdalambdai <- lapply(datai, function(x) x$n * x$Deltatilde)
    ilambdalambda <- Reduce("+",ilambdalambdai)
    
    
    
    ldethati <- sapply(datai, function(x) as.numeric(determinant(x$Deltahat,logarithm=TRUE)$modulus) )
    ldettildei <- sapply(datai, function(x)as.numeric(determinant(x$Deltatilde,logarithm=TRUE)$modulus) )
    ldetC1 <- as.numeric(determinant(C1,logarithm=TRUE)$modulus)
    ldeti <- as.numeric(determinant(ilambdalambda,logarithm=TRUE)$modulus)
    
    w <- sum(ns*ldethati) - sum(ns*ldettildei)
    
    
    
    ## numerator
    # ntri <- sapply(datai,function(x) x$n * t(x$nu) %*% x$Deltatilde %*% x$nu)
    # ntr <- Reduce("+",ntri)
    num <- as.numeric( (d / 2) * log( sum(sapply(datai,function(x) x$n * t(x$nu) %*% x$Deltatilde %*% x$nu))) )
    
    ## denominator 
    # dtri <- sapply(datai, function(x) x$n * t(x$nu) %*% x$Deltahat %*% x$nu)
    # dtr <- Reduce("+",dtri)
    den <- as.numeric( log( sum( sapply(datai, function(x) x$n * t(x$nu) %*% x$Deltahat %*% x$nu) ) ) )
    
    ## information matrix 
    fact <- ( (p + 2) / 2 ) * sum(ldethati - ldettildei) 
    
    nuisance <- (0.5*(ldetC1 - ldeti))
    
    
    #gamma
    loggamma <- ( num - den - ( (d / 2) - 1) * log(w) + fact +nuisance )
    
    wstar1 <- w * (1 - loggamma / w) ^ 2
    wstar2 <- w - 2 * loggamma
    
    p1 <- pchisq(wstar1,d,lower=FALSE)
    p2 <- pchisq(wstar2,d,lower=FALSE)
    
    list(wstar1=wstar1,p.value1=p1,wstar2=wstar2,p.value2=p2,loggamma=loggamma)
  }
  
  
  
  lrt <- function(k,Si,Si.tilde,ns)
  { # first-order LRT from marginal likelihood based on S
    p <- sapply(Si,ncol)[1]
    k <- length(ns)
    d <- p * (k-1)
    
    ldeti.hat <- sapply(Si, function(x) as.numeric(determinant(x,logarithm=TRUE)$modulus) )
    ldeti.tilde <- sapply(Si.tilde, function(x)as.numeric(determinant(x,logarithm=TRUE)$modulus) )
    
    W <- sum(ns*ldeti.tilde) - sum(ns*ldeti.hat) 
    
    pvalue <- pchisq(W,d,lower=FALSE)
    
    list(Woss=W,pvalue=pvalue)
  }
  
  
  
  BehrensFisher <- function(ybari,Si.unbias,ns)
  {
    S1.unbias <- Si.unbias[[1]]
    S2.unbias <- Si.unbias[[2]]
    S.unbias <- S1.unbias / ns[1] + S2.unbias / ns[2]
    p <- ncol(S1.unbias)
    y1bar <- ybari[[1]]
    y2bar <- ybari[[2]]
    n1 <- ns[1]
    n2 <- ns[2]
    T2star <- drop(  t(y1bar - y2bar) %*% solve(S.unbias) %*% (y1bar - y2bar) )
    df2.u <-  sum(diag( S.unbias %*% S.unbias )) + ( sum(diag( S.unbias )) )^2
    df2.d1 <- (sum(diag(1/n1^2 * S1.unbias %*% S1.unbias )) + (sum(diag(S1.unbias/n1)))^2)/(n1-1)
    df2.d2 <- (sum(diag(1/n2^2 * S2.unbias %*% S2.unbias )) + (sum(diag(S2.unbias/n2)))^2)/(n2-1)
    
    df2 <- df2.u / (df2.d1 + df2.d2) - p +1
    df1 <- p
    BF <- df2 * T2star / df1 / (df2+p-1) 
    
    if(BF > 0){
      BF.pvalue <- pf(BF, df1, df2 ,lower=FALSE)
    }else{
      BF.pvalue <- pf(BF, df1, df2,lower=TRUE)
    }
    
    list(BF=BF,BF.pvalue=BF.pvalue,df1=df1,df2=df2)
  }
  
  
  ## negative profile log-likelihood
  nprofile <- function(mu,data)
  {
    sum(sapply(data, function(x){
      x$n /2 * log(1+t(x$ybar -mu) %*% x$Delta %*% (x$ybar - mu) ) 
    }))
  }
  
  if (trace) print(data.id[[sim]]$id)
  
  y <- as.list(numeric(length(ns)))
  for (kk in 1:length(ns)) {
    y[[kk]] <- rmvnorm(ns[kk],mean=Mus[kk,],sigma=Sigmas[[kk]])
  }
  
  ###############
  p <- sapply(y,ncol)[1]
  ns <- sapply(y, nrow)
  n <- sum(ns)
  k <- length(ns)
  df <- p*(k-1)
  
  ## MLE
  ybari <- lapply(y, colMeans)
  # names(ybari) <- 1:k
  
  Si.unbias <- lapply(y, cov)
  Si <- lapply(y, function(x) (nrow(x)-1)/(nrow(x))*cov(x) )
  Deltai <- lapply(Si,solve)
  # names(Si) <- 1:k
  
  
  datai <- list()
  for (i in 1:k) {
    datai[[i]] <- list(ybar=ybari[[i]],S=Si[[i]],Delta=Deltai[[i]],n = ns[i])
  }
  
  
  mu.result <-  optim(par=runif(p), fn=nprofile,
                      method = 'BFGS',
                      data =datai)
  
  mu.tilde <- mu.result$par
  #####################
  ### here 
  ####################
  Si.tilde <- lapply(datai, function(x){
    x$S + (x$ybar - mu.tilde) %*% t(x$ybar - mu.tilde)
  })
  
  
  Deltatildei <- lapply(Si.tilde,solve)
  nui <- lapply(ybari,function(x) x-mu.tilde)
  Ai0 <- lapply(nui,function(x) x %*% t(x))
  
  datai <- list()
  for (i in 1:k) {
    datai[[i]] <- list(ybar=ybari[[i]],S=Si[[i]],S.tilde=Si.tilde[[i]],
                       Deltahat=Deltai[[i]],Deltatilde=Deltatildei[[i]],
                       nu=nui[[i]],A0 = Ai0[[i]],n = ns[i])
  }
  
  
  # first-order likelihood ratio statistic (W) 
  firstorder <- lrt(k=k,Si=Si,Si.tilde=Si.tilde,ns=ns)
  W <- firstorder$Woss
  W_P <- firstorder$pvalue
  
  
  # directional p-values
  ########################
  
  
  D_P <- directional.p(k=k,mu.tilde=mu.tilde,ns=ns,datai=datai,plot = FALSE)
  W_D <- qchisq(D_P, df,lower.tail = FALSE)
  # Skovgaard's p-values
  ########################
  
  app <- W.skovgaard(k=k,mu.tilde=mu.tilde,ns=ns,datai=datai)
  W_star <- app$wstar1
  W_star_P <- app$p.value1
  W_star2 <- app$wstar2
  W_star2_P <- app$p.value2
  loggamma <- app$loggamma
  
  if(k==2){
    behrensfisher <- BehrensFisher(ybari=ybari,Si.unbias=Si.unbias,ns=ns)
    BF <- behrensfisher$BF
    BF_P <- behrensfisher$BF.pvalue
    df1 <- behrensfisher$df1
    df2 <- behrensfisher$df2
    data.frame(W = W, 
               W_star = W_star, 
               W_star2 = W_star2,
               W_D = W_D,
               BF = BF,
               W_P = W_P, 
               W_star_P = W_star_P, 
               W_star2_P = W_star2_P,
               D_P = D_P,
               BF_P = BF_P,
               df1=df1,df2=df2,
               loggamma = loggamma)
  }else{
    data.frame(W = W, 
               W_star = W_star, 
               W_star2 = W_star2,
               W_D = W_D,
               W_P = W_P, 
               W_star_P = W_star_P, 
               W_star2_P = W_star2_P,
               D_P = D_P,
               loggamma = loggamma)
  }
  
  
}


############################################
####' linear relationship between p and n 
####' 
#############################################

library(parallel)
library(mvtnorm)
k <- 2
Nsim <- 10000
ra <- 0.3
seed <- 123
num <- c(300,600,900,1200,1500)

for (iii in 1:length(num)) {
  set.seed(seed)
  RNGkind("L'Ecuyer-CMRG")
  ns <- rep(num[iii],k)
  
  c <- ra
  p <- c*ns[1]
  p.interest <- p*(k-1)
  
  ### same covariance matrices
  Sigmas <- as.list(numeric(k))
  sigma <- as.matrix(diag(rep(1,p)))
  for (ii in 1:k) {
    Sigmas[[ii]] <- sigma
  }
  
  
  ### different covariance matrices AR model sparse
  # Sigma1 <- matrix(rep(0,p*p),nrow=p)
  # Sigmas <- as.list(numeric(k))
  # rho <- seq(0.1,0.9,length.out=k)
  # for (ii in 1:k) {
  #   Sigmas[[ii]] <- rho[ii]^abs(col(Sigma1)-row(Sigma1))
  # }
  # 
  
  
  ### different covariance matrices dense
  Sigma1 <- matrix(rep(0,p*p),nrow=p)
  Sigmas <- as.list(numeric(k))
  rho <- seq(0.1,0.9,length.out=k)
  for (ii in 1:k) {
    Sigmas[[ii]] <- (1-rho[ii]) * as.matrix(diag(rep(1,p))) + rho[ii] * (rep(1,p)%*% t(rep(1,p)))
  }
  
  
  group <- factor(rep(1:k,each=p))
  mu <- 0
  psi <- rep(0,(k-1)) #group number -1
  X <- model.matrix(~group)
  mean <- as.vector(X%*%c(mu,psi))
  Mus <- matrix(mean,ncol=p,byrow = T )
  
  data.id <- as.list(numeric(Nsim))
  for (jj in 1:Nsim){
    data.id[[jj]] <- list(id=jj)
  }
  
  res <-  mclapply(1:Nsim, FUN=normal.sim.single,
                   Mus=Mus, Sigmas=Sigmas, ns=ns,
                   data.id=data.id,trace=TRUE,
                   mc.cores = 10,mc.set.seed = TRUE)
  
  if(k==2){
    normal.result <- list(W = NULL, W_star = NULL, W_star2 = NULL,
                          W_D = NULL,BF = NULL,
                          W_P = NULL, W_star_P = NULL, W_star2_P = NULL,
                          D_P = NULL,BF_P = NULL,df1=NULL,df2=NULL,loggamma=NULL)
  }else{
    normal.result <- list(W = NULL, W_star = NULL, W_star2 = NULL,
                          W_D = NULL,
                          W_P = NULL, W_star_P = NULL, W_star2_P = NULL,
                          D_P = NULL,loggamma=NULL)
  }
  
  
  for (i in 1:length(res)){
    for (j in 1:length(normal.result)) {
      normal.result[[j]] <- cbind(normal.result[[j]],res[[i]][[j]])
    }
  }
  
  (filename=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_alpha",".RData",sep =""))
  save(normal.result=normal.result,file=filename)
  
}


