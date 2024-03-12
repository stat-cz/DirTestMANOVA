

rm(list=ls())
setwd("C:/Users/HCZ/Desktop/dirmean/example")
normal.example <- function(y)
{ 
  
  
  check.t <- function(S,S0,A,ns,t.start=1,step=0.001,trace=FALSE)
  { # check how far we could go with t in S(t)
    require(matrixcalc)
    n <- sum(ns)
    t.old <- t.new <- t.start
    control=TRUE
    while (control)
    {
      t.old <- t.new
      t.new <- t.old+step
      if (trace) print(t.old)
      St <- S0 - t.new^2 * A/n
      control <- is.positive.definite(St)
    }
    t.old
  }
  
  
  directional.p <- function(S,S0,A,ns,plot=FALSE,npts=100,step.eps=10^-6)
  {  # directional p-value
    p <- nrow(S0)
    k <- length(ns)
    d <- p*(k-1) 
    n <- sum(ns)
    
    
    # .Delta <- solve(S0)
    
    B0.inv <- solve(chol(S0))
    
    Q <- 1/n * t(B0.inv) %*% A %*% B0.inv
    lambda <- eigen(Q,only.values = TRUE)$values
    
    if( class(lambda)=="complex"||any( lambda < 0)   ){
      
      ldetS0 <- as.numeric(determinant(S0,logarithm=TRUE)$modulus)
      
      
      step <- tsup <- 1
      while (step>=step.eps)
      {
        tsup <- check.t(S,S0,A,ns,step=step,t.start=tsup)
        step=step/10
      }
      
      
      gg <- function(x)#ng=log(t^(d-1)h(t))
      {
        Sx <- S0 - x^2 * A/n
        detx <- as.numeric(determinant(Sx,logarithm=TRUE)$modulus)
        gbar <- ((d - 1) * log(x) + 0.5 * (n - p - 1 - k) * (detx - ldetS0) ) 
        gbar
      }
      
      
      
      information <- function(x){
        .Sx <- S0 - x^2 * A/n
        .Deltax <- solve(.Sx)
        
        tr <- sum(diag( 4*x^2/n^2 * .Deltax %*% A %*% .Deltax %*% A + 2/n  * .Deltax %*% A) )
        ( (d-1) / (x^2) + 0.5 * (n-p-1-k) * tr )
      }
      
      ff <- function(x)
      {
        Sx <-  S0 - x^2 * A/n
        detx <- as.numeric(determinant(Sx,logarithm=TRUE)$modulus)
        val <- exp((d - 1) * log(x) + 0.5 * (n - p - 1 - k) * (detx-ldetS0)  - gbarhat) 
        val
      }
      
    }else{
      
      tsup <- sqrt(1/max(lambda))
      
      gg <- function(x)#ng=log(t^(d-1)h(t))
      {
        ((d - 1) * log(x)
         + 0.5 * (n - p - 1-k) * sum(log(1 - x^2 * lambda)) )
        
      }
      
      information <- function(x){
        ((d - 1) / (x ^ 2)
         + 0.5 * (n - p-1-k)* sum((2*lambda + 2* x^2 *lambda^2) / (1 - x^2*lambda )^2) )
      }
      
      ff <- function(x)
      {
        # Sx <- S0 - x^2 * A/n
        # ldetx <- as.numeric(determinant(Sx,logarithm=TRUE)$modulus)
        val <- exp((d - 1) * log(x)
                   + 0.5 * (n - p - 1 -k) * sum(log(1 - x^2 * lambda))
                   - gbarhat)
        val
      }
      
      
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
    
    
    ff.v <- Vectorize(ff)
    if (plot) plot(ff.v,t_min,t_max,n=npts)
    
    # up <- integrate(ff.v,lower=1,upper=tsup,rel.tol = 500*.Machine$double.eps)
    # down <- integrate(ff.v,lower=0,upper=1,rel.tol = 500*.Machine$double.eps)
    # up$value/(down$value+up$value)
    
    
    up <- integrate(ff.v,lower=1,upper=t_max)#;up
    down <- integrate(ff.v,lower=t_min,upper=1)#;down
    up$value/(down$value+up$value)
    
    
  }
  
  W.skovgaard <- function(datai,A,S,S0,ns)
  { # W*
    # skovgaard method need to check. 
    # maybe some problems come from the inverse of observed infromation matrix.
    p <- nrow(S0)
    k <- length(ns)
    d <- p*(k-1)
    n <- sum(ns)
    Deltahat <- solve(S)
    Deltatilde <- solve(S0)
    
    
    ldethat <- as.vector(determinant(Deltahat,logarithm=TRUE)$modulus)
    ldettilde <- as.vector(determinant(Deltatilde,logarithm=TRUE)$modulus)
    
    w <- n*ldethat-n*ldettilde
    
    #gamma
    ybar <- rowSums(sapply(datai, function(x) x$n * x$ybar))/n
    
    num1 <- sum ( sapply(datai, function(x) x$ni * (t(x$ybar -ybar) %*% Deltatilde %*% (x$ybar -ybar) )  ) )
    .A0 <-  ( ( drop(t(ybar)%*%Deltatilde%*%ybar) * Deltatilde + Deltatilde%*%ybar%*%t(ybar)%*%Deltatilde) )
    
    num2 <- ( Reduce("+",lapply(datai,function(x) x$ni * t(x$ybar - ybar) ) ) 
              %*% (.A0/n) %*%  
                Reduce("+",lapply(datai,function(x) x$ni * (x$ybar - ybar) ) ) )
    
    ntr <- num1 + num2    
    num <- (d/2) * log(ntr)  
    
    
    dtr <- sum (diag (Deltatilde %*% A)) 
    den <- log (dtr)
    
    fact <- ((p+1+k)/2)*(ldethat-ldettilde)
    
    loggamma <- (num-den-((d/2)-1)*log(w)+fact)
    
    wstar1 <- w*(1-loggamma/w)^2
    wstar2 <- w-2*loggamma
    
    p1 <- pchisq(wstar1,d,lower=FALSE)
    p2 <- pchisq(wstar2,d,lower=FALSE)
    
    list(wstar1=wstar1,p.value1=p1,wstar2=wstar2,p.value2=p2,loggamma=loggamma)
  }
  
  
  lrt <- function(S,S0,ns)
  { # first-order LRT from marginal likelihood based on S
    p <- nrow(S)
    # d <- p*(2-1)
    k <- length(ns)
    d <- p * (k-1)
    n <- sum(ns)
    
    det0 <- as.numeric(determinant(S0,logarithm=TRUE)$modulus)
    det <- as.numeric(determinant(S,logarithm=TRUE)$modulus)
    
    
    W <- n*det0-n*det
    
    pvalue <- pchisq(W,d,lower=FALSE)
    
    list(Woss=W,pvalue=pvalue)
  }
  
  
  HotellingT2 <- function(ybari,B,ns)
  {
    S.unbias <- 1/(sum(ns)-2) * B
    n1 <- ns[1]
    n2 <- ns[2]
    y1bar <- ybari[[1]]
    y2bar <- ybari[[2]]
    p <- ncol(B)
    T2 <- drop( (n1 * n2) / (n1 + n2) * t(y1bar - y2bar) %*% solve(S.unbias) %*% (y1bar - y2bar) )
    HT <- (n1+n2-2-p+1) * T2 / p / (n1+n2-2) 
    
    if(HT > 0){
      HT.pvalue <- pf(HT, p, (n1+n2-2-p+1),lower=FALSE)
    }else{
      HT.pvalue <- pf(HT, p, (n1+n2-2-p+1),lower=TRUE)
    }
    
    list(HT=HT,HT.pvalue=HT.pvalue)
  }
  

  
  #mle
  # first order
  ###############
  p <- sapply(y,ncol)[1]
  ns <- sapply(y, nrow)
  n <- sum(ns)
  k <- length(ns)
  df <- p*(k-1)
  
  ## MLE
  ybari <- lapply(y, colMeans)
  Si <- lapply(y, function(x) (nrow(x)-1)/(nrow(x))*cov(x) )
  
  datai <- list()
  for (i in 1:k) {
    datai[[i]] <- list(ybar=ybari[[i]],S=Si[[i]],ni = ns[i])
  }
  
  ybar <- rowSums(sapply(datai, function(x) x$ni * x$ybar))/n
  B <- Reduce("+",lapply(datai,function(x) x$ni * x$S ) )
  A <- Reduce("+",lapply(datai, function(x) x$ni * (x$ybar-ybar) %*% t(x$ybar - ybar)  ) )
  
  S <- 1/n * B
  S0 <- 1/n*(A+B)
  
  
  
  # first-order likelihood ratio statistic (W) 
  firstorder <- lrt(S=S,S0=S0,ns)
  W <- firstorder$Woss
  W_P <- firstorder$pvalue
  
  
  # directional p-values
  ########################
  
  D_P <- directional.p(S,S0,A,ns,plot = FALSE)
  W_D <- qchisq(D_P, df,lower.tail = FALSE)
  # Skovgaard's p-values
  ########################
  
  app <- W.skovgaard(datai,A,S,S0,ns)
  W_star <- app$wstar1
  W_star_P <- app$p.value1
  W_star2 <- app$wstar2
  W_star2_P <- app$p.value2
  loggamma <- app$loggamma
  
  ##bartlett correction
  rho <- 1 - 1/n * (1/2*p + 1/2*k +1) 
  
  W_BC <- rho*W
  W_BC_P <- pchisq(W_BC,df,lower.tail = FALSE)
  
  
  if(k==2){
    hotelling <- HotellingT2(ybari,B,ns)
    HT <- hotelling$HT
    HT_P <- hotelling$HT.pvalue
    
    data.frame(W = W, 
               W_star = W_star, 
               W_star2 = W_star2,
               W_BC = W_BC,
               W_D = W_D,
               HT = HT,
               W_P = W_P, 
               W_star_P = W_star_P, 
               W_star2_P = W_star2_P,
               W_BC_P = W_BC_P,
               D_P = D_P,
               HT_P = HT_P,
               loggamma = loggamma)
  }else{
    ## He et al. 2020; central limit theorem
    mun <- n/2 * ( (n-p-k-1/2) * log( (n-1-p) * (n-k) / (n-p-k) / (n-1) ) + (k-1) * log((n-1-p)/(n-1)) + p* log( (n-k)/(n-1) ) ) 
    sigma_square <- 1/2 * (log(1- p/(n-1)) - log(1- p/(n-k)) )
    HE <- ( -W / 2 - mun ) / n / sqrt(sigma_square)
    HE_P <- pnorm(HE,lower.tail = FALSE)
    
    data.frame(W = W, 
               W_star = W_star, 
               W_star2 = W_star2,
               W_BC = W_BC,
               W_D = W_D,
               W_P = W_P, 
               W_star_P = W_star_P, 
               W_star2_P = W_star2_P,
               W_BC_P = W_BC_P,
               D_P = D_P,
               HE = HE,
               HE_P = HE_P,
               loggamma = loggamma)
  }
  
}

normal.sim.single <- function(sim,Mus,Sigmas,ns,data.id,trace)
{ 
  
  
  check.t <- function(S,S0,A,ns,t.start=1,step=0.001,trace=FALSE)
  { # check how far we could go with t in S(t)
    require(matrixcalc)
    n <- sum(ns)
    t.old <- t.new <- t.start
    control=TRUE
    while (control)
    {
      t.old <- t.new
      t.new <- t.old+step
      if (trace) print(t.old)
      St <- S0 - t.new^2 * A/n
      control <- is.positive.definite(St)
    }
    t.old
  }
  
  
  directional.p <- function(S,S0,A,ns,plot=FALSE,npts=100,step.eps=10^-6)
  {  # directional p-value
    p <- nrow(S0)
    k <- length(ns)
    d <- p*(k-1) 
    n <- sum(ns)
    
    
    # .Delta <- solve(S0)
    
    B0.inv <- solve(chol(S0))
    
    Q <- 1/n * t(B0.inv) %*% A %*% B0.inv
    lambda <- eigen(Q,only.values = TRUE)$values
    
    if( class(lambda)=="complex"||any( lambda < 0) ){
      
      ldetS0 <- as.numeric(determinant(S0,logarithm=TRUE)$modulus)
      
      
      step <- tsup <- 1
      while (step>=step.eps)
      {
        tsup <- check.t(S,S0,A,ns,step=step,t.start=tsup)
        step=step/10
      }
      
      
      gg <- function(x)#ng=log(t^(d-1)h(t))
      {
        Sx <- S0 - x^2 * A/n
        detx <- as.numeric(determinant(Sx,logarithm=TRUE)$modulus)
        gbar <- ((d - 1) * log(x) + 0.5 * (n - p - 1 - k) * (detx - ldetS0) ) 
        gbar
      }
      
      
      
      information <- function(x){
        .Sx <- S0 - x^2 * A/n
        .Deltax <- solve(.Sx)
        
        tr <- sum(diag( 4*x^2/n^2 * .Deltax %*% A %*% .Deltax %*% A + 2/n  * .Deltax %*% A) )
        ( (d-1) / (x^2) + 0.5 * (n-p-1-k) * tr )
      }
      
      ff <- function(x)
      {
        Sx <-  S0 - x^2 * A/n
        detx <- as.numeric(determinant(Sx,logarithm=TRUE)$modulus)
        val <- exp((d - 1) * log(x) + 0.5 * (n - p - 1 - k) * (detx-ldetS0)  - gbarhat) 
        val
      }
      
    }else{
      
      tsup <- sqrt(1/max(lambda))
      
      gg <- function(x)#ng=log(t^(d-1)h(t))
      {
        ((d - 1) * log(x)
         + 0.5 * (n - p - 1-k) * sum(log(1 - x^2 * lambda)) )
        
      }
      
      information <- function(x){
        ((d - 1) / (x ^ 2)
         + 0.5 * (n - p-1-k)* sum((2*lambda + 2* x^2 *lambda^2) / (1 - x^2*lambda )^2) )
      }
      
      ff <- function(x)
      {
        # Sx <- S0 - x^2 * A/n
        # ldetx <- as.numeric(determinant(Sx,logarithm=TRUE)$modulus)
        val <- exp((d - 1) * log(x)
                   + 0.5 * (n - p - 1 -k) * sum(log(1 - x^2 * lambda))
                   - gbarhat)
        val
      }
      
      
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
    
    
    ff.v <- Vectorize(ff)
    if (plot) plot(ff.v,t_min,t_max,n=npts)
    
    # up <- integrate(ff.v,lower=1,upper=tsup,rel.tol = 500*.Machine$double.eps)
    # down <- integrate(ff.v,lower=0,upper=1,rel.tol = 500*.Machine$double.eps)
    # up$value/(down$value+up$value)
    
    
    up <- integrate(ff.v,lower=1,upper=t_max)#;up
    down <- integrate(ff.v,lower=t_min,upper=1)#;down
    up$value/(down$value+up$value)
    
    
  }
  
  W.skovgaard <- function(datai,A,S,S0,ns)
  { # W*
    # skovgaard method need to check. 
    # maybe some problems come from the inverse of observed infromation matrix.
    p <- nrow(S0)
    k <- length(ns)
    d <- p*(k-1)
    n <- sum(ns)
    Deltahat <- solve(S)
    Deltatilde <- solve(S0)
    
    
    ldethat <- as.vector(determinant(Deltahat,logarithm=TRUE)$modulus)
    ldettilde <- as.vector(determinant(Deltatilde,logarithm=TRUE)$modulus)
    
    w <- n*ldethat-n*ldettilde
    
    #gamma
    ybar <- rowSums(sapply(datai, function(x) x$n * x$ybar))/n
    
    num1 <- sum ( sapply(datai, function(x) x$ni * (t(x$ybar -ybar) %*% Deltatilde %*% (x$ybar -ybar) )  ) )
    .A0 <-  ( ( drop(t(ybar)%*%Deltatilde%*%ybar) * Deltatilde + Deltatilde%*%ybar%*%t(ybar)%*%Deltatilde) )
    
    num2 <- ( Reduce("+",lapply(datai,function(x) x$ni * t(x$ybar - ybar) ) ) 
              %*% (.A0/n) %*%  
                Reduce("+",lapply(datai,function(x) x$ni * (x$ybar - ybar) ) ) )
    
    ntr <- num1 + num2    
    num <- (d/2) * log(ntr)  
    
    
    dtr <- sum (diag (Deltatilde %*% A)) 
    den <- log (dtr)
    
    fact <- ((p+1+k)/2)*(ldethat-ldettilde)
    
    loggamma <- (num-den-((d/2)-1)*log(w)+fact)
    
    wstar1 <- w*(1-loggamma/w)^2
    wstar2 <- w-2*loggamma
    
    p1 <- pchisq(wstar1,d,lower=FALSE)
    p2 <- pchisq(wstar2,d,lower=FALSE)
    
    list(wstar1=wstar1,p.value1=p1,wstar2=wstar2,p.value2=p2,loggamma=loggamma)
  }
  
  
  lrt <- function(S,S0,ns)
  { # first-order LRT from marginal likelihood based on S
    p <- nrow(S)
    # d <- p*(2-1)
    k <- length(ns)
    d <- p * (k-1)
    n <- sum(ns)
    
    det0 <- as.numeric(determinant(S0,logarithm=TRUE)$modulus)
    det <- as.numeric(determinant(S,logarithm=TRUE)$modulus)
    
    
    W <- n*det0-n*det
    
    pvalue <- pchisq(W,d,lower=FALSE)
    
    list(Woss=W,pvalue=pvalue)
  }
  
  
  HotellingT2 <- function(ybari,B,ns)
  {
    S.unbias <- 1/(sum(ns)-2) * B
    n1 <- ns[1]
    n2 <- ns[2]
    y1bar <- ybari[[1]]
    y2bar <- ybari[[2]]
    p <- ncol(B)
    T2 <- drop( (n1 * n2) / (n1 + n2) * t(y1bar - y2bar) %*% solve(S.unbias) %*% (y1bar - y2bar) )
    HT <- (n1+n2-2-p+1) * T2 / p / (n1+n2-2) 
    
    if(HT > 0){
      HT.pvalue <- pf(HT, p, (n1+n2-2-p+1),lower=FALSE)
    }else{
      HT.pvalue <- pf(HT, p, (n1+n2-2-p+1),lower=TRUE)
    }
    
    list(HT=HT,HT.pvalue=HT.pvalue)
  }
  
  
  if (trace) print(data.id[[sim]]$id)
  
  y <- as.list(numeric(length(ns)))
  for (kk in 1:length(ns)) {
    y[[kk]] <- rmvnorm(ns[kk],mean=Mus[kk,],sigma=Sigmas[[kk]])
  }
  
  
  #mle
  # first order
  ###############
  p <- sapply(y,ncol)[1]
  ns <- sapply(y, nrow)
  n <- sum(ns)
  k <- length(ns)
  df <- p*(k-1)
  
  ## MLE
  ybari <- lapply(y, colMeans)
  Si <- lapply(y, function(x) (nrow(x)-1)/(nrow(x))*cov(x) )
  
  datai <- list()
  for (i in 1:k) {
    datai[[i]] <- list(ybar=ybari[[i]],S=Si[[i]],ni = ns[i])
  }
  
  ybar <- rowSums(sapply(datai, function(x) x$ni * x$ybar))/n
  B <- Reduce("+",lapply(datai,function(x) x$ni * x$S ) )
  A <- Reduce("+",lapply(datai, function(x) x$ni * (x$ybar-ybar) %*% t(x$ybar - ybar)  ) )
  
  S <- 1/n * B
  S0 <- 1/n*(A+B)
  
  
  
  # first-order likelihood ratio statistic (W) 
  firstorder <- lrt(S=S,S0=S0,ns)
  W <- firstorder$Woss
  W_P <- firstorder$pvalue
  
  
  # directional p-values
  ########################
  
  D_P <- directional.p(S,S0,A,ns,plot = FALSE)
  W_D <- qchisq(D_P, df,lower.tail = FALSE)
  # Skovgaard's p-values
  ########################
  
  app <- W.skovgaard(datai,A,S,S0,ns)
  W_star <- app$wstar1
  W_star_P <- app$p.value1
  W_star2 <- app$wstar2
  W_star2_P <- app$p.value2
  loggamma <- app$loggamma
  
  ##bartlett correction
  rho <- 1 - 1/n * (1/2*p + 1/2*k +1) 
  
  W_BC <- rho*W
  W_BC_P <- pchisq(W_BC,df,lower.tail = FALSE)
  
  
  if(k==2){
    hotelling <- HotellingT2(ybari,B,ns)
    HT <- hotelling$HT
    HT_P <- hotelling$HT.pvalue
    
    data.frame(W = W, 
               W_star = W_star, 
               W_star2 = W_star2,
               W_BC = W_BC,
               W_D = W_D,
               HT = HT,
               W_P = W_P, 
               W_star_P = W_star_P, 
               W_star2_P = W_star2_P,
               W_BC_P = W_BC_P,
               D_P = D_P,
               HT_P = HT_P,
               loggamma = loggamma)
  }else{
    ## He et al. 2020; central limit theorem
    mun <- n/2 * ( (n-p-k-1/2) * log( (n-1-p) * (n-k) / (n-p-k) / (n-1) ) + (k-1) * log((n-1-p)/(n-1)) + p* log( (n-k)/(n-1) ) ) 
    sigma_square <- 1/2 * (log(1- p/(n-1)) - log(1- p/(n-k)) )
    HE <- ( -W / 2 - mun ) / n / sqrt(sigma_square)
    HE_P <- pnorm(HE,lower.tail = FALSE)
    
    data.frame(W = W, 
               W_star = W_star, 
               W_star2 = W_star2,
               W_BC = W_BC,
               W_D = W_D,
               W_P = W_P, 
               W_star_P = W_star_P, 
               W_star2_P = W_star2_P,
               W_BC_P = W_BC_P,
               D_P = D_P,
               HE = HE,
               HE_P = HE_P,
               loggamma = loggamma)
  }
  
}


######################
### example ##########
#####################
library(car)
data(Pottery)
data1 <- Pottery
level <- levels(data1$Site)
X1 <- data1[data1$Site==level[1],][,-1]
X2 <- data1[data1$Site==level[2],][,-1]
X3 <- data1[data1$Site==level[3],][,-1]
X4 <- data1[data1$Site==level[4],][,-1]


y <- list(y1=X1,y2=X2,y3=X3, y4=X4)



normal.example(y)

yAL <- list(y1=X1,y4=X4)
normal.example(yAL)

yAI <- list(y1=X1,y3=X3)
normal.example(yAI)

yAC <- list(y1=X1,y2=X2)
normal.example(yAC)# this will false since sum ni is not greater than p + k+1

yCL <- list(y2=X2,y4=X4)
normal.example(yCL)

### 
# yCI <- list(y2=X2,y3=X3)
# normal.example(yCI) # this will false since sum ni is not greater than p + k+1

#### simulation base on dataset 
mu1 <- colMeans(X1)
mu2 <- colMeans(X2)
mu3 <- colMeans(X3)
mu4 <- colMeans(X4)
mu <- colMeans(data1[,-1])
Sigma <- cov(data1[,-1])

plot(mu1,type = "o",pch=19,ylim=c(0,18))
points(mu3,pch=17)
lines(mu3)

points(mu2,pch=18)
lines(mu2)

points(mu4,pch=15)
lines(mu4)

library(parallel)
library(mvtnorm)
Nsim <- 10000
p <- sapply(y,ncol)[1]
ns <- sapply(y, nrow)
n <- sum(ns)
k <- length(ns)
df <- p*(k-1)
seed <- 123
set.seed(seed)
RNGkind("L'Ecuyer-CMRG")

p.interest <- p*(k-1)

Sigmas <- as.list(numeric(k))
for (ii in 1:k) {
  Sigmas[[ii]] <- Sigma
}

# Mus <- mu

Mus <- matrix(mu,ncol=p,nrow=k,byrow = T)

data.id <- as.list(numeric(Nsim))
for (jj in 1:Nsim){
  data.id[[jj]] <- list(id=jj)
}

res <-  mclapply(1:Nsim, FUN=normal.sim.single,
                 Mus=Mus, Sigmas=Sigmas, ns=ns,
                 data.id=data.id,trace=TRUE,
                 mc.cores = 1,mc.set.seed = TRUE)

if(k==2){
  normal.result.example <- list(W = NULL, W_star = NULL, W_star2 = NULL,
                        W_BC = NULL,W_D = NULL,HT = NULL,
                        W_P = NULL, W_star_P = NULL, W_star2_P = NULL,
                        W_BC_P = NULL,D_P = NULL,HT_P = NULL,loggamma=NULL)
}else{
  normal.result.example <- list(W = NULL, W_star = NULL, W_star2 = NULL,
                        W_BC = NULL,W_D = NULL,W_P = NULL, 
                        W_star_P = NULL, 
                        W_star2_P = NULL,
                        W_BC_P = NULL,
                        D_P = NULL,
                        HE = NULL,
                        HE_P = NULL,
                        loggamma = NULL)
  # list(W = NULL, W_star = NULL, W_star2 = NULL,
  #                     W_BC = NULL,W_D = NULL,
  #                     W_P = NULL, W_star_P = NULL, W_star2_P = NULL,
  #                     W_BC_P = NULL,D_P = NULL,loggamma=NULL)
}

for (i in 1:length(res)){
  for (j in 1:length(normal.result.example)) {
    normal.result.example[[j]] <- cbind(normal.result.example[[j]],res[[i]][[j]])
  }
}



(filename=paste("HDnormal_example_n",sum(ns),"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same",".RData",sep =""))
save(normal.result.example=normal.result.example,file=filename)

#############################################
setwd("C:/Users/HCZ/Desktop/dirmean/example")
size.plot <- function(out,level=0.05){
  ###############################
  ##corrected significance level
  #############################
  
  D_c <-  quantile( out$D_P, probs = level)
  LRT_c <- quantile(out$W_P, probs = level)
  sko1_c <- quantile(out$W_star_P, probs = level)
  sko2_c <- quantile(out$W_star2_P, probs = level)
  BC_c <- quantile(out$W_BC_P, probs = level)
  HE_c <- quantile( (1-out$HE_P), probs = level)
  
  
  
  D_s <-  mean( out$D_P < level)
  LRT_s <- mean(out$W_P < level)
  sko1_s <- mean(out$W_star_P < level)
  sko2_s <- mean(out$W_star2_P < level)
  BC_s <- mean(out$W_BC_P < level)
  HE_s <- mean((1-out$HE_P) < level)
  
  
  #result
  # corrected.alpha <- cbind(D_c,LRT_c,BC_c,sko1_c,sko2_c)
  # size <- cbind(D_s,LRT_s,BC_s,sko1_s,sko2_s)
  # rownames(corrected.alpha) <- c("corrected.alpha")
  # rownames(size) <- c("size")
  # colnames(size) <- c("DT","LRT","BC","sko1","sko2")
  # ###conbine the result
  # rbind( size,corrected.alpha)
  
  #result
  corrected.alpha <- cbind(D_c,HE_c,LRT_c,BC_c,sko1_c,sko2_c)
  size <- cbind(D_s,HE_s,LRT_s,BC_s,sko1_s,sko2_s)
  rownames(corrected.alpha) <- c("corrected.alpha")
  rownames(size) <- c("size")
  colnames(size) <- c("DT","CLT","LRT","BC","sko1","sko2")
  ###conbine the result
  rbind( size,corrected.alpha)
  
}


ns <- c(5 ,2,5,14 )
k <- 4
p <- 5
Nsim <- 10000
(filename=paste("HDnormal_example_n",sum(ns),"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same",".RData",sep =""))
load(file = filename)
size.plot( out=normal.result.example, level = 0.05)



W_P <- normal.result.example$W_P
D_P <- normal.result.example$D_P
BC_P <- normal.result.example$W_BC_P
SKO1_P <- normal.result.example$W_star_P
windows()
par(mar=c(4.5,4.6,0.5,0.5),mfrow=c(1,2))

pdf(paste("section1_same_example_n_",n,"_p",p,"_k",k,"_LRT",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(W_P,xlab = "",ylab="",main = "",cex.axis=2.5,
     ylim=c(0,2000),col = "lightblue")
# hist(W_P,xlab = "",ylab="",main = "",cex.axis=2.5,density=T)
dev.off()

pdf(paste("section1_same_example_n_",n,"_p",p,"_k",k,"_DT",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(D_P,xlab = "",ylab="",main = "",cex.axis=2.5,col = "lightblue")
dev.off()


pdf(paste("section1_same_example_n_",n,"_p",p,"_k",k,"_BC",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(BC_P,xlab = "",ylab="",main = "",cex.axis=2.5)
dev.off()


pdf(paste("section1_same_example_n_",n,"_p",p,"_k",k,"_SKO1",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(SKO1_P,xlab = "",ylab="",main = "",cex.axis=2.5)
dev.off()

