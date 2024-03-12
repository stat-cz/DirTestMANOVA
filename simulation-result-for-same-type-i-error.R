
rm(list=ls())
#######################################################################
#######################################################################
setwd("C:/Users/HCZ/Desktop/dirmean/same/null")

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





######################################################
##' various ratio 
##' fixed sample size 
##' moderate dimension setting
##' with p  =  floor(ni^epsilon)
#######################################################

plot.size <- function(n,case,k=3,epsilons=c(seq(6,22,2))/24,Nsim=10000){
  res <- as.list(numeric(length(epsilons)))
  for (i in 1:length(epsilons)) {
    ns <- rep(n,k)
    
    p <- floor(ns[1]^epsilons[i])
    
    # p.interest <- p*(k-1)
    # load(file = paste("HDnormal_n",n1,"_k",k,"_p",p,"_Nsim",10000,"_meanvector_same",".RData",sep =""))
    
    if(case =="normal"){
      # load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he_sum",".RData",sep =""))
      
      load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
      res[[i]] <- size.plot( out=normal.result, level = 0.05)
    }   
    
    if(case =="normalsum"){
      p <- floor(sum(ns)^epsilons[i])#n1*ratio[i]
      load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he_sum",".RData",sep =""))
      # load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
      res[[i]] <- size.plot( out=normal.result, level = 0.05)
    }   
    
    
    
    
    
    if(case =="t"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
      res[[i]] <- size.plot( out=t.result, level = 0.05)
    }   
    if(case =="skewn"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
      res[[i]] <- size.plot( out=skewn.result, level = 0.05)
    }   
    if(case =="skewt"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
      res[[i]] <- size.plot( out=skewt.result, level = 0.05)
    }   
    if(case =="laplace"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
      res[[i]] <- size.plot( out=laplace.result, level = 0.05)
    }   
    
    
    # load(file = paste("HDnormal_n1",n1,"_n2",n2,"_n3",n3,"_p",p,"_Nsim",10000,"_meanvector_same_group3",".RData",sep =""))
  }
  
  # res
  
  
  
  result <- array(NA,dim = c(2*length(epsilons)+1, ncol(res[[1]])) )
  for (i in 1:length(epsilons)) {
    result[i,] <- res[[i]][1,]
  }
  for (j in (length(epsilons)+2):(2*length(epsilons)+1) ) {
    result[j,] <- res[[j-(length(epsilons)+1) ]][2,]
  }
  
  library(MASS)
  
  
  result <- round(result,digits = 3)
  rownames(result) <- c(floor((ns[1])^epsilons), "NA", floor((ns[1])^epsilons));
  if(case =="normalsum") {rownames(result) <- c(floor(sum(ns)^epsilons), "NA", floor(sum(ns)^epsilons));}
  colnames(result) <- c("DT","CLT","LRT","BC","Sko1","Sko2")
  # result
  library(xtable)
  print(xtable(result,digits = 3))
  
  
  library(reshape2)
  library(ggplot2)
  library(ggsci) 
  # colnames(results) <- c("solid","dashed","dotted","longdash")
  results <- result[1:length(epsilons),]
  dataplot <- melt(results)
  dataplot$Var1 <-  fractions(epsilons)
  colnames(dataplot) <- c("tau","method","typeIerror")
  
  ## create the plot
  # if(j==1) y_lab = 'Corrected power'
  # else y_lab = ''
  
  y_lab = "Empirical size"
  if(case == "normal"&& (n==500 || n==1000) ) y_lab = ""
  if(case == "skewn" || case=="laplace" ) y_lab = ""
  
  # if(j==4|j==5 |j==6) x_lab = expression(p/n)
  # else x_lab = ''
  
  x_lab = expression(tau)
  
  
  cn <- seq(6,22,2)/24
  break2 <- (ns[1]^cn)
  p <- floor(ns[1]^cn)
  if(case=="normalsum") {p <- floor(sum(ns)^cn);break2 <- sum(ns)^cn}#n1*ratio[i]
  
  legend.label <- "none"
  yaxistext <- element_blank()
  if(case == "t"|| case=="normalsum"){legend.label <- c(.1,.75);yaxistext <- element_text(size=25)}
  
  if(case=="normal" && n=="100"){legend.label <- c(.1,.75);yaxistext <- element_text(size=25)}
  if(k > 3){legend.label <- c(.1,.75);yaxistext <- element_text(size=25)}
  plotsize <- ggplot(dataplot, aes(x=tau, y=typeIerror, col=method)) +
    geom_hline(yintercept=0.05, linetype="dashed", color = "grey",linewidth=1)+
    # geom_vline(xintercept = 2/3, linetype="dotted", 
    #            color = "grey", linewidth=1.5)+
    # geom_vline(xintercept = 4/5, linetype="dotted",
    #            color = "grey", linewidth=1.5)+
    geom_line(aes(linetype=method),linewidth=1) + 
    geom_point(aes(shape=method),size=2)  + 
    labs(x = x_lab, y = y_lab)+
    #ylab('Corrected power')+xlab(expression(lambda))+
    scale_color_d3()+
    # scale_color_aaas()+
    scale_x_continuous(breaks= round(cn,3),
                       sec.axis = sec_axis(function(x){if(case=="normalsum"){
                         sum(ns)^x
                       }else ns[1]^x}, 
                                           breaks = break2,
                                           labels = p,
                                           name = expression(p) )) +
    # ylim(c(0,0.2))+
    coord_cartesian(ylim = c(0,.15)) +
    theme_bw()+
    # theme_minimal() +
    theme(plot.title = element_text(size=25,face = "italic"),
          axis.text.x=element_text(size=25),
          axis.text.y = yaxistext,
          axis.title= element_text(size = 25),
          legend.title = element_blank(),
          legend.text = element_text(size=25),
          # legend.text = element_blank(),
          # legend.position = "bottom"
          legend.position = legend.label
    )
  
  print(plotsize)
  
  # ggsave(file=paste("Empirical_type_I_error_ni",n,".png",sep=""),plot=plotsize,width=10,height=10)
  # ggsave(file=paste("Empirical_type_I_error_ni",n,".pdf",sep=""),plot=plotsize,width=10,height=10)
  
  ### for robustness
  ggsave(file=paste("Empirical_type_I_error_",case,"_ni",n,"_k",k,".png",sep=""),plot=plotsize,width=10,height=7)
  ggsave(file=paste("Empirical_type_I_error_",case,"_ni",n,"_k",k,".pdf",sep=""),plot=plotsize,width=10,height=7)
  
}

k <- 3
num <- c(100,500,1000,1500)
epsilons <- c(seq(6,22,2))/24
Nsim <- 10000
n <- 500
case <- "t"



windows()
plot.size(n=100,case="normal",k=30,epsilons=c(seq(6,22,2))/24,Nsim=10000)



######################################################
##' various ratio 
##' fixed sample size 
##' high dimension setting
##' with p  =  kappa * n
##' update 2023-11-21
#######################################################

plot.size.high <- function(n,case,k=3,kappa =  c(0.05,seq(1,10,by=2)/10),Nsim=10000){
  res <- as.list(numeric(length(kappa)))
  for (i in 1:length(kappa)) {
    ns <- rep(n,k)
    
    p <- ns[1] * kappa[i]
    
    # p.interest <- p*(k-1)
    # load(file = paste("HDnormal_n",n1,"_k",k,"_p",p,"_Nsim",10000,"_meanvector_same",".RData",sep =""))
    
    if(case =="normal"){
      # load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he_sum",".RData",sep =""))
      
      load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same",".RData",sep =""))
      res[[i]] <- size.plot( out=normal.result, level = 0.05)
    }   
    
    if(case =="normalsum"){
      p <- floor(sum(ns)^epsilons[i])#n1*ratio[i]
      load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he_sum",".RData",sep =""))
      # load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
      res[[i]] <- size.plot( out=normal.result, level = 0.05)
    }   
    
    
    if(case=="extreme"){
      load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_extreme",".RData",sep =""))
      # load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
      res[[i]] <- size.plot( out=normal.result, level = 0.05)
      
    }
    
    
    if(case =="t"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same",".RData",sep =""))
      res[[i]] <- size.plot( out=t.result, level = 0.05)
    }   
    if(case =="skewn"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same",".RData",sep =""))
      res[[i]] <- size.plot( out=skewn.result, level = 0.05)
    }   
    if(case =="skewt"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same",".RData",sep =""))
      res[[i]] <- size.plot( out=skewt.result, level = 0.05)
    }   
    if(case =="laplace"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same",".RData",sep =""))
      res[[i]] <- size.plot( out=laplace.result, level = 0.05)
    }   
    
    
    # load(file = paste("HDnormal_n1",n1,"_n2",n2,"_n3",n3,"_p",p,"_Nsim",10000,"_meanvector_same_group3",".RData",sep =""))
  }
  
  # res
  
  
  
  result <- array(NA,dim = c(2*length(kappa)+1, ncol(res[[1]])) )
  for (i in 1:length(kappa)) {
    result[i,] <- res[[i]][1,]
  }
  for (j in (length(kappa)+2):(2*length(kappa)+1) ) {
    result[j,] <- res[[j-(length(kappa)+1) ]][2,]
  }
  
  library(MASS)
  
  
  result <- round(result,digits = 3)
  rownames(result) <- c(ns[1]*kappa, "NA", ns[1]*kappa);
  #if(case =="extreme") {rownames(result) <- c((sum(ns)*kappa), "NA", (sum(ns)*kappa));}
  colnames(result) <- c("DT","CLT","LRT","BC","Sko1","Sko2")
  # result
  library(xtable)
  print(xtable(result,digits = 3))
  
  
  library(reshape2)
  library(ggplot2)
  library(ggsci) 
  # colnames(results) <- c("solid","dashed","dotted","longdash")
  results <- result[1:length(kappa),]
  dataplot <- melt(results)
  dataplot$Var1 <-  kappa
  colnames(dataplot) <- c("kappa","method","typeIerror")
  
  ## create the plot
  # if(j==1) y_lab = 'Corrected power'
  # else y_lab = ''
  
  y_lab = "Empirical size"
  if(case == "normal"&& (n==500 || n==1000) ) y_lab = ""
  if(case == "skewn" || case=="laplace" ) y_lab = ""
  
  # if(j==4|j==5 |j==6) x_lab = expression(p/n)
  # else x_lab = ''
  
  x_lab = expression(kappa)
  
  
  cn <-c(0.05,seq(1,10,2)/10)
  if(case=="extreme") cn <- c(1,1.5,2,2.5,2.9)
  # p <- floor(ns[1]^cn)
  # if(case=="normalsum") p <- floor(sum(ns)^cn)#n1*ratio[i]
  
  legend.label <- "none"
  yaxistext <- element_blank()
  if(case == "t"|| case=="extreme"){legend.label <- c(.1,.75);yaxistext <- element_text(size=25)}
  
  if(case=="normal" && n=="100"){legend.label <- c(.1,.75);yaxistext <- element_text(size=25)}
  if(k > 3){legend.label <- c(.1,.75);yaxistext <- element_text(size=25)}
  plotsize <- ggplot(dataplot, aes(x=kappa, y=typeIerror, col=method)) +
    geom_hline(yintercept=0.05, linetype="dashed", color = "grey",linewidth=1)+
    # geom_vline(xintercept = 2/3, linetype="dotted", 
    #            color = "grey", linewidth=1.5)+
    # geom_vline(xintercept = 4/5, linetype="dotted",
    #            color = "grey", linewidth=1.5)+
    geom_line(aes(linetype=method),linewidth=1) + 
    geom_point(aes(shape=method),size=2)  + 
    labs(x = x_lab, y = y_lab)+
    #ylab('Corrected power')+xlab(expression(lambda))+
    scale_color_d3()+
    # scale_color_aaas()+
    scale_x_continuous(breaks= round(cn,3), #if cn is (0.05,0.1,0.3,....) we should delete the first element.
                       sec.axis = sec_axis(~.*n,breaks = cn*n,
                                           name = expression(p))) +
    # ylim(c(0,0.2))+ ,
    # sec.axis = sec_axis(c(1,2,3,4,5,6,7,8,9),
    #                     name = 'Test positive rate')
    coord_cartesian(ylim = c(0,.15)) +
    theme_bw()+
    # theme_minimal() +
    theme(plot.title = element_text(size=25,face = "italic"),
          axis.text.x=element_text(size=25),
          axis.text.y = yaxistext,
          axis.title= element_text(size = 25),
          legend.title = element_blank(),
          legend.text = element_text(size=25),
          # legend.text = element_blank(),
          # legend.position = "bottom"
          legend.position = legend.label
    )
  
  print(plotsize)
  
  # ggsave(file=paste("Empirical_type_I_error_ni",n,".png",sep=""),plot=plotsize,width=10,height=10)
  # ggsave(file=paste("Empirical_type_I_error_ni",n,".pdf",sep=""),plot=plotsize,width=10,height=10)
  
  ### for robustness
  if(case=="extreme"){
    ggsave(file=paste("Empirical_type_I_error_",case,"_ni",n,"_k",k,"_extreme.png",sep=""),plot=plotsize,width=10,height=7)
    ggsave(file=paste("Empirical_type_I_error_",case,"_ni",n,"_k",k,"_extreme.pdf",sep=""),plot=plotsize,width=10,height=7)
    
  }else{
    ggsave(file=paste("Empirical_type_I_error_",case,"_ni",n,"_k",k,"_high.png",sep=""),plot=plotsize,width=10,height=7)
    ggsave(file=paste("Empirical_type_I_error_",case,"_ni",n,"_k",k,"_high.pdf",sep=""),plot=plotsize,width=10,height=7)
    
  }
  
}

k <- 3
num <- c(100,500,1000,1500)
kappa <- seq(1,10,2)/10
Nsim <- 10000
n <- 100
case <- "normal"



windows()
plot.size.high(n=100,case="t",k=3,kappa =  c(0.05,seq(1,10,by=2)/10),Nsim=10000)

plot.size.high(n=100,case="extreme",k=3,kappa = c(1,1.5,2,2.5,2.9),Nsim=10000)

#############################################################################

res <- as.list(numeric(length(epsilons)))
for (i in 1:length(epsilons)) {
  ns <- rep(n,k)

  p <- floor(ns[1]^epsilons[i])
  
  # p.interest <- p*(k-1)
  # load(file = paste("HDnormal_n",n1,"_k",k,"_p",p,"_Nsim",10000,"_meanvector_same",".RData",sep =""))
 
  if(case =="normal"){
    # load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he_sum",".RData",sep =""))
    
    load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
    res[[i]] <- size.plot( out=normal.result, level = 0.05)
  }   
  
  if(case =="normalsum"){
    p <- floor(sum(ns)^epsilons[i])#n1*ratio[i]
    load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he_sum",".RData",sep =""))
    # load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
    res[[i]] <- size.plot( out=normal.result, level = 0.05)
  }   
  
  
  
  if(case =="t"){
    load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
    res[[i]] <- size.plot( out=t.result, level = 0.05)
  }   
  if(case =="skewn"){
    load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
    res[[i]] <- size.plot( out=skewn.result, level = 0.05)
  }   
  if(case =="skewt"){
    load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
    res[[i]] <- size.plot( out=skewt.result, level = 0.05)
  }   
  if(case =="laplace"){
    load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
    res[[i]] <- size.plot( out=laplace.result, level = 0.05)
  }   
  
  
  # load(file = paste("HDnormal_n1",n1,"_n2",n2,"_n3",n3,"_p",p,"_Nsim",10000,"_meanvector_same_group3",".RData",sep =""))
}

# res



result <- array(NA,dim = c(2*length(epsilons)+1, ncol(res[[1]])) )
for (i in 1:length(epsilons)) {
  result[i,] <- res[[i]][1,]
}
for (j in (length(epsilons)+2):(2*length(epsilons)+1) ) {
  result[j,] <- res[[j-(length(epsilons)+1) ]][2,]
}

library(MASS)


result <- round(result,digits = 3)
rownames(result) <- c(floor((ns[1])^epsilons), "NA", floor((ns[1])^epsilons));
if(case =="normalsum") {rownames(result) <- c(floor(sum(ns)^epsilons), "NA", floor(sum(ns)^epsilons));}
colnames(result) <- c("DT","CLT","LRT","BC","Sko1","Sko2")
result
library(xtable)
xtable(result,digits = 3)


library(reshape2)
library(ggplot2)
library(ggsci) 
# colnames(results) <- c("solid","dashed","dotted","longdash")
results <- result[1:length(epsilons),]
dataplot <- melt(results)
dataplot$Var1 <-  fractions(epsilons)
colnames(dataplot) <- c("tau","method","typeIerror")

## create the plot
# if(j==1) y_lab = 'Corrected power'
# else y_lab = ''

y_lab = "Empirical size"
if(case == "normal"&& (n==500 || n==1000) ) y_lab = ""
if(case == "skewn" || case=="laplace" ) y_lab = ""

# if(j==4|j==5 |j==6) x_lab = expression(p/n)
# else x_lab = ''

x_lab = expression(tau)


cn <- seq(6,22,2)/24

legend.label <- "none"
if(case == "t"|| case=="normalsum") legend.label <- c(.1,.85)
if(case=="normal" && n=="100") legend.label <- c(.1,.85)
if(k > 3) legend.label <- c(.1,.85)
plotsize <- ggplot(dataplot, aes(x=tau, y=typeIerror, col=method)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "grey",linewidth=1)+
  # geom_vline(xintercept = 2/3, linetype="dotted", 
  #            color = "grey", linewidth=1.5)+
  # geom_vline(xintercept = 4/5, linetype="dotted",
  #            color = "grey", linewidth=1.5)+
  geom_line(aes(linetype=method),linewidth=1) + 
  geom_point(aes(shape=method),size=2)  + 
  labs(x = x_lab, y = y_lab)+
  #ylab('Corrected power')+xlab(expression(lambda))+
  scale_color_d3()+
  # scale_color_aaas()+
  scale_x_continuous(breaks= round(cn,3)) +
  # ylim(c(0,0.2))+
  coord_cartesian(ylim = c(0,.2)) +
  theme_bw()+
  # theme_minimal() +
  theme(plot.title = element_text(size=25,face = "italic"),
        axis.text=element_text(size=25),
        axis.title= element_text(size = 25),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        # legend.text = element_blank(),
        # legend.position = "bottom"
        legend.position = legend.label
        )

plotsize

# ggsave(file=paste("Empirical_type_I_error_ni",n,".png",sep=""),plot=plotsize,width=10,height=10)
# ggsave(file=paste("Empirical_type_I_error_ni",n,".pdf",sep=""),plot=plotsize,width=10,height=10)

### for robustness
ggsave(file=paste("Empirical_type_I_error_",case,"_ni",n,"_k",k,".png",sep=""),plot=plotsize,width=10,height=7)
ggsave(file=paste("Empirical_type_I_error_",case,"_ni",n,"_k",k,".pdf",sep=""),plot=plotsize,width=10,height=7)



############################
### histogram of p-value
############################
### k = 30
k <- 30
ni <- 100
ns <- rep(ni,k)
Nsim <- 10000
p <- floor(ni^(22/24))

load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he",".RData",sep =""))
size.plot( out=normal.result, level = 0.05)

W_P <- normal.result$W_P
D_P <- normal.result$D_P
BC_P <- normal.result$W_BC_P

pdf(paste("section1_same_he_ni_",ni,"_p",p,"_k",k,"_LRT",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(W_P,xlab = "",ylab="",main = "",cex.axis=2.5,ylim = c(0,1500))
dev.off()

pdf(paste("section1_same_he_ni_",ni,"_p",p,"_k",k,"_BC",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(BC_P,xlab = "",ylab="",main = "",cex.axis=2.5)
dev.off()

pdf(paste("section1_same_he_ni_",ni,"_p",p,"_k",k,"_DT",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(D_P,xlab = "",ylab="",main = "",cex.axis=2.5)
dev.off()

#### k=3 ni=100
k <- 3
ni <- 100
ns <- rep(ni,k)
Nsim <- 10000
p <- floor(sum(ns)^(23/24))

load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he_sum",".RData",sep =""))
size.plot( out=normal.result, level = 0.05)

W_P <- normal.result$W_P
D_P <- normal.result$D_P
BC_P <- normal.result$W_BC_P

pdf(paste("section1_same_he_sum_ni_",ni,"_p",p,"_k",k,"_LRT",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(W_P,xlab = "",ylab="",main = "",cex.axis=2.5)
dev.off()

pdf(paste("section1_same_he_sum_ni_",ni,"_p",p,"_k",k,"_BC",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(BC_P,xlab = "",ylab="",main = "",cex.axis=2.5)
dev.off()

pdf(paste("section1_same_he_sum_ni_",ni,"_p",p,"_k",k,"_DT",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(D_P,xlab = "",ylab="",main = "",cex.axis=2.5)
dev.off()


#### k=3 ni=500
k <- 3
ni <- 500
ns <- rep(ni,k)
Nsim <- 10000
p <- floor(sum(ns)^(23/24))

load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_same_he_sum",".RData",sep =""))
size.plot( out=normal.result, level = 0.05)

W_P <- normal.result$W_P
D_P <- normal.result$D_P

pdf(paste("section1_same_he_sum_ni_",ni,"_p",p,"_k",k,"_LRT",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(W_P,xlab = "",ylab="",main = "",cex.axis=2.5)
dev.off()

pdf(paste("section1_same_he_sum_ni_",ni,"_p",p,"_k",k,"_DT",".pdf",sep=""),width=10,height=10)
par(mar=c(4.5,4.6,0.5,0.5))# (bottom,left,top,right)
hist(D_P,xlab = "",ylab="",main = "",cex.axis=2.5)
dev.off()


