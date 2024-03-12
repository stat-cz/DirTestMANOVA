

##################      result for full model #######################
setwd("C:/Users/HCZ/Desktop/dirmean/full/null")
#######################################################################
#######################################################################

size.plot <- function(out,level=0.05){
  ###############################
  ##corrected significance level
  #############################
  
  D_c <-  quantile( out$D_P, probs = level)
  LRT_c <- quantile(out$W_P, probs = level)
  sko1_c <- quantile(out$W_star_P, probs = level)
  sko2_c <- quantile(out$W_star2_P, probs = level)
  # BC_c <- quantile(out$W_BC_P, probs = level)
  if(k==2) BF_c <- quantile(out$BF_P, probs = level) else BF_c <- NA
  
  
  D_s <-  mean( out$D_P < level)
  LRT_s <- mean(out$W_P < level)
  sko1_s <- mean(out$W_star_P < level)
  sko2_s <- mean(out$W_star2_P < level)
  # BC_s <- mean(out$W_BC_P < level)
  if(k==2) BF_s <- mean(out$BF_P < level) else BF_s <- NA
  
  #result
  corrected.alpha <- cbind(D_c,BF_c,LRT_c,sko1_c,sko2_c)
  size <- cbind(D_s,BF_s,LRT_s,sko1_s,sko2_s)
  rownames(corrected.alpha) <- c("corrected.alpha")
  rownames(size) <- c("size")
  colnames(size) <- c("DT","BF","LRT","sko1","sko2")
  ###conbine the result
  rbind( size,corrected.alpha)
}





######################################################
##' various ratio 
##' fixed sample size 
##' with p  =  ceilling(ni^epsilon)
#######################################################

plot.size0 <- function(n,case,k=2,epsilons =c(seq(6,22,1))/24,Nsim = 10000){
  if(k==5 && n==1000) epsilons <- c(seq(6,20,1))/24
  res <- as.list(numeric(length(epsilons)))
  for (i in 1:length(epsilons)) {
    ns <- rep(n,k)
    p <- ceiling(ns[1]^epsilons[i])#n1*ratio[i]
    p.interest <- p*(k-1)
    # # load(file = paste("HDnormal_n",n1,"_k",k,"_p",p,"_Nsim",10000,"_meanvector_same",".RData",sep =""))
    # load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_autoregression","_meanvector_full_he",".RData",sep =""))
    # 
    # # load(file = paste("HDnormal_n1",n1,"_n2",n2,"_n3",n3,"_p",p,"_Nsim",10000,"_meanvector_same_group3",".RData",sep =""))
    # res[[i]] <- size.plot( out=normal.result, level = 0.05)
    # 
    
    if(case =="normal"){
      Nsim=10000
      if(n==1000 && (p== 422  || p==563) ) Nsim =5000
      load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_autoregression","_meanvector_full_he",".RData",sep =""))
      res[[i]] <- size.plot( out=normal.result, level = 0.05)
    }    
    if(case =="t"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
      res[[i]] <- size.plot( out=t.result, level = 0.05)
    }   
    if(case =="skewn"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
      res[[i]] <- size.plot( out=skewn.result, level = 0.05)
    }   
    if(case =="skewt"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
      res[[i]] <- size.plot( out=skewt.result, level = 0.05)
    }   
    if(case =="laplace"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
      res[[i]] <- size.plot( out=laplace.result, level = 0.05)
    }   
    
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
  rownames(result) <- c(ceiling(ns[1]^epsilons), "NA", ceiling(ns[1]^epsilons));
  colnames(result) <- c("DT","BF","LRT","Sko1","Sko2")
  # result
  if(k>2) result <- result[,-2]
  library(xtable)
  if(case=="normal"){
    print(xtable(result,digits = 3))
  }else{
    print(xtable(result[c(seq(1,length(epsilons),by=2),(length(epsilons)+1),seq((length(epsilons)+2),(2*length(epsilons)+1),by=2)),],digits = 3))
  } 
  
  
  
  library(reshape2)
  library(ggplot2)
  library(ggsci) 
  
  results <- result[1:length(epsilons),]
  # colnames(results) <- c("solid","longdash","dashed","dotted","dotdash") #“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”
  # if(k>2) colnames(results) <- c("solid","dashed","dotted","longdash") #“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”
  dataplot <- melt(results)
  dataplot$Var1 <-  fractions(epsilons)
  dataplot$p <- ceiling(n^dataplot$Var1)
  colnames(dataplot) <- c("tau","method","typeIerror","p")
  
  ## create the plot
  # if(j==1) y_lab = 'Corrected power'
  # else y_lab = ''
  
  y_lab = "Empirical size"
  if(case == "normal"&& (n==500 || n==1000) ) y_lab = ""
  if(case == "skewn" || case=="laplace" ) y_lab = ""
  # if(j==4|j==5 |j==6) x_lab = expression(p/n)
  # else x_lab = ''
  
  x_lab = expression(p)
  
  cn <- seq(6,22,2)/24;ceiling(n^cn)
  axis2 <- c(4,7,10,15,22,32,47,69)
  legend.label <- "none"
  yaxistext <- element_blank()
  
  if(case == "t") {legend.label <- c(.1,.8);yaxistext <- element_text(size=25)}
  if(case=="normal" && n=="100"){legend.label <- c(.1,.8);yaxistext <- element_text(size=25)}
  if(k > 3){legend.label <- c(.1,.8);yaxistext <- element_text(size=25)}
  # windows()
  plotsize <- ggplot(dataplot, aes(x=p, y=typeIerror, col=method)) +
    geom_hline(yintercept=0.05, linetype="dashed", color = "grey",linewidth=1)+
    # geom_vline(xintercept = 2/3, linetype="dotted", 
    #            color = "grey", linewidth=1.5)+
    # geom_vline(xintercept = 1/2, linetype="dotted", 
    #            color = "grey", linewidth=1.5)+
    # geom_vline(xintercept = 5/6, linetype="dotted", 
    #            color = "grey", linewidth=1.5)+
    geom_line(aes(linetype=method),linewidth=1) + 
    geom_point(aes(shape=method),size=2)  +
    labs(x = x_lab, y = y_lab)+
    #ylab('Corrected power')+xlab(expression(lambda))+
    scale_color_d3()+
    # scale_color_aaas()+
    scale_x_continuous(breaks= axis2,
                       sec.axis = sec_axis(function(x) log(x)/log(n),
                                           breaks = round(log(axis2)/log(n),digits = 3),
                                           labels = c(0.25,0.42,0.50,0.59,0.67,0.75,0.83,0.91),
                                           name = expression(tau) )) +
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
  
  ### for robustness
  ggsave(file=paste("Empirical_type_I_error_",case,"_full0_ni",n,"_k",k,".png",sep=""),plot=plotsize,width=10,height=7)
  ggsave(file=paste("Empirical_type_I_error_",case,"_full0_ni",n,"_k",k,".pdf",sep=""),plot=plotsize,width=10,height=7)
  
  
}


######################################################
##' various ratio 
##' fixed sample size 
##' with p  =  ceilling(ni^epsilon)
##' 2023-11-10
#######################################################

plot.size <- function(n,case,k=2,epsilons =c(seq(6,22,1))/24,Nsim = 10000){
  if(k==5 && n==1000) epsilons <- c(seq(6,20,1))/24
  res <- as.list(numeric(length(epsilons)))
  for (i in 1:length(epsilons)) {
    ns <- rep(n,k)
    p <- ceiling(ns[1]^epsilons[i])#n1*ratio[i]
    p.interest <- p*(k-1)
    # # load(file = paste("HDnormal_n",n1,"_k",k,"_p",p,"_Nsim",10000,"_meanvector_same",".RData",sep =""))
    # load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_autoregression","_meanvector_full_he",".RData",sep =""))
    # 
    # # load(file = paste("HDnormal_n1",n1,"_n2",n2,"_n3",n3,"_p",p,"_Nsim",10000,"_meanvector_same_group3",".RData",sep =""))
    # res[[i]] <- size.plot( out=normal.result, level = 0.05)
    # 
    
    if(case =="normal"){
      Nsim=10000
      if(n==1000 && (p== 422  || p==563) ) Nsim =5000
      load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_autoregression","_meanvector_full_he",".RData",sep =""))
      res[[i]] <- size.plot( out=normal.result, level = 0.05)
    }    
    if(case =="t"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
      res[[i]] <- size.plot( out=t.result, level = 0.05)
    }   
    if(case =="skewn"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
      res[[i]] <- size.plot( out=skewn.result, level = 0.05)
    }   
    if(case =="skewt"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
      res[[i]] <- size.plot( out=skewt.result, level = 0.05)
    }   
    if(case =="laplace"){
      load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
      res[[i]] <- size.plot( out=laplace.result, level = 0.05)
    }   
    
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
  rownames(result) <- c(ceiling(ns[1]^epsilons), "NA", ceiling(ns[1]^epsilons));
  colnames(result) <- c("DT","BF","LRT","Sko1","Sko2")
  # result
  if(k>2) result <- result[,-2]
  library(xtable)
  if(case=="normal"){
    print(xtable(result,digits = 3))
  }else{
    print(xtable(result[c(seq(1,length(epsilons),by=2),(length(epsilons)+1),seq((length(epsilons)+2),(2*length(epsilons)+1),by=2)),],digits = 3))
  } 
  
  
  
  library(reshape2)
  library(ggplot2)
  library(ggsci) 
  
  results <- result[1:length(epsilons),]
  # colnames(results) <- c("solid","longdash","dashed","dotted","dotdash") #“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”
  # if(k>2) colnames(results) <- c("solid","dashed","dotted","longdash") #“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”
  dataplot <- melt(results)
  dataplot$Var1 <-  fractions(epsilons)
  # dataplot$p <- ceiling(n^dataplot$Var1)
  colnames(dataplot) <- c("tau","method","typeIerror")
  # colnames(dataplot) <- c("tau","method","typeIerror","p")
  
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
  yaxistext <- element_blank()
  
  if(case == "t") {legend.label <- c(.1,.8);yaxistext <- element_text(size=25)}
  if(case=="normal" && n=="100"){legend.label <- c(.1,.8);yaxistext <- element_text(size=25)}
  if(k > 3){legend.label <- c(.1,.8);yaxistext <- element_text(size=25)}
  # windows()
  plotsize <- ggplot(dataplot, aes(x=tau, y=typeIerror, col=method)) +
    geom_hline(yintercept=0.05, linetype="dashed", color = "grey",linewidth=1)+
    # geom_vline(xintercept = 2/3, linetype="dotted", 
    #            color = "grey", linewidth=1.5)+
    # geom_vline(xintercept = 1/2, linetype="dotted", 
    #            color = "grey", linewidth=1.5)+
    # geom_vline(xintercept = 5/6, linetype="dotted", 
    #            color = "grey", linewidth=1.5)+
    geom_line(aes(linetype=method),linewidth=1) + 
    geom_point(aes(shape=method),size=2)  +
    labs(x = x_lab, y = y_lab)+
    #ylab('Corrected power')+xlab(expression(lambda))+
    scale_color_d3()+
    # scale_color_aaas()+
    scale_x_continuous(breaks= round(cn,digits=3),
                       sec.axis = sec_axis(function(x) n^x,
                                           breaks = (n^cn),
                                           labels = ceiling(n^cn),
                                           name = expression(p) )) +
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
  
  ### for robustness
  ggsave(file=paste("Empirical_type_I_error_",case,"_full0_ni",n,"_k",k,".png",sep=""),plot=plotsize,width=10,height=7)
  ggsave(file=paste("Empirical_type_I_error_",case,"_full0_ni",n,"_k",k,".pdf",sep=""),plot=plotsize,width=10,height=7)
  
  
}


k <- 2
num <- c(100,500,1000)
Nsim <- 10000

## for ni=100
epsilons <- c(seq(6,22,1))/24
n <- 100
case <- "normal"


plot.size(n=500,case="skewt",k=2,epsilons =c(seq(6,22,1))/24,Nsim = 10000)

# 
# ### 
# epsilons <- c(seq(6,22,1))/24
# n <- 500
# 
# ######## for ni = 1000
# # epsilons <- c(6:15/24, c(61,62)/96, 16:20/24)
# 
# epsilons <- 6:22/24
# n <- 1000
# 
# 
# ######## for ni = 2000
# epsilons <- c(6:16/24)
# n <- 2000
# 
# 
# ######## for ni = 4000
# epsilons <- c(12:15/24)
# n <- 4000





res <- as.list(numeric(length(epsilons)))
for (i in 1:length(epsilons)) {
  ns <- rep(n,k)
  p <- ceiling(ns[1]^epsilons[i])#n1*ratio[i]
  p.interest <- p*(k-1)
  # # load(file = paste("HDnormal_n",n1,"_k",k,"_p",p,"_Nsim",10000,"_meanvector_same",".RData",sep =""))
  # load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_autoregression","_meanvector_full_he",".RData",sep =""))
  # 
  # # load(file = paste("HDnormal_n1",n1,"_n2",n2,"_n3",n3,"_p",p,"_Nsim",10000,"_meanvector_same_group3",".RData",sep =""))
  # res[[i]] <- size.plot( out=normal.result, level = 0.05)
  # 
  
  if(case =="normal"){
    Nsim=10000
    if(n==1000 && (p== 422  || p==563) ) Nsim =5000
    load(file=paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_autoregression","_meanvector_full_he",".RData",sep =""))
    res[[i]] <- size.plot( out=normal.result, level = 0.05)
  }    
  if(case =="t"){
    load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
    res[[i]] <- size.plot( out=t.result, level = 0.05)
  }   
  if(case =="skewn"){
    load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
    res[[i]] <- size.plot( out=skewn.result, level = 0.05)
  }   
  if(case =="skewt"){
    load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
    res[[i]] <- size.plot( out=skewt.result, level = 0.05)
  }   
  if(case =="laplace"){
    load(file=paste("HD",case,"_n",ns[1],"_k",k,"_p",p,"_Nsim",Nsim,"_meanvector_full_he",".RData",sep =""))
    res[[i]] <- size.plot( out=laplace.result, level = 0.05)
  }   
  
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
rownames(result) <- c(ceiling(ns[1]^epsilons), "NA", ceiling(ns[1]^epsilons));
colnames(result) <- c("DT","BF","LRT","Sko1","Sko2")
result
if(k>2) result <- result[,-2]
library(xtable)
xtable(result[c(seq(1,length(epsilons),by=2),(length(epsilons)+1),seq((length(epsilons)+2),(2*length(epsilons)+1),by=2)),],digits = 3)
if(case=="normal") xtable(result,digits = 3)



library(reshape2)
library(ggplot2)
library(ggsci) 

results <- result[1:length(epsilons),]
# colnames(results) <- c("solid","longdash","dashed","dotted","dotdash") #“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”
# if(k>2) colnames(results) <- c("solid","dashed","dotted","longdash") #“solid”, “dashed”, “dotted”, “dotdash”, “longdash”, “twodash”
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
if(case == "t") legend.label <- c(.1,.85)
if(case=="normal" && n=="100") legend.label <- c(.1,.85)
if(k > 3) legend.label <- c(.1,.85)
# windows()
plotsize <- ggplot(dataplot, aes(x=tau, y=typeIerror, col=method)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "grey",linewidth=1)+
  # geom_vline(xintercept = 2/3, linetype="dotted", 
  #            color = "grey", linewidth=1.5)+
  # geom_vline(xintercept = 1/2, linetype="dotted", 
  #            color = "grey", linewidth=1.5)+
  # geom_vline(xintercept = 5/6, linetype="dotted", 
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

# ggsave(file=paste("Empirical_type_I_error_full_ni",n,".png",sep=""),plot=plotsize,width=10,height=10)
# 
# ggsave(file=paste("Empirical_type_I_error_full_ni",n,".pdf",sep=""),plot=plotsize,width=10,height=10)


### for robustness
ggsave(file=paste("Empirical_type_I_error_",case,"_full0_ni",n,"_k",k,".png",sep=""),plot=plotsize,width=10,height=7)
ggsave(file=paste("Empirical_type_I_error_",case,"_full0_ni",n,"_k",k,".pdf",sep=""),plot=plotsize,width=10,height=7)



####################################################################
####################################
####### for ni = 1000
#### small distance for tau
####### c(61,62,63)/96
#####################################
cn <- c(6:15/24, c(61,62)/96, 16:20/24)
ns <- rep(1000,k)
p <- ceiling(ns[1]^cn[1]);p.interest <-  p*(k-1);p
load(file = paste("HDnormal_n",ns[1],"_k",k,"_p",p,"_Nsim",10000,"_autoregression","_meanvector_full_he",".RData",sep =""))
size.plot( out=normal.result, level = 0.05)




