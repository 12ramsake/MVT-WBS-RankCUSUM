

library(MASS)
library(xtable)


#results directory
dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/WBS_SIM/"
setwd(dirr)


#simulation parameters

#N
Ns<-c(1000,2500,5000)
#Number of intervals, changed from 3 to 8
mod=100
numInts<-floor(log(Ns))*mod

#simulation size, num repetitions
sim.size=100


# 2,3.5 cps
thetas<-list(c(.333,.666), 
             c(.25,.5,.75),
             c(1/6,2/6,3/6,4/6,5/6))

ds=c(2,3,5,10)
# ds=c(2)



#95 of BB dont dchange for bonff
thresh<-1

distributions1=1:3
names(distributions1)<-c("Normal", "Cauchy", "Skew Normal")
##Create Parameter Vector

numUniqueRuns<-length(Ns)*length(thetas)*length(ds)*length(distributions1)

paramterIndices<-matrix(0,ncol=4,nrow=numUniqueRuns)
curr=0
for(i1 in 1:length(Ns)){
  for(i2 in 1:length(thetas)){
    for(i3 in 1:length(distributions1)){
      for(i4 in 1:length(ds)){
        curr=curr+1
        paramterIndices[curr,]=c(i1,i2,i3,i4)
      }
    }
  }
}







#simulation functions

#generates iid d-normal, sigmaxI Rv
normalMaster<-function(n,d,sigmaSq){mvtnorm::rmvnorm(n,sigma=diag(sigmaSq*rep(1,d)))}

norm1<-function(n,d){normalMaster(n,d,1)}
norm2<-function(n,d){normalMaster(n,d,2.5)}
norm3<-function(n,d){normalMaster(n,d,4)}
norm4<-function(n,d){normalMaster(n,d,2.25)}
norm5<-function(n,d){normalMaster(n,d,5)}

normals<-list(norm1,norm2,norm3,norm4,norm5,norm1)


#generates iid d-cauchy, sigmaxI scale
cauchyMaster<-function(n,d,sigmaSq){replicate(d,rcauchy(n,scale=sigmaSq))}


cauchy1<-function(n,d){cauchyMaster(n,d,1)}
cauchy2<-function(n,d){cauchyMaster(n,d,2.5)}
cauchy3<-function(n,d){cauchyMaster(n,d,4)}
cauchy4<-function(n,d){cauchyMaster(n,d,2.25)}
cauchy5<-function(n,d){cauchyMaster(n,d,5)}

cauchys<-list(cauchy1,cauchy2,cauchy3,cauchy4,cauchy5,cauchy1)

skewNormalMaster<-function(n,d,sigmaSq,skewParam){
  rmsn(n, dp=cp2dp(list(mean=rep(0,d), 
                        var.cov=diag(sigmaSq*rep(1,d)), 
                        gamma1=skewParam*rep(1,d)/d), "SN"))
}

skewNorm1<-function(n,d){skewNormalMaster(n,d,1,.1)}
skewNorm2<-function(n,d){skewNormalMaster(n,d,2.5,.1)}
skewNorm3<-function(n,d){skewNormalMaster(n,d,4,.1)}
skewNorm4<-function(n,d){skewNormalMaster(n,d,2.25,.1)}
skewNorm5<-function(n,d){skewNormalMaster(n,d,5,.1)}

skewNormals<-list(skewNorm1,skewNorm2,skewNorm3,skewNorm4,skewNorm5,skewNorm1)

distributions<-list(normals,cauchys,skewNormals)
names(distributions)=c("Normal","Cauchy","Skew Normal")


getScwartzCP<-function(alpha=1,results){
  
  
  applySCH<-function(models){
    #get indidvidual criteria for a set of cp
    getsSic<-function(cp,sighatSq){
      
      
      #if at least 1 cp
      if(length(cp)>=1){
        
        sSic<-(N/2)*log(sighatSq)+length(cp)*(log(N))^alpha
      }
      else{
        sSic<-(N/2)*log(sighatSq)
      }
      
      return(sSic)
    }
    
    #get ssic for all ammounts of cp
    
    sSic<-getsSic(models[[1]]$cp,models[[1]]$sigSq)
    
    for(j in 2:length(models)){
      
      sSic<-c(sSic,getsSic(models[[j]]$cp,models[[j]]$sigSq))
      
    }
    
    minVal<-which.min(sSic)
    
    return(models[[minVal]]$cp)  
    
  }
  
  
  return(lapply(results,applySCH))
}

# resultsHS_Sch<-getScwartzCP(.1,resultsHS)
# resultsSpat_Sch<-getScwartzCP(.1,resultsSpat)
# resultsMahal_Sch<-getScwartzCP(.1,resultsMahal)
# resultsMahal75_Sch<-getScwartzCP(.1,resultsMahal75)


getSummaryL<-function(results,N,theta,distt,d,numInt,thresh,depth){
  
  
  #l-lhat
  
  l<-length(theta)
  print(l)
  lhat<-lapply(results,length)
  lhat<-unlist(lhat)
  llhat<-lhat-l
  summ<-c(summary(llhat),mean(llhat^2),median(abs(llhat)))
  names(summ)=c(names(summ)[1:6],"RMSE","MAD")
  
  return(median(llhat))
  
}
getSummaryL_MSE<-function(results,N,theta,distt,d,numInt,thresh,depth){
  
  
  #l-lhat
  
  l<-length(theta)
  # print(l)
  lhat<-lapply(results,length)
  lhat<-unlist(lhat)
  llhat<-lhat-l
  
  return(mean(llhat^2))
  
}


getmedians<-Vectorize(function(alpha){
  
  resultsHS_Sch<-getScwartzCP(alpha,resultsHS)
  resultsSpat_Sch<-getScwartzCP(alpha,resultsSpat)
  resultsMahal_Sch<-getScwartzCP(alpha,resultsMahal)
  resultsMahal75_Sch<-getScwartzCP(alpha,resultsMahal75)
  
  summHSL<-getSummaryL(resultsHS_Sch,N,theta,distr,d,numInt,thresh,depth="hs")
  summSpatL<-getSummaryL(resultsSpat_Sch,N,theta,distr,d,numInt,thresh,depth="spat")
  summMahalL<-getSummaryL(resultsMahal_Sch,N,theta,distr,d,numInt,thresh,depth="mahal")
  summMahal75L<-getSummaryL(resultsMahal75_Sch,N,theta,distr,d,numInt,thresh,depth="mahal75")
  
  return(c(summHSL,summSpatL,summMahalL,summMahal75L))
})

getMSE<-Vectorize(function(alpha,all=T){
  if(all){
    resultsHS_Sch<-getScwartzCP(alpha,resultsHS)
    resultsSpat_Sch<-getScwartzCP(alpha,resultsSpat)
    resultsMahal_Sch<-getScwartzCP(alpha,resultsMahal)
    resultsMahal75_Sch<-getScwartzCP(alpha,resultsMahal75)
    
    summHSL<-getSummaryL_MSE(resultsHS_Sch,N,theta,distr,d,numInt,thresh,depth="hs")
    summSpatL<-getSummaryL_MSE(resultsSpat_Sch,N,theta,distr,d,numInt,thresh,depth="spat")
    summMahalL<-getSummaryL_MSE(resultsMahal_Sch,N,theta,distr,d,numInt,thresh,depth="mahal")
    summMahal75L<-getSummaryL_MSE(resultsMahal75_Sch,N,theta,distr,d,numInt,thresh,depth="mahal75")
    
    return(c(summHSL,summSpatL,summMahalL,summMahal75L))
  }
  else{
    resultsMahal_Sch<-getScwartzCP(alpha,resultsMahal)
    summMahalL<-getSummaryL_MSE(resultsMahal_Sch,N,theta,distr,d,numInt,thresh,depth="mahal")
    return(summMahalL)
  }
})


#MSE PLOT
hs<-c()
sp<-c()
mh<-c()
mh75<-c()
mh_nc<-c()
##N=1000

for(i in 1:36){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*mod
  dName=names(distributions1)[params[3]]
  
  dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/WBS_SIM/Scen_1/"
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_1_simsize_",sim.size,sep="")
  fileName_s<-paste(dirr,fileName,".Rda",sep="")
  load(fileName_s)
  
  
  
  alphas=seq(.0001,2,length.out = 100)
  vals<-sqrt(getMSE(alphas))
  
  hs<-rbind(hs,vals[1,])
  sp<-rbind(sp,vals[2,])
  mh<-rbind(mh,vals[3,])
  mh75<-rbind(mh75,vals[4,])
  dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/WBS_SIM/Scen_2/"
  fileName_NC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_2_simsize_",sim.size,sep="")
  fileName_NC<-paste(dirr,fileName_NC,".Rda",sep="")
  load(fileName_NC)
  mh_nc=rbind( mh_nc,sqrt(getMSE(alphas,F)))
}

a=c(colors()[617],"black","violetred")
a2<-c(a)
cols=a2[paramterIndices[,3]]
ltys=paramterIndices[,2]

alphaM<-t(matrix(alphas,nrow=100,ncol=36))
matplot(t(alphaM),t(hs),type='l',lwd=2,col=cols,lty=ltys)
matplot(t(alphaM),t(sp),type='l',lwd=2,col=cols)
matplot(t(alphaM),t(mh),type='l',lwd=2,col=cols)
matplot(t(alphaM),t(mh75),type='l',lwd=2,col=cols)
matplot(t(alphaM),t(mh_nc),type='l',lwd=2,col=cols)

# 
# sub=paramterIndices[,3]==3
# a2<-c(a,2)
# cols=a[paramterIndices[,3]]
# ltys=paramterIndices[sub,2]
# 
# alphaM<-t(matrix(alphas,nrow=100,ncol=36))
# matplot(t(alphaM)[,sub],t(hs)[,sub],type='l',lwd=2,col=cols,lty=ltys)
# matplot(t(alphaM)[,sub],t(sp)[,sub],type='l',lwd=2,col=cols)
# matplot(t(alphaM)[,sub],t(mh)[,sub],type='l',lwd=2,col=cols)
# matplot(t(alphaM)[,sub],t(mh75)[,sub],type='l',lwd=2,col=cols)



pdf(file="HS_Alpha.pdf",height=8,width=10)
matplot(t(alphaM),t(hs),type='l',lwd=2,col=cols,lty=ltys,ylab="",cex.main=2,cex.axis=2,
        xlab="",cex=1.75,ylim=c(0,8),main="N=1000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
mtext(expression(alpha),1, line=3,adj=.5,cex=3)
# abline(h=0.5,lty=1,col="black")
dev.off()
pdf(file="SPAT_Alpha.pdf",height=8,width=10)
matplot(t(alphaM),t(sp),type='l',lwd=2,col=cols,lty=ltys,ylab="",cex.main=2,cex.axis=2,
        xlab="",cex=1.75,ylim=c(0,8),main="N=1000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
mtext(expression(alpha),1, line=3,adj=.5,cex=3)
dev.off()
pdf(file="MH_Alpha.pdf",height=8,width=10)
matplot(t(alphaM),t(mh),type='l',lwd=2,col=cols,lty=ltys,ylab="",cex.main=2,cex.axis=2,
        xlab="",cex=1.75,ylim=c(0,8),main="N=1000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
mtext(expression(alpha),1, line=3,adj=.5,cex=3)
dev.off()
pdf(file="MH75_Alpha.pdf",height=8,width=10)
matplot(t(alphaM),t(mh75),type='l',lwd=2,col=cols,lty=ltys,ylab="",cex.main=2,cex.axis=2,
        xlab="",cex=1.75,ylim=c(0,8),main="N=1000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
mtext(expression(alpha),1, line=3,adj=.5,cex=3)
dev.off()

#other graphs for supplementary
pdf(file="M_1000_NC_Alpha.pdf",height=8,width=10)
matplot(t(alphaM),t(mh_nc),type='l',lwd=2,col=cols,lty=ltys,ylab="",cex.main=2,cex.axis=2,
        xlab="",cex=1.75,ylim=c(0,8),main="N=1000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
mtext(expression(alpha),1, line=3,adj=.5,cex=3)
dev.off()




#MSE PLOT
mh_oc<-c()
mh_nc<-c()
##N=1000

for(i in 37:108){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*mod
  dName=names(distributions1)[params[3]]
  
  
  
  dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/WBS_SIM/Scen_1/"
  alphas=seq(.0001,2,length.out = 100)
  fileName_OC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_1_simsize_",sim.size,sep="")
  fileName_OC<-paste(dirr,fileName_OC,".Rda",sep="")
  load(fileName_OC)
  mh_oc=rbind( mh_oc,sqrt(getMSE(alphas,F)))
  
  dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/WBS_SIM/Scen_2/"
  fileName_NC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_2_simsize_",sim.size,sep="")
  fileName_NC<-paste(dirr,fileName_NC,".Rda",sep="")
  load(fileName_NC)
  mh_nc=rbind( mh_nc,sqrt(getMSE(alphas,F)))
}


alphaM<-t(matrix(alphas,nrow=100,ncol=36))
matplot(t(alphaM),t(mh_oc[1:36,]),type='l',lwd=2,col=cols)
matplot(t(alphaM),t(mh_nc[1:36,]),type='l',lwd=2,col=cols)


matplot(t(alphaM),t(mh_oc[37:72,]),type='l',lwd=2,col=cols)
matplot(t(alphaM),t(mh_nc[37:72,]),type='l',lwd=2,col=cols)



pdf(file="M_2500_OC_Alpha.pdf",height=8,width=10)
matplot(t(alphaM),t(mh_oc[1:36,]),type='l',lwd=2,col=cols,lty=ltys,ylab="",cex.main=2,cex.axis=2,
        xlab="",cex=1.75,ylim=c(0,8),main="N=2500")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
mtext(expression(alpha),1, line=3,adj=.5,cex=3)
dev.off()
pdf(file="M_2500_NC_Alpha.pdf",height=8,width=10)
matplot(t(alphaM),t(mh_nc[1:36,]),type='l',lwd=2,col=cols,lty=ltys,ylab="",cex.main=2,cex.axis=2,
        xlab="",cex=1.75,ylim=c(0,8),main="N=2500")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
mtext(expression(alpha),1, line=3,adj=.5,cex=3)
dev.off()
pdf(file="M_5000_OC_Alpha.pdf",height=8,width=10)
matplot(t(alphaM),t(mh_oc[37:72,]),type='l',lwd=2,col=cols,lty=ltys,ylab="",cex.main=2,cex.axis=2,
        xlab="",cex=1.75,ylim=c(0,8),main="N=5000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
mtext(expression(alpha),1, line=3,adj=.5,cex=3)
dev.off()
pdf(file="M_5000_NC_Alpha.pdf",height=8,width=10)
matplot(t(alphaM),t(mh_nc[37:72,]),type='l',lwd=2,col=cols,lty=ltys,ylab="",cex.main=2,cex.axis=2,
        xlab="",cex=1.75,ylim=c(0,8),main="N=5000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
mtext(expression(alpha),1, line=3,adj=.5,cex=3)
dev.off()

