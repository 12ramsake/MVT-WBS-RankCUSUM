

library(MASS)
library(RColorBrewer)
library(stringi)
library(xtable)


dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/KW_PELT_SIMULATION/"
setwd(dirr)
Ns<-c(100,200)
sim.size=100
thetas<-list(c(.333,.666), 
             c(.25,.5,.75),
             c(1/6,2/6,3/6,4/6,5/6))
ds=c(2,3,5,10)

distributions1=1:2
names(distributions1)<-c("Normal", "Cauchy")

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



distributions<-list(1,2,3)
names(distributions)=c("Normal","Cauchy","Skew Normal")

getBestC<-function(vals,constants){
  constants[which.min(apply(vals,1,sum))]
}


getSummaryL<-function(results,theta){
  
  
  #l-lhat
  results<-lapply(results,function(x){x[x!=0]})
  l<-length(theta)
  lhat<-lapply(results,length)
  lhat<-unlist(lhat)
  llhat<-lhat-l
  
  # return(median(llhat,trim=0))
  return(mean(llhat^2,trim=0))
  
}


getDist<-Vectorize(function(){
  
  
  
  summHSL<-getSummaryL(resultsHS,theta)
  summSpatL<-getSummaryL(resultsSpat,theta)
  summMahalL<-getSummaryL(resultsMH,theta)
  summMahal75L<-getSummaryL(resultsMH75,theta)
  
  return(c(summHSL,summSpatL,summMahalL,summMahal75L))
})



constants<-c(0.001,0.01,0.1,1,3,5)
constants2<-c(seq(0.05,1.1,l=20))
constants2=round(constants2,2)
constants1<-c(seq(0.005,0.05,l=10))
constants1=round(constants1,2)
constants=c(constants1,constants2[1:2])
# constants=constants1
hs<-NULL
sp<-NULL
m<-NULL
m75<-NULL

hs<-matrix(0,ncol=numUniqueRuns,nrow=length(constants))
sp<-hs
m=hs
m75=hs

for(i in 1:numUniqueRuns){
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  dName=names(distributions1)[params[3]]
  for(j in 1:length(constants)){
    constant=constants[j]
    fileName<-paste0("PELT",N,"_",length(theta),"_",dName,"_",d,"_new_constant_C1_",constant,"SS_PELT_ranks_simsize_",sim.size,sep="")
    fileName1<-paste("CC/",fileName,".Rda",sep="")
    e=load(fileName1)
    
    resultsMH=resultsMahal
    resultsMH75=resultsMahal75
    
    vals<-getDist()
    
    hs[j,i]<-vals[1]
    sp[j,i]<-vals[2]
    m[j,i]<-vals[3]
    m75[j,i]<-vals[4]
    
  }
}


res<-paramterIndices[1:24,3]==3
res=1:24
# res=37:72
cols=paramterIndices[,3]


ltys=paramterIndices[,3]
ltys=ltys[res]

Xvals<-matrix(constants,nrow=length(constants),ncol=24)
a1<-c(colors()[617],"black","violetred")
a=a1[paramterIndices[,3]]



par(mfrow=c(1,1))

matplot(Xvals[,1:24],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:24],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:24],sqrt(hs[,res+24]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:24],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:24],sqrt(sp[,res+24]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:24],sqrt(m[,res+24]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:24],sqrt(m75[,res+24]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))

par(mfrow=c(2,2))

matplot(Xvals[,1:24],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:24],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:24],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:24],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2)


matplot(Xvals[,1:24],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,xlim=c(0.1,1),ylim=c(0,20))
matplot(Xvals[,1:24],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlim=c(0.02,1))
matplot(Xvals[,1:24],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlim=c(0.02,1))
matplot(Xvals[,1:24],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlim=c(0.02,1))


bor=paramterIndices[,3]
bor=a[bor]
ltys=paramterIndices[,2]



par(mfrow=c(2,2))
matplot(Xvals[,res],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[,res],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[,res],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
par(mfrow=c(1,1))
matplot(Xvals[,res],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[,res],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[,res],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")

getBestC(sqrt(hs),constants)
getBestC(sqrt(sp),constants)
getBestC(sqrt(m),constants)
getBestC(sqrt(m75),constants)


pdf(file="HS_C_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(hs[,res]),type='l',col=a,
        lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,xlim=c(0.001,1),
        ylim=c(0,5),ylab="",main="N=100")# abline(h=0.5,lty=1,col="black")
# title(ylab="RMSE", line=1, cex.lab=2.5)
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
legend("topleft",legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a1,cex=2.5)

legend("topright",legend=c(paste(stri_unescape_unicode('\\u2113'),"=2"),
                           paste(stri_unescape_unicode('\\u2113'),"=3"),
                           paste(stri_unescape_unicode('\\u2113'),"=5")),lty=c(1,2,3),cex=2.5)

dev.off()
pdf(file="SPAT_C_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=1000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
legend("topleft",legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a1,cex=2.5)

legend("topright",legend=c(paste(stri_unescape_unicode('\\u2113'),"=2"),
                           paste(stri_unescape_unicode('\\u2113'),"=3"),
                           paste(stri_unescape_unicode('\\u2113'),"=5")),lty=c(1,2,3),cex=2.5)
dev.off()
pdf(file="MH_C_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=1000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file="MH75_C_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=1000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()



hs<-NULL
sp<-NULL
m<-NULL
m75<-NULL


hs<-matrix(0,ncol=24,nrow=length(constants))
sp<-hs
m=hs
m75=hs

for(i in 1:24){
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  dName=names(distributions1)[params[3]]
  for(j in 1:length(constants)){
    constant=constants[j]
    dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/KW_PELT_SIMULATION/Scen_2/"
    fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_constant_C1_",constant,"NC_PELT_ranks_simsize_",sim.size,sep="")
    fileName1<-paste(dirr,fileName,".Rda",sep="")
    load(fileName1)
    
    resultsMH=resultsMahal
    resultsMH75=resultsMahal75
    
    vals<-getDist()
    
    hs[j,i]<-vals[1]
    sp[j,i]<-vals[2]
    m[j,i]<-vals[3]
    m75[j,i]<-vals[4]
    
  }
}


#as d grows, the less sensitive it is to the parameter choice.

bor=paramterIndices[,3]
bor=a[bor]
ltys=paramterIndices[,2]


getBestC(sqrt(hs),constants  )
getBestC(sqrt(sp),constants)
getBestC(sqrt(m75),constants)
getBestC(sqrt(m),constants)




Xvals<-matrix(constants,nrow=length(constants),ncol=24)
par(mfrow=c(1,1))
matplot(Xvals[,res],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,
        xlab=expression('C'[1]),ylab="RMSE",cex.lab=1.75,xlim=c(0.14,.3),
        ylim=c(0,2),main="N=1000")
legend("topleft",legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a1,cex=1.5)
legend("topright",legend=c(paste(stri_unescape_unicode('\\u2113'),"=2"),
                           paste(stri_unescape_unicode('\\u2113'),"=3"),
                           paste(stri_unescape_unicode('\\u2113'),"=5")),lty=c(1,2,3),cex=1.5)
matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[,res],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[,res],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")


pdf(file="HS_C_NC_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(hs[,res]),type='l',col=a,
        lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),
        ylim=c(0,2),ylab="",main="N=1000")# abline(h=0.5,lty=1,col="black")
# title(ylab="RMSE", line=1, cex.lab=2.5)
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
legend("topleft",legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a1,cex=2.5)

legend("topright",legend=c(paste(stri_unescape_unicode('\\u2113'),"=2"),
                           paste(stri_unescape_unicode('\\u2113'),"=3"),
                           paste(stri_unescape_unicode('\\u2113'),"=5")),lty=c(1,2,3),cex=2.5)
dev.off()
pdf(file="SPAT_C_NC_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=1000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file="MH_C_NC_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=1000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file="MH75_C_NC_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=1000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()






hs<-matrix(0,ncol=24,nrow=length(constants))
sp<-hs
m=hs
m75=hs

for(i in 37:72){
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  dName=names(distributions1)[params[3]]
  for(j in 1:length(constants)){
    constant=constants[j]
    dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/KW_PELT_SIMULATION/Scen_1/"
    fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_constant_C1_",constant,"OC_PELT_ranks_simsize_",sim.size,sep="")
    fileName1<-paste(dirr,fileName,".Rda",sep="")
    load(fileName1)
    
    resultsMH=resultsMahal
    resultsMH75=resultsMahal75
    
    vals<-getDist()
    
    hs[j,i-24]<-vals[1]
    sp[j,i-24]<-vals[2]
    m[j,i-24]<-vals[3]
    m75[j,i-24]<-vals[4]
    
  }
}

res=1:24
#
cols=paramterIndices[,4]


#
ltys=paramterIndices[,3]
ltys=ltys[res]
Xvals<-matrix(constants,nrow=length(constants),ncol=24)

bor=paramterIndices[,3]
bor=a[bor]
ltys=paramterIndices[,2]


par(mfrow=c(2,2))
matplot(Xvals,hs,type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals,sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals,sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals,sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")


getBestC(hs,constants)
getBestC(sp,constants)
getBestC(m75,constants)
getBestC(m,constants)

pdf(file="HS_C_OC_N_2500_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=2500")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)



dev.off()
pdf(file="SPAT_C_OC_N_2500_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=2500")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file="MH_C_OC_N_2500_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=2500")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file="MH75_C_OC_N_2500_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=2500")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()



hs<-matrix(0,ncol=24,nrow=length(constants))
sp<-hs
m=hs
m75=hs

for(i in 73:108){
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  dName=names(distributions1)[params[3]]
  for(j in 1:length(constants)){
    constant=constants[j]
    dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/KW_PELT_SIMULATION/Scen_1/"
    fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_constant_C1_",constant,"OC_PELT_ranks_simsize_",sim.size,sep="")
    fileName1<-paste(dirr,fileName,".Rda",sep="")
    load(fileName1)
    
    resultsMH=resultsMahal
    resultsMH75=resultsMahal75
    
    vals<-getDist()
    
    hs[j,i-72]<-vals[1]
    sp[j,i-72]<-vals[2]
    m[j,i-72]<-vals[3]
    m75[j,i-72]<-vals[4]
    
  }
}
# #C=transformC(C,5000)
Xvals<-matrix(constants,nrow=length((constants)),ncol=24)
par(mfrow=c(1,1))
matplot(Xvals[1:length(constants),],hs,type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[1:length(constants),],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[1:length(constants),],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[1:length(constants),],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")


pdf(file="HS_C_OC_N_5000_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=5000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)


dev.off()

pdf(file="SPAT_C_OC_N_5000_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=5000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file="MH_C_OC_N_5000_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=5000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file="MH75_C_OC_N_5000_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=5000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()



getBestC(hs,constants)
getBestC(sp,constants)
getBestC(m75,constants)
getBestC(m,constants)



##########NC bigger N



hs<-matrix(0,ncol=24,nrow=length(constants))
sp<-hs
m=hs
m75=hs

for(i in 37:72){
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  dName=names(distributions1)[params[3]]
  for(j in 1:length(constants)){
    constant=constants[j]
    dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/KW_PELT_SIMULATION/Scen_2/"
    fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_constant_C1_",constant,"NC_PELT_ranks_simsize_",sim.size,sep="")
    fileName1<-paste(dirr,fileName,".Rda",sep="")
    load(fileName1)
    
    resultsMH=resultsMahal
    resultsMH75=resultsMahal75
    
    vals<-getDist()
    
    hs[j,i-24]<-vals[1]
    sp[j,i-24]<-vals[2]
    m[j,i-24]<-vals[3]
    m75[j,i-24]<-vals[4]
    
  }
}


res=1:24
#
cols=paramterIndices[,4]


#
ltys=paramterIndices[,3]
ltys=ltys[res]

Xvals<-matrix(constants,nrow=length(constants),ncol=24)

bor=paramterIndices[,3]
bor=a[bor]
ltys=paramterIndices[,2]





par(mfrow=c(1,1))
matplot(Xvals,hs,type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals,sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals,sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals,sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")


getBestC(hs,constants)
getBestC(sp,constants)
getBestC(m75,constants)
getBestC(m,constants)

pdf(file="HS_C_NC_N_2500_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=2500")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()

pdf(file="SPAT_C_NC_N_2500_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=2500")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()

pdf(file="MH_C_NC_N_2500_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=2500")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()

pdf(file="MH75_C_NC_N_2500_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=2500")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()




hs<-matrix(0,ncol=24,nrow=length(constants))
sp<-hs
m=hs
m75=hs

for(i in 73:108){
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  dName=names(distributions1)[params[3]]
  for(j in 1:length(constants)){
    constant=constants[j]
    dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/KW_PELT_SIMULATION/Scen_2/"
    fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_constant_C1_",constant,"NC_PELT_ranks_simsize_",sim.size,sep="")
    fileName1<-paste(dirr,fileName,".Rda",sep="")
    load(fileName1)
    
    resultsMH=resultsMahal
    resultsMH75=resultsMahal75
    
    vals<-getDist()
    
    hs[j,i-72]<-vals[1]
    sp[j,i-72]<-vals[2]
    m[j,i-72]<-vals[3]
    m75[j,i-72]<-vals[4]
    
  }
}


Xvals<-matrix(constants,nrow=length((constants)),ncol=24)
par(mfrow=c(1,1))
matplot(Xvals[1:length(constants),],hs,type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[1:length(constants),],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[1:length(constants),],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")
matplot(Xvals[1:length(constants),],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="C",ylab="RMSE")


pdf(file="HS_C_NC_N_5000_v2.pdf",height=8,width=10)
matplot(Xvals[1:length(constants),],hs,type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=5000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)

dev.off()
pdf(file="SPAT_C_NC_N_5000_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=5000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()

pdf(file="MH_C_NC_N_5000_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=5000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()

pdf(file="MH75_C_NC_N_5000_v2.pdf",height=8,width=10)
matplot(Xvals[,res],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",
        cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=5000")
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()

