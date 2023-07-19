

library(MASS)
library(RColorBrewer)
library(stringi)
library(xtable)

dirr_save="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/PELT Codes/Plots/"


dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/KW_PELT_SIMULATION/"
setwd(dirr)
Ns<-c(100,200,1000,2500)
sim.size=100
thetas<-list(c(.333,.666), 
             c(.25,.5,.75),
             c(1/6,2/6,3/6,4/6,5/6))
ds=c(2,3,5,10)

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
constants=unique(round(constants1,2))
# constants=c(constants1,constants2[1:2])
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
    fileName<-paste0("PELT",N,"_",length(theta),"_",dName,"_",d,"_new_constant_C1_",constant,"OC_PELT_ranks_simsize_",sim.size,sep="")
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


res<-paramterIndices[1:36,3]==3
res=1:36
# res=37:72
cols=paramterIndices[,3]


ltys=paramterIndices[,3]
ltys=ltys[res]

Xvals<-matrix(constants,nrow=length(constants),ncol=36)
a1<-c(colors()[617],"black","violetred")
a=a1[paramterIndices[,3]]




par(mfrow=c(2,2))

matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',lty=ltys,lwd=2,ylim=c(0,10),col=rep(1:4,9))
matplot(Xvals[,1:36],sqrt(sp[,res+36]),type='l',lty=ltys,lwd=2,ylim=c(0,10),col=rep(1:4,9))
matplot(Xvals[,1:36],sqrt(m[,res+36*2]),type='l',lty=ltys,lwd=2,ylim=c(0,10),col=rep(1:4,9))
matplot(Xvals[,1:36],sqrt(m75[,res+36*3]),type='l',lty=ltys,lwd=2,ylim=c(0,10),col=rep(1:4,9))

matplot(Xvals[,1:36],sqrt(sp[,res]),type='l',lty=ltys,lwd=2,ylim=c(0,10),col=rep(1:4,9))
matplot(Xvals[,1:36],sqrt(sp[,res+36]),type='l',lty=ltys,lwd=2,ylim=c(0,10),col=rep(1:4,9))
matplot(Xvals[,1:36],sqrt(sp[,res+36*2]),type='l',lty=ltys,lwd=2,ylim=c(0,10),col=rep(1:4,9))
matplot(Xvals[,1:36],sqrt(sp[,res+36*3]),type='l',lty=ltys,lwd=2,ylim=c(0,10),col=rep(1:4,9))


make_dim_plot=function(mat ,leg=F,main="N=100"){
  l1=rowMeans(mat[,seq(1,36,by=4)])
  l2=rowMeans(mat[,seq(1,36,by=4)+1])
  l3=rowMeans(mat[,seq(1,36,by=4)+2])
  l4=rowMeans(mat[,seq(1,36,by=4)+3])
  rmse=cbind(l1,l2,l3,l4)
  matplot(Xvals[,res],rmse,type='l',col=1:4,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,main=main,ylim=c(0,10))
  mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
  title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
  if(leg)
    legend("topleft",legend=c("d=2","d=3","d=5","d=10"),lty=1,col=1:4,cex=2.5)
  return(rmse)
}

make_dim_plot=function(mat ,leg=F,mainn="N=100"){
  # l1=rowMeans(mat[,seq(1,36,by=4)])
  # l2=rowMeans(mat[,seq(1,36,by=4)+1])
  # l3=rowMeans(mat[,seq(1,36,by=4)+2])
  # l4=rowMeans(mat[,seq(1,36,by=4)+3])
  # rmse=cbind(l1,l2,l3,l4)
  matplot(Xvals[,1:36],mat,type='l',col=1:4,
          lty=ltys,lwd=2,xlab="",ylab="",
          cex.main=2,cex.axis=2,main=mainn,ylim=c(0,10))
  mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
  title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
  if(leg){
    legend("topleft",legend=c("d=2","d=3","d=5","d=10"),lty=1,col=1:4,cex=2.5)}
  # return(rmse)
}




# matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=1000")
par(mfrow=c(1,1))
Cairo::CairoPDF(file=paste(dirr_save,"lambda_dimension_100.pdf",sep=""),height=8,width=10)
make_dim_plot(sqrt(sp[,res]),T)
dev.off()
Cairo::CairoPDF(file=paste(dirr_save,"lambda_dimension_200.pdf",sep=""),height=8,width=10)
make_dim_plot(sqrt(sp[,res+36]),mainn="N=200")
dev.off()
Cairo::CairoPDF(file=paste(dirr_save,"lambda_dimension_1000.pdf",sep=""),height=8,width=10)
make_dim_plot(sqrt(sp[,res+36*2]),F,"N=1000")
dev.off()
Cairo::CairoPDF(file=paste(dirr_save,"lambda_dimension_2500.pdf",sep=""),height=8,width=10)
make_dim_plot(sqrt(sp[,res+36*3]),F,"N=2500")
dev.off()



par(mfrow=c(2,2))

make_dim_plot(sqrt(sp[,res]))
make_dim_plot(sqrt(sp[,res+36]))
make_dim_plot(sqrt(sp[,res+36*2]))
make_dim_plot(sqrt(sp[,res+36*3]))


matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:36],sqrt(sp[,res+36]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:36],sqrt(m[,res+36*2]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:36],sqrt(m75[,res+36*3]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))

par(mfrow=c(1,1))

matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:36],sqrt(hs[,res+36]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:36],sqrt(hs[,res+36*2]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:36],sqrt(sp[,res+36]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:36],sqrt(m[,res+36]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))
matplot(Xvals[,1:36],sqrt(m75[,res+36]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,10))

par(mfrow=c(2,2))

matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:36],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:36],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:36],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2)


matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,xlim=c(0.1,1),ylim=c(0,20))
matplot(Xvals[,1:36],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlim=c(0.02,1))
matplot(Xvals[,1:36],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlim=c(0.02,1))
matplot(Xvals[,1:36],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlim=c(0.02,1))


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


# pdf(file="HS_C_v2.pdf",height=8,width=10)
# matplot(Xvals[,res],sqrt(hs[,res]),type='l',col=a,
#         lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,xlim=c(0.001,1),
#         ylim=c(0,5),ylab="",main="N=100")# abline(h=0.5,lty=1,col="black")
# # title(ylab="RMSE", line=1, cex.lab=2.5)
# mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
# title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
# legend("topleft",legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a1,cex=2.5)
# 
# legend("topright",legend=c(paste(stri_unescape_unicode('\\u2113'),"=2"),
#                            paste(stri_unescape_unicode('\\u2113'),"=3"),
#                            paste(stri_unescape_unicode('\\u2113'),"=5")),lty=c(1,2,3),cex=2.5)
# 
# dev.off()
# pdf(file="SPAT_C_v2.pdf",height=8,width=10)
# matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=1000")
# mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
# title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
# legend("topleft",legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a1,cex=2.5)
# 
# legend("topright",legend=c(paste(stri_unescape_unicode('\\u2113'),"=2"),
#                            paste(stri_unescape_unicode('\\u2113'),"=3"),
#                            paste(stri_unescape_unicode('\\u2113'),"=5")),lty=c(1,2,3),cex=2.5)
# dev.off()
# pdf(file="MH_C_v2.pdf",height=8,width=10)
# matplot(Xvals[,res],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=1000")
# mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
# title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
# dev.off()
# pdf(file="MH75_C_v2.pdf",height=8,width=10)
# matplot(Xvals[,res],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,xlab="",ylab="",cex.main=2,cex.axis=2,xlim=c(0.14,.3),ylim=c(0,2),main="N=1000")
# mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
# title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
# dev.off()


# constants=constants1

# constants<-c(seq(0.0025,0.1,l=10))
# constants=round(constants,2); constants
# constants[1]=0.005

constants<-c(seq(0.03,0.15,l=13))
constants=round(constants,2); constants
constants



hs<-NULL
sp<-NULL
m<-NULL
m75<-NULL

hs<-matrix(0,ncol=numUniqueRuns,nrow=length(constants))
sp<-hs
m=hs
m75=hs

for(i in 1:(36*4)){
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  dName=names(distributions1)[params[3]]
  for(j in 1:length(constants)){
    constant=constants[j]
    fileName<-paste0("PELT",N,"_",length(theta),"_",dName,"_",d,"_new_constant_2_C1_",constant,"OC_PELT_ranks_simsize_",sim.size,sep="")
    fileName1<-paste("CC_New_Tuning/",fileName,".Rda",sep="")
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


res<-paramterIndices[1:36,3]==3
res=1:36
# res=37:72
cols=paramterIndices[,3]


ltys=paramterIndices[,3]
ltys=ltys[res]

Xvals<-matrix(constants,nrow=length(constants),ncol=36)
a1<-c(colors()[617],"black","violetred")
a=a1[paramterIndices[,3]]

par(mfrow=c(2,2))

matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:36],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:36],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:36],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2)



matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))


matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res+36]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res+36*2]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res+36*3]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))


matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=rep(1:4,9),lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res+36]),type='l',col=rep(1:4,9),lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res+36*2]),type='l',col=rep(1:4,9),lty=ltys,lwd=2,ylim=c(0,5))



pdf(file=paste(dirr_save,"SPAT_logn_100.pdf",sep=""),height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,
        lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,
        ylim=c(0,5),ylab="",main="N=100")# abline(h=0.5,lty=1,col="black")
# title(ylab="RMSE", line=1, cex.lab=2.5)
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
legend("topleft",legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a1,cex=2.5)

legend("topright",legend=c(paste(stri_unescape_unicode('\\u2113'),"=2"),
                           paste(stri_unescape_unicode('\\u2113'),"=3"),
                           paste(stri_unescape_unicode('\\u2113'),"=5")),lty=c(1,2,3),cex=2.5)

dev.off()
pdf(file=paste(dirr_save,"SPAT_logn_200.pdf",sep=""),height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res+36]),type='l',col=a,
        lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,
        ylim=c(0,5),ylab="",main="N=200")# abline(h=0.5,lty=1,col="black")
# title(ylab="RMSE", line=1, cex.lab=2.5)
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file=paste(dirr_save,"SPAT_logn_200.pdf",sep=""),height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res+36]),type='l',col=a,
        lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,
        ylim=c(0,5),ylab="",main="N=200")# abline(h=0.5,lty=1,col="black")
# title(ylab="RMSE", line=1, cex.lab=2.5)
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file=paste(dirr_save,"SPAT_logn_1000.pdf",sep=""),height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res+36*2]),type='l',col=a,
        lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,
        ylim=c(0,5),ylab="",main="N=1000")# abline(h=0.5,lty=1,col="black")
# title(ylab="RMSE", line=1, cex.lab=2.5)
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file=paste(dirr_save,"SPAT_logn_2500.pdf",sep=""),height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res+36*3]),type='l',col=a,
        lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,
        ylim=c(0,5),ylab="",main="N=2500")# abline(h=0.5,lty=1,col="black")
# title(ylab="RMSE", line=1, cex.lab=2.5)
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()









constants<-c(seq(0.03,0.15,l=13))
constants=round(constants,2); constants
constants
constants<-c(seq(0.03,0.2,l=17))
constants=round(constants,2); constants
constants
constants<-c(seq(0.03,0.23,l=21))
constants=round(constants,2); constants
constants
# 
# hs<-NULL
# sp<-NULL
# m<-NULL
# m75<-NULL

hs<-matrix(0,ncol=numUniqueRuns,nrow=length(constants))
sp<-hs
m=hs
m75=hs

for(i in 1:(36*4)){
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  dName=names(distributions1)[params[3]]
  for(j in 1:length(constants)){
    constant=constants[j]
    fileName<-paste0("PELT",N,"_",length(theta),"_",dName,"_",d,"_new_constant_3_C1_",constant,"OC_PELT_ranks_simsize_",sim.size,sep="")
    fileName1<-paste("CC_New_Tuning_4/",fileName,".Rda",sep="")
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


res<-paramterIndices[1:36,3]==3
res=1:36
# res=37:72
cols=paramterIndices[,3]


ltys=paramterIndices[,3]
ltys=ltys[res]

Xvals<-matrix(constants,nrow=length(constants),ncol=36)
a1<-c(colors()[617],"black","violetred")
a=a1[paramterIndices[,3]]

par(mfrow=c(2,2))

matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:36],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:36],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2)
matplot(Xvals[,1:36],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2)



matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(sp[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(m[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(m75[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))


matplot(Xvals[,1:36],sqrt(hs[,res]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res+36]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res+36*2]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res+36*3]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))





getBestC(sqrt(hs[,res]),constants)
getBestC(sqrt(hs[,res+36]),constants)
getBestC(sqrt(hs[,res+36*2]),constants)
getBestC(sqrt(hs[,36*3]),constants)

par(mfrow=c(1,1))
matplot(Xvals[,1:36],sqrt(hs[,res+36*2]),type='l',col=a,lty=ltys,lwd=2,ylim=c(0,5))
legend("topleft",legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a1,cex=2.5)

legend("topright",legend=c(paste(stri_unescape_unicode('\\u2113'),"=2"),
                           paste(stri_unescape_unicode('\\u2113'),"=3"),
                           paste(stri_unescape_unicode('\\u2113'),"=5")),lty=c(1,2,3),cex=2.5)



#Normal
par(mfrow=c(3,1))
matplot(Xvals[,1:36],sqrt(hs[,res[c(1:4,13:16,25:30)]+36*2]),type='l',lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res[c(1:4,13:16,25:30)+4]+36*2]),type='l',lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res[c(1:4,13:16,25:30)+8]+36*2]),type='l',lty=ltys,lwd=2,ylim=c(0,5))


#Normal
par(mfrow=c(3,1))
matplot(Xvals[,1:36],sqrt(hs[,res[c(1:4,13:16,25:30)]+36*1]),type='l',lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res[c(1:4,13:16,25:30)+4]+36*1]),type='l',lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res[c(1:4,13:16,25:30)+8]+36*1]),type='l',lty=ltys,lwd=2,ylim=c(0,5))

#Normal
par(mfrow=c(3,1))
matplot(Xvals[,1:36],sqrt(hs[,res[c(1:4,13:16,25:30)]]),type='l',lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res[c(1:4,13:16,25:30)+4]]),type='l',lty=ltys,lwd=2,ylim=c(0,5))
matplot(Xvals[,1:36],sqrt(hs[,res[c(1:4,13:16,25:30)+8]]),type='l',lty=ltys,lwd=2,ylim=c(0,5))

pdf(file=paste(dirr_save,"SPAT_rlogn_100.pdf",sep=""),height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res]),type='l',col=a,
        lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,
        ylim=c(0,5),ylab="",main="N=100")# abline(h=0.5,lty=1,col="black")
# title(ylab="RMSE", line=1, cex.lab=2.5)
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
legend("topleft",legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a1,cex=2.5)

legend("topright",legend=c(paste(stri_unescape_unicode('\\u2113'),"=2"),
                           paste(stri_unescape_unicode('\\u2113'),"=3"),
                           paste(stri_unescape_unicode('\\u2113'),"=5")),lty=c(1,2,3),cex=2.5)

dev.off()
pdf(file=paste(dirr_save,"SPAT_rlogn_200.pdf",sep=""),height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res+36]),type='l',col=a,
        lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,
        ylim=c(0,5),ylab="",main="N=200")# abline(h=0.5,lty=1,col="black")
# title(ylab="RMSE", line=1, cex.lab=2.5)
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file=paste(dirr_save,"SPAT_rlogn_1000.pdf",sep=""),height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res+36*2]),type='l',col=a,
        lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,
        ylim=c(0,5),ylab="",main="N=1000")# abline(h=0.5,lty=1,col="black")
# title(ylab="RMSE", line=1, cex.lab=2.5)
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()
pdf(file=paste(dirr_save,"SPAT_rlogn_2500.pdf",sep=""),height=8,width=10)
matplot(Xvals[,res],sqrt(sp[,res+36*3]),type='l',col=a,
        lty=ltys,lwd=2,xlab="",cex.main=2,cex.axis=2,
        ylim=c(0,5),ylab="",main="N=2500")# abline(h=0.5,lty=1,col="black")
# title(ylab="RMSE", line=1, cex.lab=2.5)
mtext("RMSE", 2, line=2,adj=.65, cex=2.5)
title(xlab=expression('C'[1]), line=1.5,cex.lab=2.5)
dev.off()

