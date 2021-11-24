
library(MASS)
library(RColorBrewer)
library(stringi)
library(xtable)
dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/KW_PELT_SIMULATION/"
setwd(dirr)



Ns<-c(1000)
sim.size=100
thetas<-list(.5)
ds=c(5,10,50,100,250,500)
distributions1=1
names(distributions1)<-c("Normal")
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

dim_of_change=5

getSummaryL<-function(results,theta){
  
  
  #l-lhat
  results<-lapply(results,function(x){x[x!=0]})
  l<-length(theta)
  lhat<-lapply(results,length)
  lhat<-unlist(lhat)
  llhat<-lhat-l
  
  return(llhat)
  
}


getDist<-function(){
  
  
  summSpatL<-getSummaryL(resultsSpat,theta)
  
  return(summSpatL)
}


makePlot<-function(vals){
  
  boxplot(vals,ylim=c(-6,6),border=bor,xaxt="n",yaxt="n",frame=F,cex.axis=1.25,lwd=2)
  
  mtext(paste(stri_unescape_unicode('\\u2113')," -",stri_unescape_unicode('\\u2113'),sep="") , 
        2, line=1.75, cex=2.6)
  axis(2, at=seq(-6,6,2),labels=seq(-6,6,2),line=-1)
  
  mtext(expression(hat(" ")), 2, line=3.3,adj=.46, cex=1.5)
  txt <-rep(c("2","3","5"),length.out=12)
  
  axis(1, at=-.1, tick=F, labels=stri_unescape_unicode('\\u2113'), mgp=c(0,0,0),cex.axis=2.5)
  # axis(1, at=0.1, tick=F, labels=":", mgp=c(0,0,0),cex.axis=2)
  
  axis(1, at=c(1,2,3),  tcl=0, labels=c("",txt[1],""), mgp=c(0,0,0),cex.axis=1.2)
  
  for(i in 2:length(txt))
    axis(1, at=c(1,2,3)+(i-1)*3, tcl=0, labels=c("",txt[i],""), mgp=c(0,0,0),cex.axis=1.2)
  
  for(i in 1:length(txt))
    axis(1, at=c(1,2,3)+(i-1)*3, tcl=.2, labels=c("","",""), mgp=c(0,0,0),cex.axis=1.2)
  
  #d
  # axis(1, at=-0.1, tick=F,line = 2, labels=expression(italic("d")), mgp=c(0,0,0),cex.axis=2)
  mtext(expression(italic("d")),1,line=3,cex=2.5)
  # axis(1, at=0.1, line=2,tick=F, labels=":", mgp=c(0,0,0),cex.axis=2)
  
  
  txt=c("2","3","5","10")
  axis(1, at=c(1,5,9),  tcl=0, labels=c("",txt[1],""), mgp=c(0,0,0),line=2,cex.axis=1.25)
  
  for(i in 2:length(txt))
    axis(1, at=c(1,5,9)+(i-1)*9, tcl=0, labels=c("",txt[i],""), mgp=c(0,0,0),line=2,cex.axis=1.25)
}


sp<-NULL
constant=0.18
dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/KW_PELT_SIMULATION/Scen_4/"

for(i in 1:numUniqueRuns){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  
  beta=constant*sqrt(N)+3.74
  
  
  dName=names(distributions1)[params[3]]
  
  fileName<-paste0(N,"_",length(theta),"_Normal_d_",d,"_HD3_dimchange_",dim_of_change,"_constant_",constant,"_PELT_ranks_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  load(fileName1)
  
  vals<-getDist()
  
  sp<-cbind(sp,vals)
  
  
}


sp
rmse=sqrt(apply(sp, 2, function(x){mean((x^2))}))
rmse
boxplot(sp)

signal=ds
plot(signal,rmse,type="l")


Cairo::CairoPDF(file="L_RMSE_HD.pdf",height=6.5,width=11)
plot(signal,rmse,type="l" ,frame=F,cex.axis=1.25,cex.lab=1.25,lwd=3,ylab="RMSE",xlab="Dimension")
dev.off()






getSummaryK<-function(results,theta){
  
  
  #l-lhat
  results=resultsSpat
  
  results<-lapply(results,function(x){x[x!=0]})
  dists=lapply(results, function(x){
    if(length(x)==0)
      return(NA)
    else if(length(x)==1){
      return((x-499))
    }
    else if(length(x)>1){
      return(x[which.min(abs(x-499))]-499)
    }
  })
  
  return(dists)
  
}


getDistK<-function(){
  
  
  summSpatK<-getSummaryK(resultsSpat,theta)
  
  return(summSpatK)
}








sp<-NULL

constant=0.18
dim_of_change=c(3,5,10,50,100,250,500)

dim_of_change=5

for(i in 1:numUniqueRuns){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  
  beta=constant*sqrt(N)+3.74
  
  
  dName=names(distributions1)[params[3]]
  
  fileName<-paste0(N,"_",length(theta),"_Normal_d_",d,"_HD3_dimchange_",dim_of_change,"_constant_",constant,"_PELT_ranks_simsize_",sim.size,sep="")  
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  load(fileName1)
  
  vals<-getDistK()
  
  sp<-cbind(sp,vals)
  
  
}



# as.matrix(sp)^2
mse=NULL
for(l in 1:6){
  mse=c(mse,mean(unlist(lapply(sp[,l],function(x){(x/1000)^2})),na.rm=T ))
}



plot(signal,mse,type="l")

plot(signal,sqrt(mse),type="l")


Cairo::CairoPDF(file="K_RMSE_HD.pdf",height=6.5,width=11)
plot(signal,sqrt(mse),type="l" ,frame=F,cex.axis=1.25,cex.lab=1.25,lwd=3,ylab="RMSE",xlab="Dimension")
dev.off()




################################################################################ increasing d and change 

thetas<-list(.5,c(.333,.666))

Ns=1000
##d= 2,3,5,10
ds=c(50,500)


distributions1=1
names(distributions1)<-c("Normal")
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
distributions<-list(1)
names(distributions)=c("Normal")

constant=0.18
sp<-NULL
for(i in 1:numUniqueRuns){
  
  params<-paramterIndices[i,]
  N=1000
  theta=thetas[[params[2]]]
  distr=1
  d=ds[params[4]]
  N=1000
  dName="Normal"
  
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_constant_",constant,"HD_PELT_ranks_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  load(fileName1)
  
  vals<-getDist()
  
  sp<-cbind(vals)
  
  
}


boxplot(sp)





