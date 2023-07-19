

library(MASS)
library(xtable)
library(stringi)
library(latex2exp)

dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS_SIM/CC/"
setwd(dirr)



Ns<-c(100,200)
mod=100
numInts<-floor(log(Ns))*mod
sim.size=100
thetas<-list(c(.333,.666), 
             c(.25,.5,.75),
             c(1/6,2/6,3/6,4/6,5/6))
ds=c(2,3,5,10)
thresh<-1

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
#simulation functions



distributions<-as.list(1:2)
names(distributions)=c("Normal","Cauchy")

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
    if(length(models)>1){
      for(j in 2:length(models)){
        
        sSic<-c(sSic,getsSic(models[[j]]$cp,models[[j]]$sigSq))
        
      }
    }
    minVal<-which.min(sSic)
    
    return(models[[minVal]]$cp)  
    
  }
  
  
  return(lapply(results,applySCH))
}


getSummaryL<-function(results,N,theta,distt,d,numInt,thresh,depth){
  
  
  #l-lhat
  
  l<-length(theta)
  lhat<-lapply(results,length)
  lhat<-unlist(lhat)
  llhat<-lhat-l
  
  return(llhat)
  
}
library(latex2exp)

makePlot<-function(vals){
  
  boxplot(vals,ylim=c(-6,6),border=bor,xaxt="n",yaxt="n",frame=F,cex.axis=1.25,lwd=2)
  
  mtext(paste(stri_unescape_unicode('\\u2113')," -",stri_unescape_unicode('\\u2113'),sep="") , 
        2, line=1.75, cex=2.6)
  axis(2, at=seq(-6,6,2),labels=seq(-6,6,2),line=-1)
  
  mtext(expression(hat(" ")), 2, line=3.3,adj=.46, cex=1.5)
  txt <-rep(c("2     ","3    ","5    "),length.out=12)
  
  axis(1, at=-.1, tick=F, labels=stri_unescape_unicode('\\u2113'), mgp=c(0,0,0),cex.axis=2.5)
  # axis(1, at=0.1, tick=F, labels=":", mgp=c(0,0,0),cex.axis=2)
  
  axis(1, at=c(1,2),  tcl=0, labels=c("",txt[1]), mgp=c(0,0,0),cex.axis=1.2)
  
  for(i in 2:length(txt))
    axis(1, at=c(1,2)+(i-1)*2, tcl=0, labels=c("",txt[i]), mgp=c(0,0,0),cex.axis=1.2)
  
  for(i in 1:length(txt))
    axis(1, at=c(1,2)+(i-1)*2, tcl=.2, labels=c("",""), mgp=c(0,0,0),cex.axis=1.2)
  
  #d
  # axis(1, at=-0.1, tick=F,line = 2, labels=expression(italic("d")), mgp=c(0,0,0),cex.axis=2)
  mtext(expression(italic("d")),1,line=3,cex=2.5)
  # axis(1, at=0.1, line=2,tick=F, labels=":", mgp=c(0,0,0),cex.axis=2)
  
  
  txt=c("2","3","5","10")
  axis(1, at=c(1,3.5,6),  tcl=0, labels=c("",txt[1],""), mgp=c(0,0,0),line=2,cex.axis=1.25)
  for(i in 2:length(txt))
    axis(1, at=c(1,3.5,6)+(i-1)*6, tcl=0, labels=c("",txt[i],""), mgp=c(0,0,0),line=2,cex.axis=1.25)
}

getDist<-function(alpha,all=T){
  
  if(all){
    resultsHS_Sch<-getScwartzCP(alpha,resultsHS)
    resultsSpat_Sch<-getScwartzCP(alpha,resultsSpat)
    resultsMahal_Sch<-getScwartzCP(alpha,resultsMahal)
    resultsMahal75_Sch<-getScwartzCP(alpha,resultsMahal75)
    
    summHSL<-getSummaryL(resultsHS_Sch,N,theta,distr,d,numInt,thresh,depth="hs")
    summSpatL<-getSummaryL(resultsSpat_Sch,N,theta,distr,d,numInt,thresh,depth="spat")
    summMahalL<-getSummaryL(resultsMahal_Sch,N,theta,distr,d,numInt,thresh,depth="mahal")
    summMahal75L<-getSummaryL(resultsMahal75_Sch,N,theta,distr,d,numInt,thresh,depth="mahal75")
    
    return(list(hs=summHSL,dp=summSpatL,m=summMahalL,m75=summMahal75L))
  }
  else{
    resultsMahal_Sch<-getScwartzCP(alpha,resultsMahal)
    summMahalL<-getSummaryL(resultsMahal_Sch,N,theta,distr,d,numInt,thresh,depth="mahal")
    return(summMahalL)
    
  }
}


hs<-NULL
sp<-NULL
m<-NULL
m75<-NULL
alpha=0.9

for(i in 1:numUniqueRuns){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*mod
  dName=names(distributions1)[params[3]]
  
  dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS_SIM/CC/WBS"
  fileName_OC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_SS_simsize_",sim.size,sep="")
  fileName_OC<-paste(dirr,fileName_OC,".Rda",sep="")
  load(fileName_OC)
  
  
  vals<-getDist(alpha)
  
  hs<-cbind(hs,vals[[1]])
  sp<-cbind(sp,vals[[2]])
  m<-cbind(m,vals[[3]])
  m75<-cbind(m75,vals[[4]])

}


vec=1:24
vec=25:48
hs2=hs[,order(paramterIndices[vec,4])]
sp2=sp[,order(paramterIndices[vec,4])]
m2=m[,order(paramterIndices[vec,4])]
m752=m75[,order(paramterIndices[vec,4])]

bor=paramterIndices[order(paramterIndices[vec,4]),][,3]
out=paramterIndices[order(paramterIndices[vec,4]),][,2]
dimm=paramterIndices[order(paramterIndices[vec,4]),][,4]


a<-c(colors()[617],"black","violetred")
bor=a[bor]

par(mfrow=c(1,1))
par(mfrow=c(2,2))
par(mfrow=c(1,1))

par(mar=c(0,4,0,0)+5)#sets margins of plotting area

# par(mfrow=c(2,2))
makePlot(sp2)
makePlot(hs2)
makePlot(m2)
makePlot(m752)



Cairo::CairoPDF(file="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS Codes/Plots/HS_BP_SS.pdf",height=6.5,width=11)
makePlot(hs2)
legend(0.325,6,legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a,cex=1.5)
dev.off()

Cairo::CairoPDF(file="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS Codes/Plots/SP_BP_SS.pdf",height=6.5,width=11)
makePlot(sp2)
dev.off()
Cairo::CairoPDF(file="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS Codes/Plots/M_BP_SS.pdf",height=6.5,width=11)
makePlot(m2)
dev.off()
Cairo::CairoPDF(file="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS Codes/Plots/M75_BP_SS.pdf",height=6.5,width=11)
makePlot(m752)
dev.off()

