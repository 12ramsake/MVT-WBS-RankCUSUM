

library(MASS)
library(xtable)
library(stringi)
library(latex2exp)

dirr<-""
setwd(dirr)



Ns<-c(1000,2500,5000)
mod=100
numInts<-floor(log(Ns))*mod
sim.size=100
thetas<-list(c(.333,.666), 
             c(.25,.5,.75),
             c(1/6,2/6,3/6,4/6,5/6))
ds=c(2,3,5,10)
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



distributions<-list(1:3)
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


getSummaryK<-function(results,N,theta,distt,d,numInt,thresh,depth){
  
  
  l<-length(theta)
  lhat<-lapply(results,length)
  lhat<-unlist(lhat)
  
  ##k-khat
  trueKs<-floor(N*theta)
  khat<-lapply(results,sort)
  
  #matrix of k-khat , one col for each cp
  khatk<-matrix(0,nrow=length(lhat),ncol=l)
  
  
  for(i in 1:length(results)){
    #all detected
    if(lhat[[i]]==l){
      khatk[i,]<-khat[[i]]-trueKs
    }
    else if(lhat[[i]]>l){
      spurious<-1:length(khat[[i]])
      for(j in 1:length(trueKs)){
        
        check<-abs(khat[[i]]-trueKs[j])
        ind<-which.min(check)
        
        while(sum(ind==spurious)==0){
          check[ind]=10000
          ind<-which.min(check)
          print("looping")
        }
        spurious<-spurious[spurious!=ind]
        khatk[i,j]<-khat[[i]][ind]-trueKs[j]
      }
    }
    else if(is.null(khat[[i]])){
      khat[[i]]<-rep(NA,l)
    }
    else{
      
      undett<-1:l
      for(j in 1:lhat[[i]]){
        
        
        check<-abs(khat[[i]][j]-trueKs)
        ind<-which.min(check)
        
        #if already detected
        while(sum(ind==undett)==0){
          check[ind]=10000
          ind<-which.min(check)
          print("looping2")
        }
        #remove deteced point
        undett<-undett[undett!=ind]
        khatk[i,ind]<-khat[[i]][j]-trueKs[ind]
        
      }
      #undetected
      
      if(!is.null(undett)){
        khatk[i,undett]<-rep(NA,length(undett))
      }
    }
  }
  
  
  # print(khatk/N)
  # 
  return(khatk/N)
}

getSummaryPUD<-function(results,N,theta,distt,d,numInt,thresh,depth){
  
  
  l<-length(theta)
  lhat<-lapply(results,length)
  lhat<-unlist(lhat)
  
  ##k-khat
  trueKs<-floor(N*theta)
  khat<-lapply(results,sort)
  
  #matrix of k-khat , one col for each cp
  khatk<-matrix(0,nrow=length(lhat),ncol=l)
  
  
  for(i in 1:length(results)){
    #all detected
    if(lhat[[i]]==l){
      khatk[i,]<-khat[[i]]-trueKs
    }
    else if(lhat[[i]]>l){
      spurious<-1:length(khat[[i]])
      for(j in 1:length(trueKs)){
        
        check<-abs(khat[[i]]-trueKs[j])
        ind<-which.min(check)
        
        while(sum(ind==spurious)==0){
          check[ind]=10000
          ind<-which.min(check)
          print("looping")
        }
        spurious<-spurious[spurious!=ind]
        khatk[i,j]<-khat[[i]][ind]-trueKs[j]
      }
    }
    else if(is.null(khat[[i]])){
      khat[[i]]<-rep(NA,l)
    }
    else{
      
      undett<-1:l
      for(j in 1:lhat[[i]]){
        
        
        check<-abs(khat[[i]][j]-trueKs)
        ind<-which.min(check)
        
        #if already detected
        while(sum(ind==undett)==0){
          check[ind]=10000
          ind<-which.min(check)
          print("looping2")
        }
        #remove deteced point
        undett<-undett[undett!=ind]
        khatk[i,ind]<-khat[[i]][j]-trueKs[ind]
        
      }
      #undetected
      
      if(!is.null(undett)){
        khatk[i,undett]<-rep(NA,length(undett))
      }
    }
  }
  
  pud<-apply(khatk,2,function(x){(is.na(x))})
  # print(khatk/N)
  # 
  return(pud)
}




getDist<-function(alpha,all=T){
  if(all){
    resultsHS_Sch<-getScwartzCP(alpha,resultsHS)
    resultsSpat_Sch<-getScwartzCP(alpha,resultsSpat)
    resultsMahal_Sch<-getScwartzCP(alpha,resultsMahal)
    resultsMahal75_Sch<-getScwartzCP(alpha,resultsMahal75)
    
    summHSL<-getSummaryK(resultsHS_Sch,N,theta,distr,d,numInt,thresh,depth="hs")
    summSpatL<-getSummaryK(resultsSpat_Sch,N,theta,distr,d,numInt,thresh,depth="spat")
    summMahalL<-getSummaryK(resultsMahal_Sch,N,theta,distr,d,numInt,thresh,depth="mahal")
    summMahal75L<-getSummaryK(resultsMahal75_Sch,N,theta,distr,d,numInt,thresh,depth="mahal75")
    
    return(list(hs=summHSL,dp=summSpatL,m=summMahalL,m75=summMahal75L))
  }
  else{
    resultsMahal_Sch<-getScwartzCP(alpha,resultsMahal)
    summMahalL<-getSummaryK(resultsMahal_Sch,N,theta,distr,d,numInt,thresh,depth="mahal")
    
    return(summMahalL)
    
  }
}

hs<-NULL
sp<-NULL
m<-NULL
m75<-NULL
cols=NULL

for(i in 1:36){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*mod
  dName=names(distributions1)[params[3]]
  if(d==10&&distr==1){
    fileName_OC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_1_simsize_",sim.size,sep="")
    fileName_OC<-paste(dirr,fileName_OC,'.Rda',sep="")

    
    
    
    alphas=0.9
    vals<-getDist(0.9)
    
    hs<-cbind(hs,vals[[1]])
    sp<-cbind(sp,vals[[2]])
    m<-cbind(m,vals[[3]])
    m75<-cbind(m75,vals[[4]])
    
    cols=c(cols,rep( params[3],dim(vals[[1]])[2]))
  }
}





# hs2=hs[,order(paramterIndices[1:36,4])]
# sp2=sp[,order(paramterIndices[1:36,4])]
# m2=m[,order(paramterIndices[1:36,4])]
# m752=m75[,order(paramterIndices[1:36,4])]

# cols<-rep(1,)
#1-24, 24-60, 61:end
vec=61:120
vec=1:24
vec=25:60
vec=1:ncol(hs)
hs2=hs[,vec]
sp2=sp[,vec]
m2=m[,vec]
m752=m75[,vec]
cols2<-cols[vec]
cols2="black"


par(mfrow=c(2,2))
boxplot(hs2,ylim=c(-.1,.1),border=cols2,xaxt="n",frame=F)
boxplot(sp2,ylim=c(-.1,.1),border=cols2,xaxt="n",frame=F)
boxplot(m2,ylim=c(-.1,.1),border=cols2,xaxt="n",frame=F)
boxplot(m752,ylim=c(-.1,.1),border=cols2,xaxt="n",frame=F)


par(mfrow=c(1,1),mgp=c(1.5,0,0),mar=c(0,2,0,0)+5)
makeBP<-function(vals){
  boxplot(vals,ylim=c(-.1,.1),border=1,xaxt="n",yaxt="n",frame=F,lwd=2,ylab="")
  mtext(TeX("$\\hat{k}/N-\\theta$"),2,cex=2,line=1)
  axis(2,line=-1.5,cex=1.75)
  txt <-rep(c("2","3","5"),length.out=9)
  mtext(stri_unescape_unicode('\\u2113'),1,line=2,cex=2.75)
  # axis(1.5, at=0.1, tick=F, labels=stri_unescape_unicode('\\u2113'), mgp=c(0,0,0),cex.axis=2)
  # axis(1.5, at=0.15, tick=F, labels=":", mgp=c(0,0,0),cex.axis=2)
  # 
  axis(1.5, at=c(1,1.5,2),  tcl=0, labels=c("","2",""), mgp=c(0,1,0),cex.axis=1.75)
  axis(1.5, at=c(3,4,5), tcl=0, labels=c("","3",""), mgp=c(0,1,0),cex.axis=1.75)
  axis(1.5, at=c(6,8.5,12), tcl=0, labels=c("","5",""), mgp=c(0,1,0),cex.axis=1.75)
  
}

Cairo::CairoPDF(file="kbp_HS.pdf",height=6,width=10)
makeBP(hs2)
dev.off()
Cairo::CairoPDF(file="kbp_SP.pdf",height=6,width=10)
makeBP(sp2)
dev.off()
Cairo::CairoPDF(file="kbp_MH.pdf",height=6,width=10)
makeBP(m2)
dev.off()
Cairo::CairoPDF(file="kbp_MH75.pdf",height=6,width=10)
makeBP(m752)
dev.off()
#10x7.5



##supplementary boxplots
makeBigPlot<-function(vals){
  
  boxplot(vals,ylim=c(-1,1)*0.15,border=bor,xaxt="n",frame=F,cex.axis=1.25,lwd=2)
  txt <-c("2","3","5")
  
  axis(1, at=-.1, tick=F, labels=stri_unescape_unicode('\\u2113'), mgp=c(0,0,0),cex.axis=2)
  axis(1, at=0.1, tick=F, labels=":", mgp=c(0,0,0),cex.axis=2)
  
  axis(1, at=c(1,3.5,6),  tcl=0, labels=c("",txt[1],""), mgp=c(0,0,0),cex.axis=1)
  axis(1, at=c(7,11,15),  tcl=0, labels=c("",txt[2],""), mgp=c(0,0,0),cex.axis=1)
  axis(1, at=c(16,23,30),  tcl=0, labels=c("",txt[3],""), mgp=c(0,0,0),cex.axis=1)
  
  
  
  for(i in 2:4){
    axis(1, at=c(1,3.5,6)+(i-1)*30,  tcl=0, labels=c("",txt[1],""), mgp=c(0,0,0),cex.axis=1)
    axis(1, at=c(7,11,15)+(i-1)*30,  tcl=0, labels=c("",txt[2],""), mgp=c(0,0,0),cex.axis=1)
    axis(1, at=c(16,23,30)+(i-1)*30,  tcl=0, labels=c("",txt[3],""), mgp=c(0,0,0),cex.axis=1)
  }
  
  
  # for(i in 1:length(txt))
  #   axis(1, at=c(1,2,3)+(i-1)*3, tcl=.2, labels=c("","",""), mgp=c(0,0,0),cex.axis=1.2)
  # 
  # #d
  axis(1, at=-0.1, tick=F,line = 2, labels=expression(italic("d")), mgp=c(0,0,0),cex.axis=2)
  axis(1, at=0.1, line=2,tick=F, labels=":", mgp=c(0,0,0),cex.axis=2)
  
  
  txt=c("2","3","5","10")
  axis(1, at=c(1,15,30),  tcl=0, labels=c("",txt[1],""), mgp=c(0,0,0),line=2,cex.axis=1.25)
  
  for(i in 2:length(txt))
    axis(1, at=c(1,15,30)+(i-1)*30, tcl=0, labels=c("",txt[i],""), mgp=c(0,0,0),line=2,cex.axis=1.25)
}

hs<-NULL
sp<-NULL
m<-NULL
m75<-NULL
cols=NULL
m_nc=NULL
dimss=NULL
for(i in 1:36){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*100
  dName=names(distributions1)[params[3]]
  

  fileName_OC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_1_simsize_",sim.size,".Rda",sep="")
  fileName_OC<-paste(dirr,fileName_OC,sep="")
  load(fileName_OC)
  
  
  
  alpha=0.9
  vals<-getDist(0.9)
  
  hs<-cbind(hs,vals[[1]])
  sp<-cbind(sp,vals[[2]])
  m<-cbind(m,vals[[3]])
  m75<-cbind(m75,vals[[4]])
  
  cols=c(cols,rep( params[3],dim(vals[[1]])[2]))
  dimss=c(dimss,rep( params[4],dim(vals[[1]])[2]))
  
  fileName_NC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_2_simsize_",sim.size,".Rda",sep="")
  fileName_NC<-paste(dirr,fileName_NC,sep="")
  load(fileName_NC)
  m_nc<-cbind(m_nc,getDist(alpha,F))
}
# }
hs2=hs[,order(dimss)]
sp2=sp[,order(dimss)]
m2=m[,order(dimss)]
m752=m75[,order(dimss)]
m_nc2=m_nc[,order(dimss)]
bor=paramterIndices[order(paramterIndices[1:36,4]),][,3]
# out=paramterIndices[order(paramterIndices[1:36,4]),][,2]
# dimm=paramterIndices[order(paramterIndices[1:36,4]),][,4]

a<-c(colors()[617],"black","violetred")
bor=a[bor]

makeBigPlot(hs2)
makeBigPlot(sp2)
makeBigPlot(m2)
makeBigPlot(m752)
makeBigPlot(m_nc2)

Cairo::CairoPDF(file="HS_k_BP_N_1000_OC_alpha_p9.pdf",height=6.5,width=10)
makeBigPlot(hs2)
title(main="N=1000, Scenario 1")
dev.off()
Cairo::CairoPDF(file="SP_k_BP_N_1000_OC_alpha_p9.pdf",height=6.5,width=10)
makeBigPlot(sp2)
title(main="N=1000, Scenario 1")
dev.off()
Cairo::CairoPDF(file="M_k_BP_N_1000_OC_alpha_p9.pdf",height=6.5,width=10)
makeBigPlot(m2)
title(main="N=1000, Scenario 1")
dev.off()
Cairo::CairoPDF(file="M_k_BP_N_1000_NC_alpha_p9.pdf",height=6.5,width=10)
makeBigPlot(m_nc2)
title(main="N=1000, Scenario 2")
dev.off()
Cairo::CairoPDF(file="M75_k_BP_N_1000_OC_alpha_p9.pdf",height=6.5,width=10)
makeBigPlot(m752)
title(main="N=1000, Scenario 1")
dev.off()




m_oc<-NULL
m_nc<-NULL
alpha=0.9

for(i in 37:numUniqueRuns){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*mod
  dName=names(distributions1)[params[3]]
  
  fileName_NC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_2_simsize_",sim.size,sep="")
  fileName_OC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_1_simsize_",sim.size,sep="")
  fileName_OC<-paste(dirr,fileName_OC,'.Rda',sep="")
  fileName_NC<-paste(dirr,fileName_NC,'.Rda',sep="")
  
  load(fileName_OC)
  # vals<-getDist(alpha)
  m_oc<-cbind(m_oc,getDist(alpha,F))
  
  load(fileName_NC)
  # vals<-getDist(alpha)
  m_nc<-cbind(m_nc,getDist(alpha,F))
  
}


m_oc2=m_oc[,1:120]
m_oc2=m_oc2[,order(dimss)]
m_nc2=m_nc[,1:120]
m_nc2=m_nc2[,order(dimss)]

bor=paramterIndices[order(paramterIndices[37:108,4]),][,3]
out=paramterIndices[order(paramterIndices[37:108,4]),][,2]
dimm=paramterIndices[order(paramterIndices[37:108,4]),][,4]

a<-c(colors()[617],"black","violetred")
bor=a[bor]

par(mfrow=c(1,1))
makeBigPlot(m_oc2)
makeBigPlot(m_nc2)

# makeBigPlot(m_oc2[,37:72])
# makeBigPlot(m_nc2[,37:72])



Cairo::CairoPDF(file="M_k_BP_N_2500_OC_alpha_p9.pdf",height=6.5,width=10)
makeBigPlot(m_oc2)
title(main="N=2500, Scenario 1")
dev.off()

Cairo::CairoPDF(file="M_k_BP_N_2500_NC_alpha_p9.pdf",height=6.5,width=10)
makeBigPlot(m_nc2)
title(main="N=2500, Scenario 2")
dev.off()











m_oc2=m_oc[,121:240]
m_oc2=m_oc2[,order(dimss)]
m_nc2=m_nc[,121:240]
m_nc2=m_nc2[,order(dimss)]





Cairo::CairoPDF(file="M_k_BP_N_5000_OC_alpha_p9.pdf",height=6.5,width=10)
makeBigPlot(m_oc2)
title(main="N=5000, Scenario 1")
dev.off()
Cairo::CairoPDF(file="M_k_BP_N_5000_NC_alpha_p9.pdf",height=6.5,width=10)
makeBigPlot(m_nc2)
title(main="N=5000, Scenario 2")
dev.off()


