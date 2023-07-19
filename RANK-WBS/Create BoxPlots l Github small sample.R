

library(MASS)
library(xtable)
library(stringi)
library(latex2exp)

dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS_SIM/Scen_1/"
setwd(dirr)



Ns<-c(100,200,500,1000,2500,5000)
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



distributions<-as.list(1:3)
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
m_nc=NULL
alpha=0.9

for(i in 1:numUniqueRuns){
  print(i)
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*mod
  dName=names(distributions1)[params[3]]
  
  dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS_SIM/Scen_1/"
  fileName_OC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_1_simsize_",sim.size,sep="")
  fileName_OC<-paste(dirr,fileName_OC,".Rda",sep="")
  load(fileName_OC)
  
  
  
  
  vals<-getDist(alpha)
  
  hs<-cbind(hs,vals[[1]])
  sp<-cbind(sp,vals[[2]])
  m<-cbind(m,vals[[3]])
  m75<-cbind(m75,vals[[4]])
  
  # Uncomment for scen 2
  # dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS_SIM/Scen_2/"
  # fileName_NC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_2_simsize_",sim.size,sep="")
  # fileName_NC<-paste(dirr,fileName_NC,".Rda",sep="")
  # load(fileName_NC)
  # m_nc<-cbind(m_nc,getDist(alpha,F))
  
}

v=73:108
v=37:72
v=1:36
hs2=hs[,order(paramterIndices[v,4])]
sp2=sp[,order(paramterIndices[v,4])]
m2=m[,order(paramterIndices[v,4])]
m752=m75[,order(paramterIndices[v,4])]

bor=paramterIndices[order(paramterIndices[v,4]),][,3]
out=paramterIndices[order(paramterIndices[v,4]),][,2]
dimm=paramterIndices[order(paramterIndices[v,4]),][,4]

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






Cairo::CairoPDF(file="HS_BP.pdf",height=6.5,width=11)
makePlot(hs2)
legend(0.325,6,legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a,cex=1.5)
dev.off()

Cairo::CairoPDF(file="SP_BP.pdf",height=6.5,width=11)
makePlot(sp2)
dev.off()
Cairo::CairoPDF(file="M_BP.pdf",height=6.5,width=11)
makePlot(m2)
dev.off()
Cairo::CairoPDF(file="M75_BP.pdf",height=6.5,width=11)
makePlot(m752)
dev.off()


m_oc<-NULL
# m_nc<-NULL
alpha=0.9

for(i in 37:numUniqueRuns){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*mod
  dName=names(distributions1)[params[3]]
  
  dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS_SIM/Scen_1/"
  fileName_OC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_1_simsize_",sim.size,sep="")
  fileName_OC<-paste(dirr,fileName_OC,".Rda",sep="")
  
  dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS_SIM/Scen_2/"
  fileName_NC<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_Scen_2_simsize_",sim.size,sep="")
  fileName_NC<-paste(dirr,fileName_NC,".Rda",sep="")
  
  load(fileName_OC)
  m_oc<-cbind(m_oc,getDist(alpha,F))
  
  load(fileName_NC)
  m_nc<-cbind(m_nc,getDist(alpha,F))
  
}


m_oc2=m_oc[,order(paramterIndices[37:108,4])]
m_nc2=m_nc[,order(paramterIndices[,4])]

bor=paramterIndices[order(paramterIndices[37:108,4]),][,3]
out=paramterIndices[order(paramterIndices[37:108,4]),][,2]
dimm=paramterIndices[order(paramterIndices[37:108,4]),][,4]

a<-c(colors()[617],"black","violetred")
bor=a[bor]

par(mfrow=c(1,1))
makePlot(m_nc2[,1:36])

makePlot(m_oc2[,1:36])
makePlot(m_nc2[,37:72])

makePlot(m_oc2[,37:72])
makePlot(m_nc2[,73:108])

Cairo::CairoPDF(file="M_l_BP_N_1000_NC_alpha_p9.pdf",height=6.5,width=11)
makePlot(m_nc2[,1:36])
title(main="N=1000, Scenario 2")
dev.off()

Cairo::CairoPDF(file="M_l_BP_N_2500_OC_alpha_p9.pdf",height=6.5,width=11)
makePlot(m_oc2[,1:36])
title(main="N=2500, Scenario 1")
dev.off()

Cairo::CairoPDF(file="M_l_BP_N_2500_NC_alpha_p9.pdf",height=6.5,width=11)
makePlot(m_nc2[,37:72])
title(main="N=2500, Scenario 2")
dev.off()
Cairo::CairoPDF(file="M_l_BP_N_5000_OC_alpha_p9.pdf",height=6.5,width=11)
makePlot(m_oc2[,37:72])
title(main="N=5000, Scenario 1")
dev.off()
Cairo::CairoPDF(file="M_l_BP_N_5000_NC_alpha_p9.pdf",height=6.5,width=11)
makePlot(m_nc2[,73:108])
title(main="N=5000, Scenario 2")
dev.off()


########SUBMATRIX STUFF ########SUBMATRIX STUFF ########SUBMATRIX STUFF ########SUBMATRIX STUFF ########SUBMATRIX STUFF
#sub matrix plots
########SUBMATRIX STUFF ########SUBMATRIX STUFF ########SUBMATRIX STUFF ########SUBMATRIX STUFF ########SUBMATRIX STUFF ########SUBMATRIX STUFF

#simulation parameters
#N,theta,rdata,d,numInt,thresh=1.358,data depth (MH,spatial, Tukey,MH75)

#N
Ns <- c(1000)

#simulation size, num repetitions
sim.size = 100




#xhange below
#0,2,3.5 cps
thetas <- list(c(.333, .666),
               c(.25, .5, .75),
               c(1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6))





##d= 2,3,5,10
ds = c(5)
thresh=1
d2s = c(4, 3, 2, 1)

distributions1 = 1
names(distributions1) <- c("Normal")
##Create Parameter Vector

numUniqueRuns <- length(thetas) * length(d2s)

paramterIndices <- matrix(0, ncol = 4, nrow = numUniqueRuns)
curr = 0
for (i1 in 1:length(Ns)) {
  for (i2 in 1:length(thetas)) {
    for (i3 in 1:length(d2s)) {
      for (i4 in 1:length(ds)) {
        curr = curr + 1
        paramterIndices[curr, ] = c(i1, i2, i3, i4)
      }
    }
  }
}



hs<-NULL
sp<-NULL
m<-NULL
m75<-NULL
alpha=0.9

for(i in 1:nrow(paramterIndices)){
  
  params <- paramterIndices[i, ]
  N = 1000
  theta = thetas[[params[2]]]
  # distr = distributions1[1]
  d2 = d2s[params[3]]
  d = 5
  numInt = floor(log(Ns[params[1]])) *100
  dName = "Normal"
  # distr = distributions[[1]]
  
  dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/WBS_SIM/Scen_3/"
  fileName <-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_", "thresh", thresh,"_d2_", d2, "_WBS_simsize_", sim.size, "_all", sep = "")
  fileName1 <- paste(dirr, fileName, ".Rda", sep = "")
  load(fileName1)
  
  
  
  vals<-getDist(alpha)
  
  hs<-cbind(hs,vals[[1]])
  sp<-cbind(sp,vals[[2]])
  m<-cbind(m,vals[[3]])
  m75<-cbind(m75,vals[[4]])
  
}


hs2=hs
sp2=sp
m2=m
m752=m75

bor=paramterIndices[order(paramterIndices[,4]),][,3]
out=paramterIndices[order(paramterIndices[,4]),][,2]
dimm=paramterIndices[order(paramterIndices[,4]),][,4]

a<-c(colors()[617],"black","violetred")
bor=a[bor]



makePlot2<-function(vals){
  
  boxplot(vals,ylim=c(-6,4),border=c(a,"gold"),xaxt="n",frame=F,cex.axis=1.25,lwd=2)
  legend("topright",legend=c("4","3","2","1"),
         cex=1.2,col=c(a,"gold"),lty=1,title=expression(italic(b)))
  txt <-rep(c("2","3","5"),each=1)
  mtext(paste(stri_unescape_unicode('\\u2113')," -",stri_unescape_unicode('\\u2113'),sep="") , 
        2, line=2,adj=.5, cex=2.53)
  mtext(expression(hat(" ")), 2, line=3.2,adj=.46, cex=2)
  axis(1, at=-.1, tick=F, labels=stri_unescape_unicode('\\u2113'), mgp=c(0,0,0),cex.axis=2)
  # axis(1, at=0.1, tick=F, labels=":", mgp=c(0,0,0),cex.axis=2)
  
  axis(1, at=c(1,2.5,4),  tcl=0, labels=c("",txt[1],""), mgp=c(0,0,0),cex.axis=1.2)
  
  for(i in 2:length(txt))
    axis(1, at=c(1,2.5,4)+(i-1)*4, tcl=0, labels=c("",txt[i],""), mgp=c(0,0,0),cex.axis=1.2)
  
  for(i in 1:length(txt))
    axis(1, at=c(1,2,3,4)+(i-1)*4, tcl=.2, labels=c("","","",""), mgp=c(0,0,0),cex.axis=1.2)
  
  
  mtext(stri_unescape_unicode('\\u2113'),1, line=2,adj=.5, cex=2.53)
  # #changes in dimension
  # axis(1, at=-0.1, tick=F,line = 2, labels=expression(italic("d")), mgp=c(0,0,0),cex.axis=2)
  # axis(1, at=0.1, line=2,tick=F, labels=":", mgp=c(0,0,0),cex.axis=2)
  # 
  # 
  # txt=c("2","3","5","10")
  # axis(1, at=c(1,5,9),  tcl=0, labels=c("",txt[1],""), mgp=c(0,0,0),line=2,cex.axis=1.25)
  # 
  # for(i in 2:length(txt))
  #   axis(1, at=c(1,5,9)+(i-1)*9, tcl=0, labels=c("",txt[i],""), mgp=c(0,0,0),line=2,cex.axis=1.25)
}

par(mar=c(0,2,0,2)+5)#sets margins of plotting area
par(mfrow=c(1,1))

makePlot2(hs2)
makePlot2(sp2)
makePlot2(m2)
makePlot2(m752)

# Cairo::CairoPDF("hs_wbs_sub.pdf",height=6.5,width=7)
Cairo::CairoPDF("hs_wbs_sub.pdf",height=6.5,width=7)
makePlot2(hs2)
dev.off()
Cairo::CairoPDF("sp_wbs_sub.pdf",height=6.5,width=7)
makePlot2(sp2)
dev.off()
Cairo::CairoPDF("m_wbs_sub.pdf",height=6.5,width=7)
makePlot2(m2)
dev.off()
Cairo::CairoPDF("m75_wbs_sub.pdf",height=6.5,width=7)
makePlot2(m752)
dev.off()

