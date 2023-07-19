


library(MASS)
library(RColorBrewer)
library(stringi)
library(xtable)
library(latex2exp)


dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/KW_PELT_SIMULATION/"
setwd(dirr)
# constant=0.16
# constant=0.14
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

distributions<-as.list(1:2)
names(distributions)=c("Normal","Cauchy")



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
  
  
  
  summHSL<-getSummaryL(resultsHS,theta)
  summSpatL<-getSummaryL(resultsSpat,theta)
  summMahalL<-getSummaryL(resultsMH,theta)
  summMahal75L<-getSummaryL(resultsMH75,theta)
  # return(list(sp=summSpatL))
  return(list(hs=summHSL,sp=summSpatL,m=summMahalL,m75=summMahal75L))
}


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


hs<-NULL
sp<-NULL
m<-NULL
m75<-NULL
constant=0.18

for(i in 1:numUniqueRuns){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  dName=names(distributions1)[params[3]]
  dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/KW_PELT_SIMULATION/CC_SS_Old_C/"
  fileName<-paste0("PELT",N,"_",length(theta),"_",dName,"_",d,"_constant_C1_",constant,"SS_PELT_ranks_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  load(fileName1)
  
  resultsMH=resultsMahal
  resultsMH75=resultsMahal75
  
  vals<-getDist()
  
  hs<-cbind(hs,vals[[1]])
  sp<-cbind(sp,vals$sp)
  m<-cbind(m,vals[[3]])
  m75<-cbind(m75,vals[[4]])
  
  
}




vec=25:48
# vec=1:24
# vec=1:ncol(sp)
hs2=hs[,vec]
sp2=sp[,vec]
m2=m[,vec]
m752=m75[,vec]


hs2=hs2[,order(paramterIndices[vec,4])]
sp2=sp2[,order(paramterIndices[vec,4])]
m2=m2[,order(paramterIndices[vec,4])]
m752=m752[,order(paramterIndices[vec,4])]

bor=paramterIndices[order(paramterIndices[vec,4]),][,3]
out=paramterIndices[order(paramterIndices[vec,4]),][,2]
dimm=paramterIndices[order(paramterIndices[vec,4]),][,4]

a<-c(colors()[617],"black","violetred")
bor=a[bor]


makePlot(hs2)
makePlot(sp2)
makePlot(m2)
makePlot(m752)



par(mar=c(0,4,0,0)+5)#sets margins of plotting area

par(mfrow=c(1,1))
Cairo::CairoPDF(file="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/KW_PELT_SIMULATION/Plots/HS_BP_KW_SS.pdf",height=6.5,width=11)
makePlot(hs2)
legend("topright",legend=c("Normal","Cauchy"),lty=1,col=a,cex=1.5)
dev.off()

Cairo::CairoPDF(file="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/KW_PELT_SIMULATION/Plots/SP_BP_KW_SS.pdf",height=6.5,width=11)
makePlot(sp2)
dev.off()
Cairo::CairoPDF(file="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/KW_PELT_SIMULATION/Plots/M_BP_KW_SS.pdf",height=6.5,width=11)
makePlot(m2)
dev.off()
Cairo::CairoPDF(file="C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/KW_PELT_SIMULATION/Plots/M75_BP_KW_SS.pdf",height=6.5,width=11)
makePlot(m752)
dev.off()

#N=2500
# # 
# vec=25:48
# 
# 
# 
# hs2=hs
# sp2=sp
# m2=m
# m752=m75
# 
# 
# hs2=hs2[,order(paramterIndices[vec,4])]
# sp2=sp2[,order(paramterIndices[vec,4])]
# m2=m2[,order(paramterIndices[vec,4])]
# m752=m752[,order(paramterIndices[vec,4])]
# 
# bor=paramterIndices[order(paramterIndices[vec,4]),][,3]
# out=paramterIndices[order(paramterIndices[vec,4]),][,2]
# dimm=paramterIndices[order(paramterIndices[vec,4]),][,4]
# # a<-wes_palette(n=3, name="Cavalcanti1")
# a<-brewer.pal(n = 3, name = "Dark2")
# a<-c("#cc0066","#cc6600","#003300")
# a<-c(colors()[617],"black","violetred")
# bor=a[bor]
# 
# 
# 
# makePlot(hs2)
# makePlot(sp2)
# makePlot(m2)
# makePlot(m752)
# 
# 
# 
# 
# 
# par(mfrow=c(1,1))
# Cairo::CairoPDF(file="HS_BP_KW_N_2500_OC_SS.pdf",height=6.5,width=11)
# makePlot(hs2)
# title(main="N=2500, Scenario 1")
# dev.off()
# 
# Cairo::CairoPDF(file="SP_BP_KW_N_2500_OC_SS.pdf",height=6.5,width=11)
# makePlot(sp2)
# title(main="N=2500, Scenario 1")
# dev.off()
# Cairo::CairoPDF(file="M_BP_KW_N_2500_OC_SS.pdf",height=6.5,width=11)
# makePlot(m2)
# title(main="N=2500, Scenario 1")
# dev.off()
# Cairo::CairoPDF(file="M75_BP_KW_N_2500_OC_SS.pdf",height=6.5,width=11)
# makePlot(m752)
# title(main="N=2500, Scenario 1")
# dev.off()
# 
