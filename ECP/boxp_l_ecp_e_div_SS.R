
library(MASS)
library(stringi)
library(xtable)
library(latex2exp)


dirr<-"C:/Users/12RAM/OneDrive/Documents/research/PhD Thesis/MKWC changepoint/ecp_sim"
# dirr<- "C:/Users/k3ramsay/OneDrive/Documents/research/PhD Thesis/Change Point R Codes/SimResults1/RData files/"
setwd(dirr)


Ns<-c(100,200)
# Ns<-1000
#Number of intervals, changed from 3 to 8
mod=100
numInts<-floor(log(Ns))*mod

#simulation size, num repetitions
sim.size=100



#numInts<-5000
#xhange below
#0,2,3.5 cps
thetas<-list(c(.333,.666), 
             c(.25,.5,.75),
             c(1/6,2/6,3/6,4/6,5/6))
# thetas<-list(c(.25,.5,.75),
# c(1/6,2/6,3/6,4/6,5/6))




##d= 2,3,5,10
ds=c(2,3,5,10)
# ds=c(2)



#95 of BB dont dchange for bonff
# thresh<-1.22

distributions=1:2
names(distributions)<-c("Normal", "Cauchy")
##Create Parameter Vector

numUniqueRuns<-length(Ns)*length(thetas)*length(ds)*length(distributions)

paramterIndices<-matrix(0,ncol=4,nrow=numUniqueRuns)
curr=0
for(i1 in 1:length(Ns)){
  for(i2 in 1:length(thetas)){
    for(i3 in 1:length(distributions)){
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


distributions<-list(normals,cauchys)
names(distributions)=c("Normal","Cauchy")



getSummaryL<-function(results,theta){
  
  
  #l-lhat
  
  l<-length(theta)
  lhat<-lapply(results,length)
  lhat<-unlist(lhat)
  llhat<-lhat-l
  
  return(llhat)
  
}

getDist<-function(){
  
  results=lapply(results,function(x){      x=x[-1]; x=x[-length(x)]      })
  summ<-getSummaryL(results,theta)
  return(summ)
  
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

output=NULL

for(i in 1:numUniqueRuns){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  dName=names(distributions1)[params[3]]
  
  

  fileName<-paste0("ECP_DIV",N,"_",length(theta),"_",dName,"_",d,"_ecp_div_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,"/CC/",fileName,".Rda",sep="")
  load(fileName1)
  
  vals<-getDist()
  mean(vals)
  output<-cbind(output,vals)
  
}



#order by dimension'
output2=output[,order(paramterIndices[,4])]
paramterIndices[order(paramterIndices[,4]),]


bor=paramterIndices[order(paramterIndices[1:24,4]),][,3]
out=paramterIndices[order(paramterIndices[1:24,4]),][,2]
dimm=paramterIndices[order(paramterIndices[1:24,4]),][,4]

a<-c(colors()[617],"black","violetred")
bor=a[bor]
par(mfrow=c(1,1))
# par(mfrow=c(2,2))
# par(mfrow=c(1,1))

par(mar=c(0,4,0,0)+5)#sets margins of plotting area

makePlot(output2[,1:24])
makePlot(output2[,1:24+24])
# makePlot(output2[,73:108])

# par(mfrow=c(2,2))

#order by dimension'
output2_1=output[,1:24]
output2_2=output[,1:24+24]
output2_1=output2_1[,order(paramterIndices[1:24,4])]
output2_2=output2_2[,order(paramterIndices[1:24,4])]
# output2_3=output2_3[,order(paramterIndices[1:24,4])]



bor=paramterIndices[order(paramterIndices[1:24,4]),][,3]
out=paramterIndices[order(paramterIndices[1:24,4]),][,2]
dimm=paramterIndices[order(paramterIndices[1:24,4]),][,4]

a<-c(colors()[617],"black","violetred")
bor=a[bor]
par(mfrow=c(1,1))
# par(mfrow=c(2,2))
# par(mfrow=c(1,1))

par(mar=c(0,4,0,0)+5)#sets margins of plotting area

makePlot(output2_1)
makePlot(output2_2)
# makePlot(output2_3)

# par(mfrow=c(2,2))

Cairo::CairoPDF(file="Plots/ecp_div_BP_100_SS.pdf",height=6.5,width=11)
makePlot(output2_1)
title(main="N=100",cex.main=2.5)
dev.off()


Cairo::CairoPDF(file="Plots/ecp_div_BP_200_SS.pdf",height=6.5,width=11)
makePlot(output2_2)
# title(main="N=200",cex.main=2.5)
dev.off()

# Cairo::CairoPDF(file="ecp_div_BP_5000.pdf",height=6.5,width=11)
# makePlot(output2_3)
# title(main="N=5000",cex.main=2.5)
# dev.off()





