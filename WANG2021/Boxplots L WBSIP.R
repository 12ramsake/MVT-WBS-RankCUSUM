
library(MASS)
library(stringi)
library(xtable)

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
distributions1=1:3
names(distributions1)<-c("Normal", "Cauchy", "Skew Normal")
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


distributions<-list(1:3)
names(distributions)=c("Normal","Cauchy","Skew Normal")



getSummaryL<-function(results,theta){
  
  
  #l-lhat
  
  l<-length(theta)
  lhat<-lapply(results,length)
  lhat<-unlist(lhat)
  llhat<-lhat-l
  
  return(llhat)
  
}

getDist<-function(thresh_bs=20,thresh_wbs=20){
  
  results_wbs=lapply(results, function(x){x$WBSIP[x$WBSIP[,2]>thresh_wbs,1]})
  results_bs=lapply(results, function(x){x$BSOP[x$BSOP[,2]>thresh_bs,1]})
  summ_wbs<-getSummaryL(results_wbs,theta)
  summ_bs<-getSummaryL(results_bs,theta)
  return(list(  wbs=summ_wbs,bs=summ_bs))
  
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



output_wbs=NULL
output_bs=NULL

for(i in 1:36){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*mod
  dName=names(distributions1)[params[3]]
  
  
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_WANG_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  load(fileName1)
  
  thresh_bs=37
  thresh_wbs=561
  

  
  vals<-getDist(thresh_bs,thresh_wbs)
  
  mean(vals$bs)
  mean(vals$wbs)
  
  # 
  # #minimize wrt mse
  # 
  # fnn<-Vectorize(function(x){median(abs(getDist(thresh_bs=thresh_bs,thresh_wbs=x)$wbs))})
 #  fnn<-Vectorize(function(x){mean(abs(getDist(thresh_bs=thresh_bs,thresh_wbs=x)$wbs)^2)})
 #  
 #  # curve(fnn(x),1,1000)
 #  # param_wbs=optimise(fnn,c(1,1000))$minimum
 #  # if(param_wbs>999)
 #  #   param_wbs=optimise(fnn,c(200,800))$minimum
 #  # print(param_wbs)
 # param_grid=seq(0,2000,by=0.1)
 #  # optimise(fnn,c(1,1000))
 # param_wbs=param_grid[which.min(fnn(param_grid))][1]
 #  # 
 #  # fnn2<-Vectorize(function(x){median(abs(getDist(thresh_bs=x,thresh_wbs=param_wbs)$bs))})
 #  fnn2<-Vectorize(function(x){mean(abs(getDist(thresh_bs=x,thresh_wbs=param_wbs)$bs)^2)})
 #  
 #  # plot(fnn2(1:100))
 #  # param_bs=optimise(fnn2,c(1,1000))$minimum
 #  param_bs=param_grid[which.min(fnn2(param_grid))][1]
 #  # 
 #  # 
 #  vals<-getDist( param_bs,  param_wbs)
 #  
  
  output_wbs<-cbind(output_wbs,vals$wbs)
  output_bs<-cbind(output_bs,vals$bs)
  
  
  
}




output_wbs2=output_wbs[,order(paramterIndices[1:36,4])]
output_bs2=output_bs[,order(paramterIndices[1:36,4])]

bor=paramterIndices[order(paramterIndices[1:36,4]),][,3]
out=paramterIndices[order(paramterIndices[1:36,4]),][,2]
dimm=paramterIndices[order(paramterIndices[1:36,4]),][,4]

a<-c(colors()[617],"black","violetred")
bor=a[bor]
par(mfrow=c(1,1))

par(mar=c(0,4,0,0)+5)#sets margins of plotting area



Cairo::CairoPDF(file="WBSIP_BP.pdf",height=6.5,width=11)
makePlot(output_wbs2)
# legend(0.325,6,legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a,cex=1.5)
dev.off()

makePlot(output_bs2)

Cairo::CairoPDF(file="BSOP_BP.pdf",height=6.5,width=11)
makePlot(output_bs2)
title(main="N=1000, BSOP",cex.main=2.5)
# legend(0.325,6,legend=c("Normal","Cauchy","Skew Normal"),lty=1,col=a,cex=1.5)
dev.off()




colMeans(abs(output_wbs2))
colMeans(abs(output_bs2))










