
library(MASS)
library(stringi)
library(xtable)
dirr<-""
setwd(dirr)


Ns<-c(1000,2500,5000)
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







distributions<-list(1:3)
names(distributions)=c("Normal","Cauchy","Skew Normal")






getSummaryK<-function(results,N,theta){
  
  
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
    else if(length(khat[[i]])==0){
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


getDist<-function(){
  
  summ<-getSummaryK(results,N,theta)
  return(summ)
  
}


output=NULL
cols=NULL

for(i in 1:36){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*mod
  dName=names(distributions1)[params[3]]
  if(d==3&&distr==2){
    fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_kcp3o_simsize_",sim.size,sep="")
    fileName1<-paste(dirr,fileName,".Rda",sep="")
    load(fileName1)
    
    
    vals<-getDist()
    
    output<-cbind(output,vals)
    
    cols=c(cols,params[3])
  }
}

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
makeBP(output)

save_dirr=""
Cairo::CairoPDF(file=paste(save_dirr,"ecp_kcp3o_kpb.pdf",sep=""),height=6,width=10)
makeBP(output)
dev.off()









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


output<-NULL
cols=NULL
dimss=NULL
for(i in 1:36){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  dName=names(distributions1)[params[3]]
  
  # if(d==10&&distr==2){
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_kcp3o_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  load(fileName1)
  
  
  vals<-getDist()
  
  output<-cbind(output,vals)
  
  cols=c(cols,rep( params[3],dim(vals)[2]))
  #reorder by dim for colors
  dimss=c(dimss,rep( params[4],dim(vals)[2]))
  # }
}



bor=paramterIndices[order(paramterIndices[1:36,4]),][,3]
dimss

a<-c(colors()[617],"black","violetred")
bor=a[bor]


Cairo::CairoPDF(file=paste(save_dirr,"kbp_ecp_kcp3o_1000.pdf",sep=""),height=7.5,width=14)
makeBigPlot(output[,order(dimss)])
title("N=1000, Scenario 1, KCP3O")
dev.off()










#########N=2500

output<-NULL
cols=NULL
dimss=NULL
for(i in 37:72){
  
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  dName=names(distributions1)[params[3]]
  
  # if(d==10&&distr==2){
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_kcp3o_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  load(fileName1)
  
  
  vals<-getDist()
  
  output<-cbind(output,vals)
  
  cols=c(cols,rep( params[3],dim(vals)[2]))
  #reorder by dim for colors
  dimss=c(dimss,rep( params[4],dim(vals)[2]))
  # }
}



bor=paramterIndices[order(paramterIndices[1:36,4]),][,3]
dimss

a<-c(colors()[617],"black","violetred")
bor=a[bor]


Cairo::CairoPDF(file=paste(save_dirr,"kbp_ecp_kcp3o_2500.pdf",sep=""),height=7.5,width=14)
makeBigPlot(output[,order(dimss)])
title("N=2500, Scenario 1, KCP3O")
dev.off()










