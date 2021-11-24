library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
library(fMultivar)
library(MASS)
library(sde)
library(rrcov)
library(xtable)
library(sn)
library(ddalpha)
library(doParallel)

dirr<-"/u/k3ramsay/ResearchDocuments/output/KW_PELT_SIMULATION/Scen_4/"
dirr2<-"/u/k3ramsay/ResearchDocuments/output/KW_PELT_SIMULATION/Scen_4/"
setwd(dirr)

Rcpp::sourceCpp('PELT_CPP.cpp')
Ns<-c(1000)
sim.size=100

################################################# ################################################# ################################################# 
################################################# ################################################# ################################################# 
################################################# ################################################# ################################################# 
################################################# Increasing Dimension ################################################# 
################################################# ################################################# ################################################# 
################################################# ################################################# ################################################# 


thetas<-list(.5,c(.333,.666))


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

#generates iid d-normal, sigmaxI Rv
normalMaster<-function(n,d,sigmaSq){mvtnorm::rmvnorm(n,sigma=diag(sigmaSq*rep(1,d)))}

norm1<-function(n,d){normalMaster(n,d,1)}
norm2<-function(n,d){normalMaster(n,d,2.5)}
norm3<-function(n,d){normalMaster(n,d,4)}
norm4<-function(n,d){normalMaster(n,d,2.25)}
norm5<-function(n,d){normalMaster(n,d,5)}

normals<-list(norm1,norm2,norm3,norm4,norm5,norm1)


distributions<-list(normals)
names(distributions)=c("Normal")



for(i in 1:numUniqueRuns){
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  dName=names(distributions1)[params[3]]
  
  
  
  
  distr=distributions[[distr]]
  
  
  #run one repetition
  runSim<-function(N,theta,rdata,d,depth="hs"){
    
    #simulate data set for one repitition
    simData<-function(N,theta,rdata,d){
      
      #cp locations
      if(!is.null(theta))
        locations<-c(1,floor(theta*N),N)
      else
        locations<-N
      
      #data
      dat<-matrix(0,nrow=N,ncol=d)
      
      for(i in 2:(length(theta)+2))
        dat[locations[i-1]:locations[i],]<-rdata[[i-1]](length(locations[i-1]:locations[i]),d)
      
      return(dat)
    }
    
    #get depth
    #spat, hs, mahal,mahal75
    getDepths<-function(data,depth){
      
      if(depth=="spat"){
        ts=depth.spatial(data,data)
      }
      else if(depth=="hs"){
        ts=depth.halfspace(data,data)
      }
      else if(depth=="mahal"){
        ts=depth.Mahalanobis(data,data)
      }
      else if(depth=="mahal75"){
        ts=depth.Mahalanobis(data,data, "MCD")
      }
      else{
        ts=NULL
        print("bad depth specification")
      }
      
      return(ts)
    }
    
    testData<-simData(N,theta,rdata,d)
    
    return(getDepths(testData,depth))
  }
  
  
  pkgs<-c("mvtnorm" ,  "fMultivar" ,"MASS"     , "sde"      , "rrcov"  ,   "sn" ,       "ddalpha"  )
  # no_cores<-detectCores()-1
  
  no_cores<-detectCores()-10
  no_cores<-100
  
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  depthsSpat=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="spat")}})
  
  errorsp=inherits(depthsSpat, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  stopCluster(cl)
  
  
  
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_","HD_Depths_simsize_",sim.size,sep="")
  
  
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  
  save(depthsSpat,file=fileName1)
  
  closeAllConnections()
  print(i)
}



constant=0.18


for(i in 1:numUniqueRuns){
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  d=ds[params[4]]
  distr=distributions1[[params[3]]]
  dName=names(distributions1)[params[3]]
  distr=distributions[[distr]]
  
  beta=constant*sqrt(N)+3.74
  
  
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_","HD_Depths_simsize_",sim.size,sep="")
  
  
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  load(fileName1)
  
  
  
  ranksSpat<-lapply(depthsSpat,rank)
  
  
  resultsSpat<-lapply(ranksSpat, PELT_T,beta=beta)
  
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_constant_",constant,"HD_PELT_ranks_simsize_",sim.size,sep="")
  fileName1<-paste(dirr2,fileName,".Rda",sep="")
  save(resultsSpat,file=fileName1)
  
  closeAllConnections()
}










################################################# ################################################# ################################################# 
################################################# ################################################# ################################################# 
################################################# ################################################# ################################################# 
################################################# Increasing Sparsity ################################################# 
################################################# ################################################# ################################################# 
################################################# ################################################# ################################################# 



sim.size=100

ds=c(5,10,50,100,250,500)


theta<-0.5
dim_of_change=5
N=1000
dName="Normal"
distr=1


for(d in ds){
  
  
  
  norm1<-function(n,d){mvtnorm::rmvnorm(n,sigma=diag(rep(1,d)))}
  norm2<-function(n,d){mvtnorm::rmvnorm(n,sigma=diag(c(2.5*rep(1,dim_of_change),rep(1,d-dim_of_change))))}
  distr=list(norm1,norm2)
  
  
  #run one repetition
  runSim<-function(N,theta,rdata,d,depth="hs"){
    
    #simulate data set for one repitition
    simData<-function(N,theta,rdata,d){
      
      #cp locations
      if(!is.null(theta))
        locations<-c(1,floor(theta*N),N)
      else
        locations<-N
      
      #data
      dat<-matrix(0,nrow=N,ncol=d)
      
      for(i in 2:(length(theta)+2))
        dat[locations[i-1]:locations[i],]<-rdata[[i-1]](length(locations[i-1]:locations[i]),d)
      
      return(dat)
    }
    
    #get depth
    #spat, hs, mahal,mahal75
    getDepths<-function(data,depth){
      
      if(depth=="spat"){
        ts=depth.spatial(data,data)
      }
      else if(depth=="hs"){
        ts=depth.halfspace(data,data)
      }
      else if(depth=="mahal"){
        ts=depth.Mahalanobis(data,data)
      }
      else if(depth=="mahal75"){
        ts=depth.Mahalanobis(data,data, "MCD")
      }
      else{
        ts=NULL
        print("bad depth specification")
      }
      
      return(ts)
    }
    
    testData<-simData(N,theta,rdata,d)
    
    return(getDepths(testData,depth))
  }
  
  
  pkgs<-c("mvtnorm" ,  "fMultivar" ,"MASS"     , "sde"      , "rrcov"  ,   "sn" ,       "ddalpha"  )
  # no_cores<-detectCores()-1
  
  no_cores<-detectCores()-10
  no_cores<-100
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  
  
  depthsSpat=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="spat")}})
  
  errorsp=inherits(depthsSpat, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  stopCluster(cl)
  
  
  
  fileName<-paste0(N,"_",length(theta),"_Normal_d_",d,"_","HD3_dimchange_",dim_of_change,"_Depths_simsize_",sim.size,sep="")
  
  
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  
  save(depthsSpat,file=fileName1)
  
  closeAllConnections()
  print(dim_of_change)
  print(d)
  
}



################PELT
constant=0.18
sim.size=100

ds=c(5,10,50,100,250,500)


theta<-0.5
dim_of_change=5
N=1000
dName="Normal"
distr=1

beta=constant*sqrt(N)+3.74


for(d in ds){
  
  
  
  
  fileName<-paste0(N,"_",length(theta),"_Normal_d_",d,"_","HD3_dimchange_",dim_of_change,"_Depths_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  load(fileName1)
  
  
  
  ranksSpat<-lapply(depthsSpat,rank)
  
  
  resultsSpat<-mclapply(ranksSpat, PELT_T,beta=beta, mc.cores =  no_cores)
  
  fileName<-paste0(N,"_",length(theta),"_Normal_d_",d,"_HD3_dimchange_",dim_of_change,"_constant_",constant,"_PELT_ranks_simsize_",sim.size,sep="")
  fileName1<-paste(dirr2,fileName,".Rda",sep="")
  save(resultsSpat,file=fileName1)
  
  closeAllConnections()
}


