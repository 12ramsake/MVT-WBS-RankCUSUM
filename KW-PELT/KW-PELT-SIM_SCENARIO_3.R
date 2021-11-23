

#packages
library(mvtnorm)
library(fMultivar)
library(MASS)
library(sde)
library(rrcov)
library(xtable)
library(sn)
library(ddalpha)
library(doParallel)
library(doRNG)
library(Rcpp)
library(RcppArmadillo)
#results dirrectory
dirr<-"/u/k3ramsay/ResearchDocuments/output/KW_PELT_SIMULATION/Scen_3/"
dirr2<-"/u/k3ramsay/ResearchDocuments/output/KW_PELT_SIMULATION/Scen_3/"
setwd(dirr)


#simulation parameters
Ns <- c(1000)

#simulation size, num repetitions
sim.size = 100




#xhange below
#0,2,3.5 cps
thetas <- list(c(.333, .666),
               c(.25, .5, .75),
               c(1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6))


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








#simulation functions

#generates iid d-normal, sigmaxI Rv
normalMaster<-function(n,d,sigmaSq,d2){mvtnorm::rmvnorm(n,sigma=diag(c(sigmaSq*rep(1,d2),rep(1,d-d2))))}

norm1<-function(n,d,d2){normalMaster(n,d,1,d2)}
norm2<-function(n,d,d2){normalMaster(n,d,2.5,d2)}
norm3<-function(n,d,d2){normalMaster(n,d,4,d2)}
norm4<-function(n,d,d2){normalMaster(n,d,2.25,d2)}
norm5<-function(n,d,d2){normalMaster(n,d,5,d2)}

normals<-list(norm1,norm2,norm3,norm4,norm5,norm1)




for(i in 1:numUniqueRuns){
  params<-paramterIndices[i,]
  N=1000
  theta=thetas[[params[2]]]
  d2=d2s[params[3]]
  d=5
  dName="Normal"
  distr=normals
  
  
  #run one repetition
  runSim<-function(N,theta,rdata,d,depth="hs",d2){
    
    #simulate data set for one repitition
    simData<-function(N,theta,rdata,d,d2){
      
      #cp locations
      if(!is.null(theta))
        locations<-c(1,floor(theta*N),N)
      else
        locations<-N
      
      #data
      dat<-matrix(0,nrow=N,ncol=d)
      
      for(i in 2:(length(theta)+2))
        dat[locations[i-1]:locations[i],]<-rdata[[i-1]](length(locations[i-1]:locations[i]),d,d2)
      
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
    
    testData<-simData(N,theta,rdata,d,d2)
    
    return(getDepths(testData,depth))
  }
  
  
  pkgs<-c("mvtnorm" ,  "fMultivar" ,"MASS"     , "sde"      , "rrcov"  ,   "sn" ,       "ddalpha"  )
  
  no_cores<-detectCores()-1
  no_cores<-100
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  depthsHS=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="hs",d2)}})
  
  errorhs=inherits(depthsHS, "try-error")
  
  if(errorhs){
    print("there was an error!")
    
  }
  stopCluster(cl)
  
  closeAllConnections()
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  
  
  depthsSpat=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="spat",d2)}})
  
  errorsp=inherits(depthsSpat, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  stopCluster(cl)
  closeAllConnections()
  
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  
  registerDoParallel(cl)
  registerDoRNG(seed = 440)
  
  depthsMahal=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="mahal",d2)}})
  
  errormh=inherits(depthsMahal, "try-error")
  if(errormh){
    print("there was an error!")
    
  }
  
  
  
  
  
  stopCluster(cl)
  closeAllConnections()
  
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  
  registerDoParallel(cl)
  registerDoRNG(seed = 440)
  
  depthsMahal75=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="mahal75",d2)}})
  
  errormh75=inherits(depthsMahal75, "try-error")
  if(errormh75){
    print("there was an error!")
    
  }
  
  
  stopCluster(cl)
  
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_d2_",d2,"_Scen_3_Depths_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  
  save(depthsHS,depthsSpat,depthsMahal,depthsMahal75,file=fileName1)
  
  closeAllConnections()
  print(i/numUniqueRuns)
}



##################### Now run PELT

Rcpp::sourceCpp('PELT_CPP.cpp')













constant=0.18
for(i in 1:numUniqueRuns){
  
  params<-paramterIndices[i,]
  N=1000
  theta=thetas[[params[2]]]
  d2=d2s[params[3]]
  d=5
  dName="Normal"
  beta=constant*sqrt(N)+3.74
  
  
  
  # distr=distributions[[1]]
  
  
  
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_d2_",d2,"_Scen_3_Depths_simsize_",sim.size,sep="")
  
  
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  load(fileName1)
  
  
  
  ranksHS<-lapply(depthsHS,rank)
  ranksSpat<-lapply(depthsSpat,rank)
  ranksMahal<-lapply(depthsMahal,rank)
  ranksMahal75<-lapply(depthsMahal75,rank)
  
  
  # Rank_PELT(ranksHS[[1]],beta=10)
  # beta=15
  resultsHS<-PELT_T(ranksHS[[1]],beta=beta)
  no_cores<-100
  
  cl <- makeCluster(no_cores,type="FORK")   
  registerDoParallel(cl) 
  
  resultsHS=try({foreach(i=1:length(ranksSpat),.packages= pkgs,.noexport = c("PELT_T")) %dopar% {PELT_T(ranksHS[[i]],beta=beta)}})
  
  errorsp=inherits(resultsSpat, "try-error")
  if(errorsp){
    print("there was an error! HS")
    
  }
  
  registerDoSEQ()
  stopCluster(cl)
  closeAllConnections()
  
  
  cl <- makeCluster(no_cores,type="FORK")   
  registerDoParallel(cl) 
  
  resultsSpat=try({foreach(i=1:length(ranksSpat),.packages= pkgs,.noexport = c("PELT_T")) %dopar% {PELT_T(ranksSpat[[i]],beta=beta)}})
  
  errorsp=inherits(resultsSpat, "try-error")
  if(errorsp){
    print("there was an error! Spat")
    
  }
  
  registerDoSEQ()
  stopCluster(cl)
  closeAllConnections()
  
  
  cl <- makeCluster(no_cores,type="FORK")   
  registerDoParallel(cl) 
  
  resultsMH=try({foreach(i=1:length(ranksSpat),.packages= pkgs,.noexport = c("PELT_T")) %dopar% {PELT_T(ranksMahal[[i]],beta=beta)}})
  
  errorsp=inherits(resultsSpat, "try-error")
  if(errorsp){
    print("there was an error! MH")
    
  }
  
  registerDoSEQ()
  stopCluster(cl)
  closeAllConnections()
  
  cl <- makeCluster(no_cores,type="FORK")   
  registerDoParallel(cl) 
  
  resultsMH75=try({foreach(i=1:length(ranksSpat),.packages= pkgs,.noexport = c("PELT_T")) %dopar% {PELT_T(ranksMahal75[[i]],beta=beta)}})
  
  errorsp=inherits(resultsSpat, "try-error")
  if(errorsp){
    print("there was an error! MH75")
    
  }
  
  registerDoSEQ()
  stopCluster(cl)
  closeAllConnections()
  
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_constant_C1_",constant,"_d2_",d2,"_Scen_3_PELT_ranks_simsize_",sim.size,sep="")
  fileName1<-paste(dirr2,fileName,".Rda",sep="")
  save(resultsHS,resultsSpat,resultsMH,resultsMH75,file=fileName1)
  
  closeAllConnections()
  
}














