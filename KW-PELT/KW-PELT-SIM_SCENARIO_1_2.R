

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


#results directory
dirr<-""
dirr2=""
#PELT results , after running pelt on depths
setwd(dirr)


#simulation parameters

#N
Ns<-c(100,200,1000,2500,5000)
#simulation size, num repetitions
sim.size=100
#2,3.5 cps
thetas<-list(c(.333,.666), 
             c(.25,.5,.75),
             c(1/6,2/6,3/6,4/6,5/6))
##d= 2,3,5,10
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

skewNormalMaster<-function(n,d,sigmaSq,skewParam){
  rmsn(n, dp=cp2dp(list(mean=rep(0,d), 
                        var.cov=diag(sigmaSq*rep(1,d)), 
                        gamma1=skewParam*rep(1,d)/d), "SN"))
}

skewNorm1<-function(n,d){skewNormalMaster(n,d,1,.1)}
skewNorm2<-function(n,d){skewNormalMaster(n,d,2.5,.1)}
skewNorm3<-function(n,d){skewNormalMaster(n,d,4,.1)}
skewNorm4<-function(n,d){skewNormalMaster(n,d,2.25,.1)}
skewNorm5<-function(n,d){skewNormalMaster(n,d,5,.1)}

skewNormals<-list(skewNorm1,skewNorm2,skewNorm3,skewNorm4,skewNorm5,skewNorm1)

distributions<-list(normals,cauchys,skewNormals)
names(distributions)=c("Normal","Cauchy","Skew Normal")




#######Simulates the data and computes the depth values, then saves them to a file. 
for(i in 1:108){
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
  
  no_cores<-100
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  
  #  set.seed(746383)
  depthsHS=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="hs")}})
  
  errorhs=inherits(depthsHS, "try-error")
  
  if(errorhs){
    print("there was an error!")
    
  }
  stopCluster(cl)
  
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  
  #  set.seed(746383)
  depthsSpat=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="spat")}})
  
  errorsp=inherits(depthsSpat, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  stopCluster(cl)
  
  
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  
  registerDoParallel(cl)
  registerDoRNG(seed = 440)
  #set.seed(746383)
  depthsMahal=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="mahal")}})
  
  errormh=inherits(depthsMahal, "try-error")
  if(errormh){
    print("there was an error!")
    
  }
  
  
  
  
  
  stopCluster(cl)
  
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  
  registerDoParallel(cl)
  registerDoRNG(seed = 440)
  # set.seed(746383)
  depthsMahal75=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="mahal75")}})
  
  errormh75=inherits(depthsMahal75, "try-error")
  if(errormh75){
    print("there was an error!")
    
  }
  
  
  stopCluster(cl)
  #OC stands for scenario 1
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_","OC_Depths_simsize_",sim.size,sep="")
  
  
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  
  save(depthsHS,depthsSpat,depthsMahal,depthsMahal75,file=fileName1)
  
  closeAllConnections()
  print(i)
}






#######################NOW RUN PELT 




Rcpp::sourceCpp('PELT_CPP.cpp')



constants<-c(0.1,0.14,0.16,0.18,0.2,0.24,0.3)

n_cores=100
n_cores=detectCores()-4
n_cores



  # print(beta)
  
for(i in 1:numUniqueRuns){
  for(constant in constants){
    params<-paramterIndices[i,]
    N=Ns[params[1]]
    beta=constant*sqrt(N)+3.74
    theta=thetas[[params[2]]]
    distr=distributions1[[params[3]]]
    d=ds[params[4]]
    
    dName=names(distributions1)[params[3]]
    distr=distributions[[distr]]
    
    
    
    fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_","OC_Depths_simsize_",sim.size,sep="")
    
    
    fileName1<-paste(dirr,fileName,".Rda",sep="")
    load(fileName1)
    
    
    
    ranksHS<-lapply(depthsHS,rank)
    ranksSpat<-lapply(depthsSpat,rank)
    ranksMahal<-lapply(depthsMahal,rank)
    ranksMahal75<-lapply(depthsMahal75,rank)
    
    
    # Rank_PELT(ranksHS[[1]],beta=10)
    # beta=15
    resultsHS<-mclapply(ranksHS, PELT_T,beta=beta,mc.cores=n_cores)
    resultsSpat<-mclapply(ranksSpat, PELT_T,beta=beta,mc.cores=n_cores)
    resultsMahal<-mclapply(ranksMahal, PELT_T,beta=beta,mc.cores=n_cores)
    resultsMahal75<-mclapply(ranksMahal75, PELT_T,beta=beta,mc.cores=n_cores)
    
    fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_constant_C1_",constant,"OC_PELT_ranks_simsize_",sim.size,sep="")
    fileName1<-paste(dirr2,fileName,".Rda",sep="")
    save(resultsHS,resultsSpat,resultsMahal,resultsMahal75,file=fileName1)
    
    closeAllConnections()
    print(i/numUniqueRuns+(which(constant==constants)-1)/length(constants)/numUniqueRuns)

  }
}




###################################################### Scenario 2  ###################################################### 

#simulation functions

#generates iid d-normal, sigmaxI Rv
normalMaster<-function(n,d,sigmaSq){mvtnorm::rmvnorm(n,sigma=diag(sigmaSq*rep(1,d)))}

norm1<-function(n,d){normalMaster(n,d,1)}
norm2<-function(n,d){normalMaster(n,d,3)}
norm3<-function(n,d){normalMaster(n,d,5)}
norm4<-function(n,d){normalMaster(n,d,3)}
norm5<-function(n,d){normalMaster(n,d,5)}

normals<-list(norm1,norm2,norm3,norm4,norm5,norm1)


#generates iid d-cauchy, sigmaxI scale
cauchyMaster<-function(n,d,sigmaSq){replicate(d,rcauchy(n,scale=sigmaSq))}


cauchy1<-function(n,d){cauchyMaster(n,d,1)}
cauchy2<-function(n,d){cauchyMaster(n,d,3)}
cauchy3<-function(n,d){cauchyMaster(n,d,5)}
cauchy4<-function(n,d){cauchyMaster(n,d,3)}
cauchy5<-function(n,d){cauchyMaster(n,d,5)}

cauchys<-list(cauchy1,cauchy2,cauchy3,cauchy4,cauchy5,cauchy1)

skewNormalMaster<-function(n,d,sigmaSq,skewParam){
  rmsn(n, dp=cp2dp(list(mean=rep(0,d), 
                        var.cov=diag(sigmaSq*rep(1,d)), 
                        gamma1=skewParam*rep(1,d)/d), "SN"))
}

skewNorm1<-function(n,d){skewNormalMaster(n,d,1,.1)}
skewNorm2<-function(n,d){skewNormalMaster(n,d,3,.1)}
skewNorm3<-function(n,d){skewNormalMaster(n,d,5,.1)}
skewNorm4<-function(n,d){skewNormalMaster(n,d,3,.1)}
skewNorm5<-function(n,d){skewNormalMaster(n,d,5,.1)}

skewNormals<-list(skewNorm1,skewNorm2,skewNorm3,skewNorm4,skewNorm5,skewNorm1)

distributions<-list(normals,cauchys,skewNormals)
names(distributions)=c("Normal","Cauchy","Skew Normal")


for(i in 1:108){
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
  no_cores<-detectCores()-5
  
  # no_cores<-100
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  
  #  set.seed(746383)
  depthsHS=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="hs")}})
  
  errorhs=inherits(depthsHS, "try-error")
  
  if(errorhs){
    print("there was an error!")
    
  }
  stopCluster(cl)
  
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  registerDoParallel(cl) 
  registerDoRNG(seed = 440)
  
  #  set.seed(746383)
  depthsSpat=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="spat")}})
  
  errorsp=inherits(depthsSpat, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  stopCluster(cl)
  
  
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  
  registerDoParallel(cl)
  registerDoRNG(seed = 440)
  #set.seed(746383)
  depthsMahal=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="mahal")}})
  
  errormh=inherits(depthsMahal, "try-error")
  if(errormh){
    print("there was an error!")
    
  }
  
  
  
  
  
  stopCluster(cl)
  
  
  cl <-  makeCluster(no_cores,type="FORK")  
  
  
  registerDoParallel(cl)
  registerDoRNG(seed = 440)
  # set.seed(746383)
  depthsMahal75=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="mahal75")}})
  
  errormh75=inherits(depthsMahal75, "try-error")
  if(errormh75){
    print("there was an error!")
    
  }
  
  
  stopCluster(cl)
  
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_","NC_Depths_simsize_",sim.size,sep="")
  
  
  fileName1<-paste(dirr,fileName,".RData",sep="")
  
  save(depthsHS,depthsSpat,depthsMahal,depthsMahal75,file=fileName1)
  
  closeAllConnections()
  print(i)
}


#######################NOW RUN PELT 




Rcpp::sourceCpp('PELT_CPP.cpp')



constants<-c(0.1,0.14,0.16,0.18,0.2,0.24,0.3)

n_cores=100
n_cores=detectCores()-4
n_cores



#Scenario 2, for NEW changes

for(i in 1:numUniqueRuns){
  for(constant in constants){
    
    params<-paramterIndices[i,]
    N=Ns[params[1]]
    beta=constant*sqrt(N)+3.74
    theta=thetas[[params[2]]]
    distr=distributions1[[params[3]]]
    d=ds[params[4]]
    
    dName=names(distributions1)[params[3]]
    distr=distributions[[distr]]
    
    
    
    
    fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_","NC_Depths_simsize_",sim.size,sep="")
    
    
    fileName1<-paste(dirr,fileName,".RData",sep="")
    load(fileName1)
    
    
    
    ranksHS<-lapply(depthsHS,rank)
    ranksSpat<-lapply(depthsSpat,rank)
    ranksMahal<-lapply(depthsMahal,rank)
    ranksMahal75<-lapply(depthsMahal75,rank)
    
    
    resultsHS<-mclapply(ranksHS, PELT_T,beta=beta,mc.cores=n_cores)
    resultsSpat<-mclapply(ranksSpat, PELT_T,beta=beta,mc.cores=n_cores)
    resultsMahal<-mclapply(ranksMahal, PELT_T,beta=beta,mc.cores=n_cores)
    resultsMahal75<-mclapply(ranksMahal75, PELT_T,beta=beta,mc.cores=n_cores)
    
    #Scenario 2, for NEW changes
    fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_constant_C1_",constant,"NC_PELT_ranks_simsize_",sim.size,sep="")
    fileName1<-paste(dirr2,fileName,".Rda",sep="")
    save(resultsHS,resultsSpat,resultsMahal,resultsMahal75,file=fileName1)
    
    closeAllConnections()
    print(i/numUniqueRuns/length(constants)+(which(constant==constants)-1)/length(constants))
  }
}









