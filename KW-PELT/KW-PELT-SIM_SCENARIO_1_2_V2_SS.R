# 
# args = commandArgs(trailingOnly=TRUE)
# 
# # test if there is at least two arguments: if not, return an error
# if (length(args)<1) {
#   stop("At least one arguments must be supplied ('name' (text) and 'numer' (integer) )", 
#        call.=FALSE)
# }
# 
# Ns      <- as.numeric(args[1])             r
# 
# print(paste("Processing with N: ", Ns, sep = ''))

# 
# install.packages(c('mvtnorm',
#                    'fMultivar',
#                    'MASS',
#                    'sde',
#                    'rrcov',
#                    'xtable',
#                    'sn',
#                    'ddalpha',
#                    'doParallel',
#                    'doRNG',
#                    'Rcpp',
#                    'RcppArmadillo'))

# /tmp/RtmpuL1uEw/downloaded_packages


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

# dir.create('PELT')
#results directory
dirr<-"PELT"
dirr2="PELT"
#PELT results , after running pelt on depths
setwd(dirr)


#simulation parameters

#N
Ns<-c(100,200,1000,2500,5000)
Ns=c(100,200)
#simulation size, num repetitions
sim_size=100


#2,3,5 change-points
thetas<-list(c(.333,.666), 
             c(.25,.5,.75),
             c(1/6,2/6,3/6,4/6,5/6))


#dimension
ds=c(2,3,5,10)

#threshold parameter
thresh<-1

#distributions simulated from
distributions<-c("Normal", "Cauchy")



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






#run simulation
no_cores<-Sys.getenv("SLURM_CPUS_PER_TASK") 

cl <-  makeCluster(no_cores,type="FORK")

registerDoParallel(cl) 
registerDoRNG(seed = 440)
pkgs<-c("mvtnorm" ,  "fMultivar" ,"MASS"     , "sde"      , "rrcov"  ,   "sn" ,       "ddalpha"  )

for(j in 1:numUniqueRuns){
  
  #set paramters
  params<-paramterIndices[j,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions[[params[3]]]
  d=ds[params[4]]

  
  
  
  
  #simulate data set for one repitition
  simData<-function(N,theta,dist,d){
    
    #simulation functions, used to generate data, simulate the change size
    change_sizes=sample(1:10,5,FALSE)
    
    #generates iid d-normal, sigmaxI Rv
    normalMaster<-function(n,d,sigmaSq){mvtnorm::rmvnorm(n,sigma=diag(sigmaSq*rep(1,d)))}
    
    norm1<-function(n,d){normalMaster(n,d,change_sizes[1])}
    norm2<-function(n,d){normalMaster(n,d,change_sizes[2])}
    norm3<-function(n,d){normalMaster(n,d,change_sizes[3])}
    norm4<-function(n,d){normalMaster(n,d,change_sizes[4])}
    norm5<-function(n,d){normalMaster(n,d,change_sizes[5])}
    
    normals<-list(norm1,norm2,norm3,norm4,norm5,norm1)
    
    
    #generates iid d-cauchy, sigmaxI scale
    cauchyMaster<-function(n,d,sigmaSq){replicate(d,rcauchy(n,scale=sigmaSq))}
    
    
    cauchy1<-function(n,d){cauchyMaster(n,d,change_sizes[1])}
    cauchy2<-function(n,d){cauchyMaster(n,d,change_sizes[2])}
    cauchy3<-function(n,d){cauchyMaster(n,d,change_sizes[3])}
    cauchy4<-function(n,d){cauchyMaster(n,d,change_sizes[4])}
    cauchy5<-function(n,d){cauchyMaster(n,d,change_sizes[5])}
    
    cauchys<-list(cauchy1,cauchy2,cauchy3,cauchy4,cauchy5,cauchy1)
    
    distributions<-list(normals,cauchys)
    
    #cp locations
    if(!is.null(theta))
      locations<-c(1,floor(theta*N),N)
    else
      locations<-N
    
    #data
    dat<-matrix(0,nrow=N,ncol=d)
    if(dist=="Normal"){
      rdata=normals
    }
    else{
      rdata=cauchys
    }
    for(i in 2:(length(theta)+2))
      dat[locations[i-1]:locations[i],]<-rdata[[i-1]](length(locations[i-1]:locations[i]),d)
    
    return(dat)
  }
  
  
  
  
  #run one repetition
  runSim<-function(N,theta,dist,d,depth="hs"){
    
   
    
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
    
    testData<-simData(N,theta,dist,d)
    
    return(getDepths(testData,depth))
  }
  

  
  #  set.seed(746383)
  depthsHS=try({foreach(i=1:sim_size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="hs")}})
  
  errorhs=inherits(depthsHS, "try-error")
  
  if(errorhs){
    print("there was an error!")
    
  }
  # stopCluster(cl)
  # 
  # 
  # cl <-  makeCluster(no_cores,type="FORK")  
  # 
  # registerDoParallel(cl) 
  # registerDoRNG(seed = 440)
  
  #  set.seed(746383)
  depthsSpat=try({foreach(i=1:sim_size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="spat")}})
  
  errorsp=inherits(depthsSpat, "try-error")
  if(errorsp){
    print("there was an error!")
    
  }
  
  # stopCluster(cl)
  
  
  # 
  # cl <-  makeCluster(no_cores,type="FORK")  
  # 
  # 
  # registerDoParallel(cl)
  # registerDoRNG(seed = 440)
  #set.seed(746383)
  depthsMahal=try({foreach(i=1:sim_size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="mahal")}})
  
  errormh=inherits(depthsMahal, "try-error")
  if(errormh){
    print("there was an error!")
    
  }
  
  
  
  
  # 
  # stopCluster(cl)
  # 
  # 
  # cl <-  makeCluster(no_cores,type="FORK")  
  # 
  # 
  # registerDoParallel(cl)
  # registerDoRNG(seed = 440)
  # set.seed(746383)
  depthsMahal75=try({foreach(i=1:sim_size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,depth="mahal75")}})
  
  errormh75=inherits(depthsMahal75, "try-error")
  if(errormh75){
    print("there was an error!")
    
  }
  
  
  # stopCluster(cl)
  #OC stands for scenario 1
  fileName<-paste0(N,"_",length(theta),"_",distr,"_",d,"_","OC_SS_Depths_simsize_",sim_size,sep="")
  
  
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  
  save(depthsHS,depthsSpat,depthsMahal,depthsMahal75,file=fileName1)
  
  # closeAllConnections()
  print(j)
}

stopCluster(cl)
closeAllConnections()

# 


