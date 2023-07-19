

#packages
library(mvtnorm)
library(fMultivar)
library(MASS)
library(sde)
library(rrcov)
library(xtable)
library(sn)
library(doParallel)
library(doRNG)
library(ecp)

#results dirrectory
dirr<- "ECP_DIV"
setwd(dirr)


#simulation parameters
#N
Ns<-c(100,200)
sim.size=100
thetas<-list(c(.333,.666), 
             c(.25,.5,.75),
             c(1/6,2/6,3/6,4/6,5/6))
ds=c(2,3,5,10)
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



#parallel running
#necessary packages
pkgs<-c("mvtnorm" ,  "fMultivar" ,"MASS"     , "sde"      , "rrcov"  ,   "sn"  ,'ecp' )

no_cores<-Sys.getenv("SLURM_CPUS_PER_TASK") 

cl <-  makeCluster(no_cores,type="FORK")

registerDoParallel(cl) 
registerDoRNG(seed = 440)



#run simulation
for(i in 1:nrow(paramterIndices)){
  
  #set paramters
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions[[params[3]]]
  d=ds[params[4]]
  dName=distr
  
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
  runSim<-function(N,theta,dist,d){
    

    #sim data
    testData<-simData(N,theta,dist,d)
    
    ecp_div=ecp::e.divisive(testData)$estimates# long 16.4 seconds
    return(ecp_div)
  }
  
  

  
  results=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d)}})
  
  error=inherits(results, "try-error")
  if(error){
    print("there was an error!")
    
  }
  

  
  
  
  
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_ecp_div_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  
  save(results,file=fileName1)
  
 
  print(i/numUniqueRuns)
}

registerDoSEQ()
stopCluster(cl)
closeAllConnections()





