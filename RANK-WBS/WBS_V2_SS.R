

#Seed
# set.seed(746383)
#packages
library(mvtnorm)
library(fMultivar)
library(MASS)
library(sde)
# library(depth)
library(rrcov)
library(xtable)
library(sn)
library(mrfDepth)
library(ddalpha)
library(doParallel)
library(doRNG)

#results dirrectory
dirr<- "WBS"
setwd(dirr)


#simulation parameters

#N- samples sizes
Ns<-c(100,200)

#Number of intervals, mod*logN intervals sampled
mod=100
numInts<-floor(log(Ns))*mod


#simulation size, num repetitions of the experiment
sim.size=100


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

for(i in 1:numUniqueRuns){
  
  #set paramters
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions[[params[3]]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*mod
  
  
  
  
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
  runSim<-function(N,theta,dist,d,numInt=10,thresh=1.3584,depth="hs"){
  
    
    
    ##wild binary segmentation
    #calculate the test statistic
    #spat, hs, mahal,mahal75 are the depth parameters
    testStat<-function(range,data,depth){
      
      if(depth=="spat"){
        ts=testStatSpat(range,data)
      }
      else if(depth=="hs"){
        ts=testStatHs(range,data)
      }
      else if(depth=="mahal"){
        ts=testStatMahal(range,data)
      }
      else if(depth=="mahal75"){
        ts=testStatMahal75(range,data)
      }
      else{
        ts=NULL
        print("bad depth specification")
      }
      
      return(ts)
    }
    
    #test cusum from depth values
    getStatFromDepths<-function(depths,N){
      ranks<-rank(depths,ties.method = "random")
      expected.val<-(N-1)/2
      std.dev<-sqrt((N^2-1)/12)  
      cusum<-cumsum(N^(-0.5)*(ranks-expected.val)/std.dev)
      return(abs(cusum)[1:(length(cusum)-1)])
    }
    
    testStatHs<-function(range,data){
      if((range[2]-range[1])>(ncol(data)+1)){
        range<-range[1]:range[2]
        N<-nrow(data[range,])
        depths<-depth.halfspace(data[range,],data[range,])
        return(getStatFromDepths(depths,N))
      }
      else
        return(-10)
    }
    
    testStatSpat<-function(range,data){
      if((range[2]-range[1])>(ncol(data)+1)){
        range<-range[1]:range[2]
        N<-nrow(data[range,])
        depths<-depth.spatial(data[range,],data[range,])
        
        return(getStatFromDepths(depths,N))
      }
      else
        return(-10)
    }
    
    testStatMahal75<-function(range,data){
      
      if((range[2]-range[1])>(ncol(data)*2)){
        range<-range[1]:range[2]
        
        N<-nrow(data[range,])
        depths<-depth.Mahalanobis(data[range,],data[range,], "MCD")
        
        return(getStatFromDepths(depths,N))
      }
      else
        return(-10)
    }
    
    testStatMahal<-function(range,data){
      if((range[2]-range[1])>(ncol(data)*2)){
        range<-range[1]:range[2]
        
        N<-nrow(data[range,])
        depths<-depth.Mahalanobis(data[range,],data[range,])
        return(getStatFromDepths(depths,N))
      }
      else
        return(-10)
    }
    
    
    #returns indices of the intervals selected, M is the number of intervals
    getIntervals<-function(indices,M){
      
      ints<-t(replicate(M,sort(sample(indices,2))))
      diffs<-(ints[,2]-ints[,1])==1
      if(any(diffs)){
        ints[diffs,]=getIntervals(indices,sum(diffs))
        return(ints)
      }
      else{
        return(ints)
      }
    }
    
    #checks if an interval is a sub of another
    checkIfSubInterval<-function(sub,super){
      return(sub[1]>=super[1]&&sub[2]<=super[2])
    }
    
    #let the set of intervals be a matrix with 2 columns
    
    WBS<-function(intervals,s,e,threshold,data){
      
      
      if((e-s)<1)
        return(NULL)
      
      else{
        #intervals contained in s,e
        Mes<-which(apply(intervals,1,checkIfSubInterval,super=c(s,e)))
        
        
        if(length(Mes)>1){
          Xtilde.abs<-apply(intervals[Mes,],1,testStat,data=data,depth=depth)
          
          
          #weird bug where x tilde comes back as matrix
          if(!is.null(dim(Xtilde.abs))){
            if(dim(Xtilde.abs)[2]>2){
              print("dim ")
              print(dim(Xtilde.abs))
            }
            Xtilde.absT<-Xtilde.abs
            Xtilde.abs<-list()
            for(i in 1:dim(Xtilde.absT)[2])
              Xtilde.abs<-append(Xtilde.abs,list(Xtilde.absT[,i]))
          }
          
          
          bs<-lapply(Xtilde.abs,which.max)
          m0<-which.max(lapply(Xtilde.abs,max))
          b0<-bs[[m0]]+intervals[Mes[m0],1]-1
          
          maxX<-Xtilde.abs[[m0]][bs[[m0]]]
          
        }
        
        
        else if(length(Mes)==1){
          
          Xtilde.abs<-testStat(intervals[Mes,],data=data,depth=depth)
          bs<-which.max(Xtilde.abs)
          m0<-1
          b0<-bs[[m0]]+intervals[Mes[m0],1]-1
          maxX<-max(Xtilde.abs)
          
        }
        
        else{
          
          return(NULL)
        }
        
      }
      
      if(maxX>threshold){
        return(rbind(c(b0,maxX),WBS(intervals,s,b0,threshold,data),WBS(intervals,b0+1,e,threshold,data)))
      }
      
      else
        return(NULL)
    }
    
    #schwartz criteria for choosing the threshold
    #returns all the different models which are the selected change-points and the sigmahatsquared value 
    #(see paper, section on choosing threshold)
    getScwartzCriterias<-function(cp,data,depth){
      
      #get depths
      #make obs univariate
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
      #get indidvidual criteria for a set of cp
      getsSic<-function(cp,data,depth){
        
        depths<-rank(getDepths(data,depth),ties.method = "random")
        # depths<-getDepths(data,depth)
        
        #if at least 1 cp
        if(length(cp)>=1){
          
          
          
          
          N<-length(depths)
          indicies<-cbind(c(1,sort(cp)),c(sort(cp),N+1))
          
          getGroups<-function(vec,dat){return(dat[vec[1]:(vec[2]-1)])}
          
          breaks<-apply(indicies,1,getGroups,dat=depths)
          
          absSum<-lapply(breaks,function(x){sum((x-mean(x))^2)})
          sighatSq<-sum(unlist(absSum))/N
          
          
          return(sighatSq)
        }
        else{
          
          
          N<-length(depths)
          sighatSq<- mean((depths-mean(depths))^2)
          return(sighatSq)
          
          
        }
        
      }
      
      
      sSic<-getsSic(NULL,data,depth)
      
      abc<-list("cp"=NULL,"sigSq"=sSic)
      
      models=list(abc)
      
      for(i in 1:length(cp)){
        
        sSic<-getsSic(cp[1:i],data,depth)
        
        abc$cp=cp[1:i]
        abc$sigSq=sSic
        
        models=append(models,list(abc))
      }
      
      return(models)
      
    }
    
    #sim data
    testData<-simData(N,theta,dist,d)
    
    s.test<-1
    e.test<-nrow(testData)
    #get intervals
    intervals<-getIntervals(1:e.test,numInt)
    #runBS
    cp<-WBS(intervals ,2,e.test,thresh,testData)
    
    ##schwartz modification
    #get the possible models, select one with minimum GSC after
    if(!is.null(cp))
      cp2<-getScwartzCriterias(cp[order(cp[,2],decreasing = T),1],testData,depth)
    
    else
      cp2<-list(list("cp"=NULL,"sigSq"=1))
    
    return(cp2)
  }
  
  
  #parallel running
  #necessary packages
  pkgs<-c("mvtnorm" ,  "fMultivar" ,"MASS"     , "sde"      , "rrcov"  ,   "sn" ,       "ddalpha"  )

  ############################################# HS ############################################# #############################################
  
  

  resultsHS=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,numInt,thresh=thresh,depth="hs")}})
  
  errorhs=inherits(resultsHS, "try-error")
  if(errorhs){
    print("there was an error! HS")
    
  }
  

  
  
  ############################################# Spat ############################################# #############################################
  

  resultsSpat=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,numInt,thresh=thresh,depth="spat")}})
  
  errorsp=inherits(resultsSpat, "try-error")
  if(errorsp){
    print("there was an error! Spat")
    
  }
  
  

  
  
  ############################################# M ############################################# #############################################
  

  resultsMahal=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,numInt,thresh=thresh,depth="mahal")}})
  
  errormh=inherits(resultsMahal, "try-error")
  if(errormh){
    print("there was an error! M")
    
  }
  

  
  ############################################# M75 ############################################# #############################################

  
  resultsMahal75=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,numInt,thresh=thresh,depth="mahal75")}})
  
  errormh75=inherits(resultsMahal75, "try-error")
  if(errormh75){
    print("there was an error! M75")
    
  }
  

  
  fileName<-paste0(N,"_",length(theta),"_",distr,"_",d,"_",numInt,"_","thresh",thresh,"_WBS_Ranks_SS_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  
  save(resultsHS,resultsSpat,resultsMahal,resultsMahal75,file=fileName1)
  print(i/numUniqueRuns)

}


stopCluster(cl)
closeAllConnections()

# 

