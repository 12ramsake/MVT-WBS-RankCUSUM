#Seed
set.seed(746383)
#packages
require(mvtnorm)
require(fMultivar)
require(MASS)
require(sde)
require(rrcov)
require(xtable)
require(sn)
require(ddalpha)
require(doParallel)
##wild binary segmentation
#calculate the test statistic
#spat, hs, mahal,mahal75
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

getStatFromDepths<-function(depths,N){
  ranks<-rank(depths,ties.method = "random")
  expected.val<-(N+1)/2
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


#returns indices of the intervals
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

checkIfSubInterval<-function(sub,super){
  return(sub[1]>=super[1]&&sub[2]<=super[2])
}

#let the set of intervals be a matrix with 2 columns

WBS<-function(intervals,s,e,threshold,data,depth="mahal"){
  
  
  # sig.level=sig.level/2
  # threshold=qBB(1-sig.level)$root
  
  if((e-s)<1)
    return(NULL)
  
  else{
    #intervals contained in s,e
    Mes<-which(apply(intervals,1,checkIfSubInterval,super=c(s,e)))
    
    
    if(length(Mes)>1){
      Xtilde.abs<-apply(intervals[Mes,],1,testStat,data=data,depth=depth)
      
      # print("Xtilde ")
      # print(dim(Xtilde.abs))
      
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
      
      # print("intervals ")
      # print(intervals)
      # print("Mes ")
      # print(Mes)
      # print("m0 ")
      # print(m0)
      
      
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
    # sig.level=sig.level/2
    return(rbind(c(b0,maxX),WBS(intervals,s,b0,threshold,data),WBS(intervals,b0+1,e,threshold,data)))
  }
  
  else
    return(NULL)
}

#schwartz criteria for choosing the threshold
getScwartzCriterias<-function(cp,data,depth){
  
  
  cp<-cp[order(cp[,2],decreasing = T),]
  
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
      # absSum<-lapply(breaks,function(x){abs(x-median(x))})
      
      
      # sighatSq<-mean(unlist(absSum),trim=.125)
      
      
      return(sighatSq)
      # sSic<-(N/2)*log(sighatSq)+length(cp)*(log(N))^alpha
    }
    else{
      
      
      N<-length(depths)
      sighatSq<- mean((depths-mean(depths))^2)
      # sighatSq<- mean(abs(depths-median(depths)),trim = .125)
      return(sighatSq)
      # sSic<-(N/2)*log(sighatSq)
      
    }
    
    # return(sSic)
  }
  
  #get ssic for all ammounts of cp
  
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
  
  # minVal<-which.min(sSic)
  # 
  # if(minVal==1)
  #   return(NULL)
  # else
  #   return(cp[1:(minVal-1)])
  
}

applySCH<-function(models,N,alpha){
  #get indidvidual criteria for a set of cp
  getsSic<-function(cp,sighatSq){
    
    
    #if at least 1 cp
    if(length(cp)>=1){
      
      sSic<-(N/2)*log(sighatSq)+length(cp)*(log(N))^alpha
    }
    else{
      sSic<-(N/2)*log(sighatSq)
    }
    
    return(sSic)
  }
  
  #get ssic for all ammounts of cp
  
  sSic<-getsSic(models[[1]]$cp,models[[1]]$sigSq)
  
  for(j in 2:length(models)){
    
    sSic<-c(sSic,getsSic(models[[j]]$cp,models[[j]]$sigSq))
    
  }
  
  minVal<-which.min(sSic)
  
  return(models[[minVal]]$cp)  
  
}
