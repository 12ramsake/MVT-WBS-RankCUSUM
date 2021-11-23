

##change-point simulation for wang2021 paper

library(mvtnorm)
library(fMultivar)
library(MASS)
library(sde)
library(rrcov)
library(xtable)
library(sn)
library(doParallel)
library(doRNG)

#results dirrectory
dirr<- ""
setwd(dirr)


#simulation parameters
#N
Ns<-c(1000,2500,5000)

#simulation size, num repetitions
sim.size=100
mod=100


#xhange below
#0,2,3.5 cps
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




#run simulation
for(i in 1:108){
  
  #set paramters
  params<-paramterIndices[i,]
  N=Ns[params[1]]
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  distr=distributions[[distr]]
  d=ds[params[4]]
  numInt=floor(log(Ns[params[1]]))*mod
  dName=names(distributions1)[params[3]]
  
  
  
  
  #run one repetition
  runSim<-function(N,theta,rdata,d,numInt=10){
    
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
    
    
    #functions for BSOP
    
    
    St=function(t,s,e,data){
      # print(paste0("t ",t))
      # print(paste0("s ",s))
      # print(paste0("e ",e))
      if(is.null(dim(data)))
        data=matrix(data,ncol=1)
      
      term1=sqrt((e-t)/((e-s)*(t-s)))*(t(data[(s+1):t,])%*%(data[(s+1):t,]))
      term2=sqrt((t-s)/((e-t)*(e-s)))*(t(data[(t+1):e,])%*%data[(t+1):e,])
      return(term1-term2)
    }
    
    
    
    BSOP=function(s,e,data,tau=1,cps=NULL){
      
      FLAG=0
      red=data[(s+1):e]
      p=ncol(data)
      n=nrow(data)
      n_se=(e-s)
      
      if(n_se>(2*p*log(n)+1) & FLAG==0){
        t_range=ceiling(s+p*log(n)):floor(e-p*log(n))
        set_op_norms=sapply(t_range,function(t){norm(St(t,s,e,data),type="2")} )
        a=max(set_op_norms)
        if(a<=tau){
          FLAG=1
          return(NULL)
        }
        else{
          b=t_range[which.max(set_op_norms)]
          cps=rbind(cps,c(b,a))
          # cps=c(cps, BSOP(s,b-1,data,tau,cps))
          # cps=c(cps, BSOP(b-1,e,data,tau,cps))
          return(rbind(cps, BSOP(s,b-1,data,tau,cps), BSOP(b-2,e,data,tau,cps)))
        }
      }
      else
        return(NULL)
      
    }
    
    
    
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
    
    
    PC1=function(Xt,alpha,beta){
      p=ncol(Xt)
      n=nrow(Xt)
      if((beta-alpha)>2*p*log(n)+1){
        
        t_vals=ceiling(alpha+p*log(n)):floor(beta-p*log(n))
        dm=t_vals[which.max(sapply(t_vals,function(t){norm(St(t,alpha,beta,Xt),type="2")} ))]
        um=eigen(St(dm,alpha,beta,Xt))$vectors[,1]
        
      }
      else{um=rep(0,p)}
      return(um)
    }
    
    PC=function(Wt,intervals){
      M=nrow(intervals)
      ums=NULL
      for(i in 1:M){
        ums=rbind(ums,PC1(Wt,intervals[i,1],intervals[i,2]))
      }
      return(ums)
    }
    #let the set of intervals be a matrix with 2 columns
    
    WBSIP<-function(data,s,e,intervals,tau){
      
      
      # sig.level=sig.level/2
      # threshold=qBB(1-sig.level)$root
      p=ncol(data)
      n=nrow(data)
      Wt=data[seq(1,nrow(data),by=2),]
      Xt=data[seq(2,nrow(data),by=2),]
      M=nrow(intervals)
      #u has M rows
      u=PC(Wt,intervals)
      # 
      # s=floor(s/2)
      # e=floor(e/2)
      # intervals2=floor(intervals)/2
      #M by n
      Ytu=u%*%t(Xt)
      Ytu=Ytu^2
      
      if((e-s)<(2*p*log(n/2)+1))
        return(NULL)
      else{
        #intervals contained in s,e
        # Mes<-which(apply(intervals2,1,checkIfSubInterval,super=c(s,e)))
        left_endpoint=sapply(intervals[,1],function(x){max(x,s)})
        right_endpoint=sapply(intervals[,2],function(x){min(x,e)})
        Mes=which((right_endpoint-left_endpoint)>=(2*log(n/2)+1))
        
        if(length(Mes)>1){
          am=rep(-1,M)
          bm=rep(-1,M)
          for(j in Mes){
            t_vals=ceiling(left_endpoint[j]+log(n/2)):floor(right_endpoint[j]-log(n/2))
            
            candidate_ys<-sapply(t_vals,function(t){abs(St(t,left_endpoint[j],right_endpoint[j],Ytu[j,]))} )
            
            mm=which.max(candidate_ys)
            
            bm[j]=t_vals[mm[1]]
            am[j]=candidate_ys[mm[1]]
          }  
          
          m=which.max(am)
          
          
          
          if(am[m[1]]>tau){
            # sig.level=sig.level/2
            return(rbind(c(bm[m[1]]*2,am[m[1]]),
                         WBSIP(data,s,bm[m[1]],intervals,tau),
                         WBSIP(data,bm[m[1]]+1,e,intervals,tau)))
          }
          else
            return(NULL)
        }
        else
          return(NULL)
      }
    }
    
    
    

    
    
    # unique(BSOP(p*log(n),n-p*log(n),data,50*sqrt(n)^.9))
    p=d
    n=N
    #sim data
    testData<-simData(N,theta,rdata,d)
    
    s.test<-1
    e.test<-N
    #get intervals
    # intervals<-getIntervals(1:e.test,numInt)
    intervals=getIntervals(0:floor(n/2),numInt)
    big_enough=function(i){(i[2]-i[1])>(2*p*log(n/2)+1)}
    intervals=intervals[apply(intervals, 1, big_enough),]
    count=1
    while(nrow(intervals)<numInt&&count<100){
      intervals2=getIntervals(0:floor(n/2),numInt)
      intervals=rbind(intervals,intervals2[apply(intervals2, 1, big_enough),])
      count=count+1
    }
    intervals=intervals[1:numInt,]
    
    
    #runBS
    # WBSIP(testData,p*log(n/2)+1,n/2-p*log(n/2)+1,intervals,1)
    cp_WBSIP<-WBSIP(testData,p*log(n/2)+1,n/2-p*log(n/2)+1,intervals,1)
    cp_BSOP<-BSOP(p*log(n),n-p*log(n),testData,1)
    
    return(list(WBSIP=cp_WBSIP,BSOP=cp_BSOP))
  }
  
  
  #parallel running
  #necessary packages
  pkgs<-c("mvtnorm" ,  "fMultivar" ,"MASS"     , "sde"      , "rrcov"  ,   "sn"   )
  
  
  no_cores<-detectCores()-1
  
  
  cl <- makeCluster(no_cores,type="FORK")  
  
  registerDoParallel(cl) 
  
  #for windows
  #clusterExport(cl=cl,c("normalMaster","cauchyMaster","skewNormalMaster"))
  
  registerDoRNG(seed = 440)
  
  
  results=try({foreach(i=1:sim.size,.packages= pkgs) %dopar% {runSim(N,theta,distr,d,numInt)}})
  
  errorhs=inherits(results, "try-error")
  if(errorhs){
    print("there was an error!")
    
  }
  
  registerDoSEQ()
  stopCluster(cl)
  closeAllConnections()
  
  
  
  
  fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_WANG_simsize_",sim.size,sep="")
  fileName1<-paste(dirr,fileName,".Rda",sep="")
  
  save(results,file=fileName1)
  
  closeAllConnections()
  print(i/36)
}








