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
# 
# checkIfSubInterval<-function(sub,super){
#   return(sub[1]>=super[1]&&sub[2]<=super[2])
# }






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




intervals=getIntervals(0:floor(n/2),100)
#check if the interval is big enough

big_enough=function(i){(i[2]-i[1])>(2*p*log(n/2)+1)}
intervals=intervals[apply(intervals, 1, big_enough),]


p=2
n=120
kappa=norm(diag(rep(9,2)),type="2")
B=5
Delta=0.075*n
tau=Delta*kappa*n^.4

data=rbind(replicate(p,rnorm(n/2)),replicate(p,rnorm(n/2,0,10)))

WBSIP(data,p*log(n/2)+1,n/2-p*log(n/2)+1,intervals,20)

