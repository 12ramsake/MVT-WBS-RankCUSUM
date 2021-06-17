
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

  # FLAG=0
  red=data[(s+1):e]
  p=ncol(data)
  n=nrow(data)
  n_se=(e-s)
  
  if(n_se>(2*p*log(n)+1)){
    
    t_range=ceiling(s+p*log(n)):floor(e-p*log(n))
    set_op_norms=sapply(t_range,function(t){norm(St(t,s,e,data),type="2")} )
    a=max(set_op_norms)
    
    
    if(a<=tau){
      # FLAG=1
      return(NULL)
    }
    else{
      b=t_range[which.max(set_op_norms)]
      cps=rbind(cps,c(b,a))
      # cps=c(cps, BSOP(s,b-1,data,tau,cps))
      # cps=c(cps, BSOP(b-1,e,data,tau,cps))
      return(rbind(cps, BSOP(s,b-1,data,tau), BSOP(b-2,e,data,tau)))
    }
  }
  else
    return(NULL)
  
}




data=rbind(replicate(2,rnorm(60)),replicate(2,rnorm(60,0,10)))
p=2
n=120
# BSOP(p*log(n),n-p*log(n),data,50*sqrt(n)^.9)
# unique()

BSOP(p*log(n),n-p*log(n),data,300)






binary.segmentation(data_M,alpha=.05,power_enhancement=TRUE,M_threshold=0.05)



