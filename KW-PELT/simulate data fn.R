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
