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
# Ns=100
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






#######################NOW RUN PELT 

# 
# 
# 
Rcpp::sourceCpp('/scratch/k3ramsay/Codes MKWC/PELT_CPP.cpp')
# 
# 
# constants<-c(0.001,0.01,0.1,1,3,5)
constants<-c(seq(0.03,0.15,l=13))
constants=round(constants,2); constants
constants
pkgs<-c("mvtnorm" ,  "fMultivar" ,"MASS"     , "sde"      , "rrcov"  ,   "sn" ,       "ddalpha"  )
# no_cores<-detectCores()-1

no_cores<-Sys.getenv("SLURM_CPUS_PER_TASK") 

cl <-  makeCluster(no_cores,type="FORK")

registerDoParallel(cl) 
registerDoRNG(seed = 440)

#   
for(j in 1:numUniqueRuns){
  params<-paramterIndices[j,]
  N=Ns[params[1]]
  # beta=constant*sqrt(N)+3.74
  
  theta=thetas[[params[2]]]
  distr=distributions1[[params[3]]]
  d=ds[params[4]]
  
  dName=names(distributions1)[params[3]]
  distr=distributions[[distr]]
  
  
  for(constant in constants){
    fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_new_constant_2_C1_",constant,"OC_PELT_ranks_simsize_",sim.size,sep="")
    fileName1<-paste(dirr2,fileName,".Rda",sep="")
    if(!file.exists(fileName1)){
      beta=constant*log(N)*sqrt(N)
  
  
      fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_","OC_Depths_simsize_",sim.size,sep="")
  
  
      fileName1<-paste(dirr,fileName,".Rda",sep="")
      load(fileName1)
  
  
  
      ranksHS<-lapply(depthsHS,rank)
      ranksSpat<-lapply(depthsSpat,rank)
      ranksMahal<-lapply(depthsMahal,rank)
      ranksMahal75<-lapply(depthsMahal75,rank)
  
  
      # Rank_PELT(ranksHS[[1]],beta=10)
      # beta=15
      resultsHS<-mclapply(ranksHS, PELT_T,beta=beta,mc.cores=no_cores)
      resultsSpat<-mclapply(ranksSpat, PELT_T,beta=beta,mc.cores=no_cores)
      resultsMahal<-mclapply(ranksMahal, PELT_T,beta=beta,mc.cores=no_cores)
      resultsMahal75<-mclapply(ranksMahal75, PELT_T,beta=beta,mc.cores=no_cores)
  
      fileName<-paste0(N,"_",length(theta),"_",dName,"_",d,"_new_constant_2_C1_",constant,"OC_PELT_ranks_simsize_",sim.size,sep="")
      fileName1<-paste(dirr2,fileName,".Rda",sep="")
      save(resultsHS,resultsSpat,resultsMahal,resultsMahal75,file=fileName1)
    }
    print(j/numUniqueRuns+(which(constant==constants)-1)/length(constants)/numUniqueRuns)

  }
}

stopCluster(cl)
closeAllConnections()
# 
# 
# 


