#You need Rcpp and the file PELT_CPP.cpp in your working directory
Rcpp::sourceCpp('PELT_CPP.cpp')

# dat is the data
# beta is the penalty
# depth is the depth function
PELT_R=function(dat,beta,depth='spat'){
  
  if(depth=='hs')
    depths<-ddalpha::depth.halfspace(dat,dat)
  if(depth=='spat')
    depths<-ddalpha::depth.spatial(dat,dat)
  if(depth=='m')
    depths<-ddalpha::depth.Mahalanobis(dat,dat)
  if(depth=='m75')
    depths<-ddalpha::depth.Mahalanobis(dat,dat,"MCD")
  
  ranks=rank(depths,ties.method = 'random')
  return(PELT_T(ranks,beta))
}

#Example run 
set.seed(440)
test_data=rbind(replicate(2,rnorm(200)),replicate(2,rnorm(200,5)),replicate(2,rnorm(200,0.2)))

beta_default=3.74+0.18*sqrt(nrow(test_data))
PELT_R(test_data,beta_default,depth='spat')
