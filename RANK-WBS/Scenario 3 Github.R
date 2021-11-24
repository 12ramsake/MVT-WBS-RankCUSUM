
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



dirr<- "/u/k3ramsay/ResearchDocuments/output/WBS_SIMULATION/"
setwd(dirr)


#simulation parameters

#N
Ns <- c(1000)

#simulation size, num repetitions
sim.size = 100

#0,2,3.5 cps
thetas <- list(c(.333, .666),
               c(.25, .5, .75),
               c(1 / 6, 2 / 6, 3 / 6, 4 / 6, 5 / 6))


thresh=1


d2s = c(4, 3, 2, 1)

ds=5
##Create Parameter Vector

numUniqueRuns <- length(thetas) * length(d2s)

paramterIndices <- matrix(0, ncol = 4, nrow = numUniqueRuns)
curr = 0
for (i1 in 1:length(Ns)) {
  for (i2 in 1:length(thetas)) {
    for (i3 in 1:length(d2s)) {
      for (i4 in 1:length(ds)) {
        curr = curr + 1
        paramterIndices[curr, ] = c(i1, i2, i3, i4)
      }
    }
  }
}







#simulation functions

#generates iid d-normal, sigmaxI Rv
normalMaster <-
  function(n, d, sigmaSq, d2) {
    mvtnorm::rmvnorm(n, sigma = diag(c(sigmaSq * rep(1, d2), rep(1, d - d2))))
  }

norm1 <- function(n, d, d2) {
  normalMaster(n, d, 1, d2)
}
norm2 <- function(n, d, d2) {
  normalMaster(n, d, 2.5, d2)
}
norm3 <- function(n, d, d2) {
  normalMaster(n, d, 4, d2)
}
norm4 <- function(n, d, d2) {
  normalMaster(n, d, 2.25, d2)
}
norm5 <- function(n, d, d2) {
  normalMaster(n, d, 5, d2)
}

normals <- list(norm1, norm2, norm3, norm4, norm5, norm1)



for (i in 1:numUniqueRuns) {
  params <- paramterIndices[i, ]
  N = 1000
  theta = thetas[[params[2]]]
  d2 = d2s[params[3]]
  d = 5
  numInt = floor(log(N)) *100
  dName = "Normal"
  distr = normals
  #repeat dirr for parallel
  
  
  #run one repetition
  runSim <-function(N,
                    theta,
                    rdata,
                    d,
                    numInt = 10,
                    thresh = 1.3584,
                    depth = "hs",
                    d2) {
    #simulate data set for one repitition
    simData <- function(N, theta, rdata, d, d2) {
      #cp locations
      if (!is.null(theta))
        locations <- c(1, floor(theta * N), N)
      else
        locations <- N
      
      #data
      dat <- matrix(0, nrow = N, ncol = d)
      
      for (i in 2:(length(theta) + 2))
        dat[locations[i - 1]:locations[i], ] <-
        rdata[[i - 1]](length(locations[i - 1]:locations[i]), d, d2)
      
      return(dat)
    }
    
    
    ##wild binary segmentation
    #calculate the test statistic
    #spat, hs, mahal,mahal75
    testStat <- function(range, data, depth) {
      if (depth == "spat") {
        ts = testStatSpat(range, data)
      }
      else if (depth == "hs") {
        ts = testStatHs(range, data)
      }
      else if (depth == "mahal") {
        ts = testStatMahal(range, data)
      }
      else if (depth == "mahal75") {
        ts = testStatMahal75(range, data)
      }
      else{
        ts = NULL
        print("bad depth specification")
      }
      
      return(ts)
    }
    
    getStatFromDepths <- function(depths, N) {
      ranks <- rank(depths, ties.method = "random")
      expected.val <- (N + 1) / 2
      std.dev <- sqrt((N ^ 2 - 1) / 12)
      cusum <- cumsum(N ^ (-0.5) * (ranks - expected.val) / std.dev)
      return(abs(cusum)[1:(length(cusum) - 1)])
    }
    
    testStatHs <- function(range, data) {
      if ((range[2] - range[1]) > (ncol(data) + 1)) {
        range <- range[1]:range[2]
        N <- nrow(data[range, ])
        depths <- depth.halfspace(data[range, ], data[range, ])
        return(getStatFromDepths(depths, N))
      }
      else
        return(-10)
    }
    
    testStatSpat <- function(range, data) {
      if ((range[2] - range[1]) > (ncol(data) + 1)) {
        range <- range[1]:range[2]
        N <- nrow(data[range, ])
        depths <- depth.spatial(data[range, ], data[range, ])
        
        return(getStatFromDepths(depths, N))
      }
      else
        return(-10)
    }
    
    testStatMahal75 <- function(range, data) {
      if ((range[2] - range[1]) > (ncol(data) * 2)) {
        range <- range[1]:range[2]
        
        N <- nrow(data[range, ])
        depths <- depth.Mahalanobis(data[range, ], data[range, ], "MCD")
        
        return(getStatFromDepths(depths, N))
      }
      else
        return(-10)
    }
    
    testStatMahal <- function(range, data) {
      if ((range[2] - range[1]) > (ncol(data) * 2)) {
        range <- range[1]:range[2]
        
        N <- nrow(data[range, ])
        depths <- depth.Mahalanobis(data[range, ], data[range, ])
        return(getStatFromDepths(depths, N))
      }
      else
        return(-10)
    }
    
    
    #returns indices of the intervals
    getIntervals <- function(indices, M) {
      ints <- t(replicate(M, sort(sample(indices, 2))))
      diffs <- (ints[, 2] - ints[, 1]) == 1
      if (any(diffs)) {
        ints[diffs, ] = getIntervals(indices, sum(diffs))
        return(ints)
      }
      else{
        return(ints)
      }
    }
    
    checkIfSubInterval <- function(sub, super) {
      return(sub[1] >= super[1] && sub[2] <= super[2])
    }
    
    #let the set of intervals be a matrix with 2 columns
    
    WBS <- function(intervals, s, e, threshold, data) {
      # sig.level=sig.level/2
      # threshold=qBB(1-sig.level)$root
      
      if ((e - s) < 1)
        return(NULL)
      
      else{
        #intervals contained in s,e
        Mes <-
          which(apply(intervals, 1, checkIfSubInterval, super = c(s, e)))
        
        
        if (length(Mes) > 1) {
          Xtilde.abs <- apply(intervals[Mes, ],
                              1,
                              testStat,
                              data = data,
                              depth = depth)
          
          # print("Xtilde ")
          # print(dim(Xtilde.abs))
          
          #weird bug where x tilde comes back as matrix
          if (!is.null(dim(Xtilde.abs))) {
            if (dim(Xtilde.abs)[2] > 2) {
              print("dim ")
              print(dim(Xtilde.abs))
            }
            Xtilde.absT <- Xtilde.abs
            Xtilde.abs <- list()
            for (i in 1:dim(Xtilde.absT)[2])
              Xtilde.abs <- append(Xtilde.abs, list(Xtilde.absT[, i]))
          }
          
          
          bs <- lapply(Xtilde.abs, which.max)
          m0 <- which.max(lapply(Xtilde.abs, max))
          
          # print("intervals ")
          # print(intervals)
          # print("Mes ")
          # print(Mes)
          # print("m0 ")
          # print(m0)
          
          
          b0 <- bs[[m0]] + intervals[Mes[m0], 1] - 1
          
          maxX <- Xtilde.abs[[m0]][bs[[m0]]]
          
        }
        
        
        else if (length(Mes) == 1) {
          Xtilde.abs <- testStat(intervals[Mes, ], data = data, depth = depth)
          bs <- which.max(Xtilde.abs)
          m0 <- 1
          b0 <- bs[[m0]] + intervals[Mes[m0], 1] - 1
          maxX <- max(Xtilde.abs)
          
        }
        
        else{
          return(NULL)
        }
        
      }
      
      if (maxX > threshold) {
        # sig.level=sig.level/2
        return(rbind(
          c(b0, maxX),
          WBS(intervals, s, b0, threshold, data),
          WBS(intervals, b0 + 1, e, threshold, data)
        ))
      }
      
      else
        return(NULL)
    }
    
    #schwartz criteria for choosing the threshold
    getScwartzCriterias <- function(cp, data, depth) {
      #get depths
      #make obs univariate
      getDepths <- function(data, depth) {
        if (depth == "spat") {
          ts = depth.spatial(data, data)
        }
        else if (depth == "hs") {
          ts = depth.halfspace(data, data)
        }
        else if (depth == "mahal") {
          ts = depth.Mahalanobis(data, data)
        }
        else if (depth == "mahal75") {
          ts = depth.Mahalanobis(data, data, "MCD")
        }
        else{
          ts = NULL
          print("bad depth specification")
        }
        
        return(ts)
      }
      #get indidvidual criteria for a set of cp
      getsSic <- function(cp, data, depth) {
        depths <- rank(getDepths(data, depth), ties.method = "random")
        # depths<-getDepths(data,depth)
        
        #if at least 1 cp
        if (length(cp) >= 1) {
          N <- length(depths)
          indicies <- cbind(c(1, sort(cp)), c(sort(cp), N + 1))
          
          getGroups <-
            function(vec, dat) {
              return(dat[vec[1]:(vec[2] - 1)])
            }
          
          breaks <- apply(indicies, 1, getGroups, dat = depths)
          
          absSum <- lapply(breaks, function(x) {
            sum((x - mean(x)) ^ 2)
          })
          sighatSq <- sum(unlist(absSum)) / N
          # absSum<-lapply(breaks,function(x){abs(x-median(x))})
          
          
          # sighatSq<-mean(unlist(absSum),trim=.125)
          
          
          return(sighatSq)
          # sSic<-(N/2)*log(sighatSq)+length(cp)*(log(N))^alpha
        }
        else{
          N <- length(depths)
          sighatSq <- mean((depths - mean(depths)) ^ 2)
          # sighatSq<- mean(abs(depths-median(depths)),trim = .125)
          return(sighatSq)
          # sSic<-(N/2)*log(sighatSq)
          
        }
        
        # return(sSic)
      }
      
      #get ssic for all ammounts of cp
      
      sSic <- getsSic(NULL, data, depth)
      
      abc <- list("cp" = NULL, "sigSq" = sSic)
      
      models = list(abc)
      
      for (i in 1:length(cp)) {
        sSic <- getsSic(cp[1:i], data, depth)
        
        abc$cp = cp[1:i]
        abc$sigSq = sSic
        
        models = append(models, list(abc))
      }
      
      return(models)
      
      # minVal<-which.min(sSic)
      #
      # if(minVal==1)
      #   return(NULL)
      # else
      #   return(cp[1:(minVal-1)])
      
    }
    
    testData <- simData(N, theta, rdata, d, d2)
    
    s.test <- 1
    e.test <- nrow(testData)
    
    intervals <- getIntervals(1:e.test, numInt)
    
    cp <- WBS(intervals , 1, e.test, thresh, testData)
    
    ##schwartz modification
    
    if (!is.null(cp))
      cp2 <-
      getScwartzCriterias(cp[order(cp[, 2], decreasing = T), 1], testData, depth)
    
    else
      cp2 <- list(list("cp" = NULL, "sigSq" = 1))
    
    return(cp2)
  }
  
  pkgs <-
    c("mvtnorm" ,
      "fMultivar" ,
      "MASS"     ,
      "sde"      ,
      "rrcov"  ,
      "sn" ,
      "ddalpha")
  
  no_cores <- detectCores() - 1
  #  no_cores<-100
  
  
  cl <- makeCluster(no_cores)
  
  registerDoParallel(cl)
  clusterExport(cl = cl, c("normalMaster"))
  registerDoRNG(seed = 440)
  
  
  resultsHS = try({
    foreach(i = 1:sim.size, .packages = pkgs) %dopar% {runSim(N,theta,distr,d,numInt,thresh = thresh,depth = "hs",d2)}
  })
  
  errorhs = inherits(resultsHS, "try-error")
  if (errorhs) {
    print("there was an error! HS")
    
  }
  
  registerDoSEQ()
  stopCluster(cl)
  closeAllConnections()
  
  ################################################################################################
  
  cl <- makeCluster(no_cores, type = "FORK", outfile = "debug.txt")
  registerDoParallel(cl)
  
  registerDoRNG(seed = 440)
  resultsSpat = try({
    foreach(i = 1:sim.size, .packages = pkgs) %dopar% {runSim(N,theta,distr,d,numInt,thresh = thresh,depth = "spat",d2)}
  })
  
  errorsp = inherits(resultsSpat, "try-error")
  if (errorsp) {
    print("there was an error! Spat")
    
  }
  
  
  registerDoSEQ()
  stopCluster(cl)
  closeAllConnections()
  
  
  ###########################################################################################
    cl <- makeCluster(no_cores)
    
    registerDoParallel(cl)
    registerDoRNG(seed = 440)
    clusterExport(cl = cl, c("normalMaster"))
    
    
    resultsMahal = try({foreach(i = 1:sim.size, .packages = pkgs) %dopar% {runSim(N,theta,distr,d,numInt,thresh = thresh,depth = "mahal",d2)}})
    
    errormh = inherits(resultsMahal, "try-error")
    if (errormh) {
      print("there was an error! M")
      
    }
    
    registerDoSEQ()
    stopCluster(cl)
    closeAllConnections()

  
    ######################################################################################################################################################################################
    cl <- makeCluster(no_cores)
    
    registerDoParallel(cl)
    registerDoRNG(seed = 440)
    clusterExport(cl = cl, c("normalMaster"))
    
    resultsMahal75 = try({foreach(i = 1:sim.size, .packages = pkgs) %dopar% {runSim(N,theta,distr,d,numInt,thresh = thresh,depth = "mahal75",d2)}})
  
    errormh75 = inherits(resultsMahal75, "try-error")
    if (errormh75) {
      print("there was an error! M75")
      
    }
    
    registerDoSEQ()
    stopCluster(cl)
    closeAllConnections()

  fileName <-paste0(N,"_",length(theta),"_",dName,"_",d,"_",numInt,"_", "thresh", thresh,"_d2_", d2, "_WBS_simsize_", sim.size, "_all", sep = "")
  
  
  fileName1 <- paste(dirr, fileName, ".Rda", sep = "")
  
  save(resultsHS, resultsSpat, resultsMahal, resultsMahal75, file = fileName1)
  
  closeAllConnections()
  print(i / numUniqueRuns)
}

