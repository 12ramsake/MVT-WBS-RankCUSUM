#packages
# library(fMultivar)
# library(MASS)
# library(depth)
# library(rrcov)
# library(mrfDepth)
# library(ddalpha)


##wild binary segmentation
#calculate the test statistic
#spat, hs, mahal,mahal75 are the depth parameters
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

#test cusum from depth values
getStatFromDepths <- function(depths, N) {
  ranks <- rank(depths, ties.method = "random")
  expected.val <- (N - 1) / 2
  std.dev <- sqrt((N ^ 2 - 1) / 12)
  cusum <- cumsum(N ^ (-0.5) * (ranks - expected.val) / std.dev)
  return(abs(cusum)[1:(length(cusum) - 1)])
}

testStatHs <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) + 1)) {
    range <- range[1]:range[2]
    N <- nrow(data[range, ])
    depths <- ddalpha::depth.halfspace(data[range, ], data[range, ])
    return(getStatFromDepths(depths, N))
  }
  else
    return(-10)
}

testStatSpat <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) + 1)) {
    range <- range[1]:range[2]
    N <- nrow(data[range, ])
    depths <- ddalpha::depth.spatial(data[range, ], data[range, ])
    
    return(getStatFromDepths(depths, N))
  }
  else
    return(-10)
}

testStatMahal75 <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) * 2)) {
    range <- range[1]:range[2]
    
    N <- nrow(data[range, ])
    depths <- ddalpha::depth.Mahalanobis(data[range, ], data[range, ], "MCD")
    
    return(getStatFromDepths(depths, N))
  }
  else
    return(-10)
}

testStatMahal <- function(range, data) {
  if ((range[2] - range[1]) > (ncol(data) * 2)) {
    range <- range[1]:range[2]
    
    N <- nrow(data[range, ])
    depths <- ddalpha::depth.Mahalanobis(data[range, ], data[range, ])
    return(getStatFromDepths(depths, N))
  }
  else
    return(-10)
}


#returns indices of the intervals selected, M is the number of intervals
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

#checks if an interval is a sub of another
checkIfSubInterval <- function(sub, super) {
  return(sub[1] >= super[1] && sub[2] <= super[2])
}



#let the set of intervals be a matrix with 2 columns

WBS <- function(intervals, s, e, threshold, data,depth) {
  if ((e - s) < 1)
    return(NULL)
  
  else{
    #intervals contained in s,e
    Mes <- which(apply(intervals, 1, checkIfSubInterval, super = c(s, e)))
    
    
    if (length(Mes) > 1) {
      Xtilde.abs <- apply(intervals[Mes, ],
                          1,
                          testStat,
                          data = data,
                          depth = depth)
      
      
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
    return(rbind(
      c(b0, maxX),
      WBS(intervals, s, b0, threshold, data,depth),
      WBS(intervals, b0 + 1, e, threshold, data,depth)
    ))
  }
  
  else
    return(NULL)
}

#get depths
#make obs univariate
getDepths <- function(data, depth) {
  if (depth == "spat") {
    ts = ddalpha::depth.spatial(data, data)
  }
  else if (depth == "hs") {
    ts = ddalpha::depth.halfspace(data, data)
  }
  else if (depth == "mahal") {
    ts = ddalpha::depth.Mahalanobis(data, data)
  }
  else if (depth == "mahal75") {
    ts = ddalpha::depth.Mahalanobis(data, data, "MCD")
  }
  else{
    ts = NULL
    print("bad depth specification")
  }
  
  return(ts)
}

#schwartz criteria for choosing the threshold
#returns all the different models which are the selected change-points and the sigmahatsquared value
#(see paper, section on choosing threshold)
getScwartzCriterias <- function(cp, data, depth) {

  #get indidvidual criteria for a set of cp
  
  
  
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
  
}

getsSic <- function(cp, data, depth) {
  depths <- rank(getDepths(data, depth), ties.method = "random")
  # depths<-getDepths(data,depth)
  
  #if at least 1 cp
  if (length(cp) >= 1) {
    N <- length(depths)
    indicies <- cbind(c(1, sort(cp)), c(sort(cp), N + 1))
    
    getGroups <- function(vec, dat) {
      return(dat[vec[1]:(vec[2] - 1)])
    }
    
    breaks <- apply(indicies, 1, getGroups, dat = depths)
    
    absSum <- lapply(breaks, function(x) {
      sum((x - mean(x)) ^ 2)
    })
    sighatSq <- sum(unlist(absSum)) / N
    
    
    return(sighatSq)
  }
  else{
    N <- length(depths)
    sighatSq <- mean((depths - mean(depths)) ^ 2)
    return(sighatSq)
    
    
  }
  
}

DWBS = function(data,
                    N,
                    d,
                    numInt = 10,
                    thresh = 1.3584,
                    depth = "hs") {
  s_test <- 1
  e_test <- nrow(data)
  
  #get intervals
  intervals <- getIntervals(1:e_test, numInt)
  
  #runWBS
  cp <- WBS(intervals , 2, e_test, thresh, data,depth)
  
  return(cp)
}


DWBS_DDT = function(data,
                    N,
                    d,
                    numInt = 10,
                    thresh = 1.3584,
                    depth = "hs") {
  s_test <- 1
  e_test <- nrow(data)
  
  #get intervals
  intervals <- getIntervals(1:e_test, numInt)
  
  #runWBS
  cp <- WBS(intervals , 2, e_test, thresh, data,depth)
  
  ## Schwartz modification
  # get the possible models, select one with minimum GSC after
  # This computes the sigma sq hats
  if (!is.null(cp))
    cp2 <-
    getScwartzCriterias(cp[order(cp[, 2], decreasing = T), 1], data, depth)
  
  else
    cp2 <- list(list("cp" = NULL, "sigSq" = 1))
  
  return(cp2)
}

set.seed(440)
test_data=rbind(replicate(2,rnorm(100)),replicate(2,rnorm(100,5)))
thresh=2

DWBS(test_data, N, d, numInt=100, thresh = thresh, depth = "spat")
DWBS_DDT(test_data, N, d, numInt=100, thresh = thresh, depth = "spat")






