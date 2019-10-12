#' Search for the optimal number of clusters
#'
#' Binary search for the optimal number of clusters for a two-arm CRT using simulations
#'
#' @param nrep the number of replications of the simulation procedure for generating the pseudo study data, default is \code{"1000"})
#' @param nt the initial number of cluster in the treatment arm to start the search
#' @param nc the initial number of cluster in the control arm to start the search
#' @param tcRatio the ratio of number of clusters in the intervention arm vs that in the control arm, default value is \code{"1"} for balanced sample size allocation 
#'                
#' @param minpower the minimum power, default at 0.8
#' @param alpha type-1 errer, default at 0.05 
#' @param increaseSamplingBy the rate to increas the number of replications (nrep) when power gets close to the minpower, default at 1 (not increase)
#' @param PermutationTest indicator for if Fisher's sharp null hypothesis is test and permutation test power is of interest, default at FALSE (No)
#' @param Npermutationtest=100 if a permutation test is of interest, the number of replications when estimating the permutation test power
#' @inheritParams simPower
#' @inheritParams simPowerPT
#' 
#' @return a list of (nt, nc, power) for the optimal design
#'
#' @examples
#'\dontrun{
#' ## distribution of cluster sizes
#' sim_cluster_size = function(N, ...){
#'   size = round(100*rnorm(N, 0, 1), 0)
#'   size[size<=0] = 1
#'   return(size)
#' }
#' ## distribution of the two potential outcomes
#' sim_potential_outcomes = function(m,...){
#'   muibar = rnorm(1, 0, 1)
#'   Y0 = rnorm(m, muibar, 10) 
#'   Y1 = Y0 + 1
#'   re = cbind(Y0, Y1)
#'   colnames(re) = list("Y0", "Y1")
#'   return(re)
#' }
#'
#' ## Test statistics: pooled difference in the two study arms
#' calc_teststat = function(W, data, ...){
#'   Y0 = data[,"Y0"]
#'   Y1 = data[,"Y1"]
#'   re = sum(Y1[W==1])/sum(W==1) - sum(Y0[W==0])/sum(W==0)
#'   return(re)
#' }
#' ## search for the optimal number of clusters 
#' CRTsearch(nrep=1e4, nt=10, nc=10, FUN_clustersize=sim_cluster_size, FUN_Ys=sim_potential_outcomes, FUN_TestStat=calc_teststat)
#'}
#' @export
CRTsearch = function(nrep=1e4, nt, nc, tcRatio=1, minpower=0.8, alpha=0.05, increaseSamplingBy=1, PermutationTest=FALSE, Npermutationtest=100, ...){
  
  ## STEP 1: ################################################################################
  ## initiate a search list
  ## searching list is the list for seatching; created at the beginning of the program run and updated at the end or each searching stage    
  if( is.null(nt) && is.null(nc) ){
    stop("missing initial values for nc and nt")
  }else if( (!is.null(nt)) && (!is.null(nc)) ){
    ## check and set up searching ratio  for (nt, nc)
    if( nt/nc != tcRatio ){
      nt = round(nt*tcRatio, 0)
      print(paste("The search starts at nc=", nc, ", nt=", nt, "according to the specificed tcRatio=", tcRatio,sep=""))                
    }
  }else{
    if(is.null(nt)){
      nt = round(nc*tcRatio, 0)
    }else{
      nc = round(nt/tcRatio, 0)
    }
  }
  search_list = cbind(nc, nt, NA)
  colnames(search_list) = list("nc", "nt", "power")
  search_list = as.data.frame(search_list)
  
  
  ## STEP 2: ################################################################################
  ## update search_list according to current power and (nc, nt) using Binary Search
  
  ## calculate power with the initial sample size
  if(PermutationTest){
    power = simPowerPT(nrep=nrep, nt=nt, nc=nc, alpha=alpha,  ...)
  }else{
    power = simPower(nrep=nrep, nt=nt, nc=nc, alpha=alpha, ...)
  }
  print(paste("start searching for optimal design at: nc:", nc, "nt:", nt, "power:", round(power,5)))
  
  ## start the searching process
  searching = TRUE
  upperbound = NULL
  lowerbound = NULL
  while(searching){
    
    if(power <= minpower){
      lowerbound =  data.frame(nc,nt,power)
      if(is.null(upperbound)){
        nt=round(lowerbound$nt*2, 0);  nc=round(lowerbound$nc*2, 0)
      }else{
        nttotal = upperbound$nt + lowerbound$nt
        nctotal = upperbound$nc + lowerbound$nc
        nt=round(nttotal/2,0);  nc=round(nctotal/2,0)
      }
      if(power >= max(0, minpower-0.2) & power <= max(1, minpower*9.5/8)){nrep=min(round(increaseSamplingBy*nrep,0), 1e6)}
      if(PermutationTest){
        power = simPowerPT(nrep=nrep, nt=nt, nc=nc, alpha=alpha, Npermutationtest=Npermutationtest, ...)
      }else{
        power = simPower(nrep=nrep, nt=nt, nc=nc, alpha=alpha, ...)
      }
    }else if(power > minpower){
      upperbound  = data.frame(nc,nt,power)
      if(is.null(lowerbound)){
        nt=round(upperbound$nt/2, 0);  nc=round(upperbound$nc/2, 0)
        nt=ifelse(nt>3, nt, 3); nc=ifelse(nt>3, nc, 3)
      }else{
        nttotal = upperbound$nt + lowerbound$nt
        nctotal = upperbound$nc + lowerbound$nc
        nt=round(nttotal/2,0);  nc=round(nctotal/2,0)
      }
      if(power >= max(0, minpower-0.2) & power <= max(1, minpower*9.5/8)){nrep=min(round(increaseSamplingBy*nrep,0), 1e6)}
      if(PermutationTest){
        power = simPowerPT(nrep=nrep, nt=nt, nc=nc, alpha=alpha, Npermutationtest=Npermutationtest, ...)
      }else{
        power = simPower(nrep=nrep, nt=nt, nc=nc, alpha=alpha, ...)
      }
    }
    
    if(power <= minpower){
      lowerbound  = data.frame(nc,nt,power)        
    }else if(power > minpower){
      upperbound  = data.frame(nc,nt,power)        
    }
    print(paste("searching for optimal design at: nc:", nc, "nt:", nt, "power:", round(power,5)))
    
    ## updete search status
    if(nt<=3 & nc<=3 & power>=minpower){
      searching=FALSE
      upperbound  = data.frame(nc,nt,power)
      print(paste("nc", nc, "nt:", nt, "power:", power, "two small sample size for a Randomization test"))        
    }else if( is.null(upperbound) |  is.null(lowerbound) ){
      searching = TRUE
    }else if( (upperbound$nt+upperbound$nc)-(lowerbound$nt+lowerbound$nc) <= ceiling(2*tcRatio) ){
      searching = FALSE
    }else{
      searching = TRUE
    }                  
  }                
  
  print("optimal design found!")
  names(upperbound) = list("nc", "nt", "power")
  print(upperbound)
  return(upperbound)
}


#' estimate hypothesis test power for a CRT using simulation
#'
#' estimate hypothesis test power for a Cluster Randomization Trail (CRT) given the number of clusters in the treatment arm and the control arm for Neyman's Null hypothesis
#'
#' @param nrep the number of replications of the simulation procedure for generating the pseudo study data, default is \code{"1000"})
#' @param nt the number of cluster in the treatment arm to start the search
#' @param nc the number of cluster in the control arm to start the search
#' @param minpower the minimum power, default at 0.8
#' @param alpha type-1 errer, default at 0.05 
#' @param uppersided if a uppersided test is of interest, default at NULL (two-sided test), FALSE if a lowersided test if of interest
#' @param FUN_TestStat user-defined function for the test statistics, should take in two arguments, W=the treatment/control assignment indicator at the individual level; data=the dataset which includes the two potential outcomes at least
#' @inheritParams simulate_CRT
#' 
#' @return a vector (length=1) for the estimated power
#'
#' @examples
#'\dontrun{
#' ## distribution of cluster sizes
#' sim_cluster_size = function(N, ...){
#'   size = round(100*rnorm(N, 0, 1), 0)
#'   size[size<=0] = 1
#'   return(size)
#' }
#' ## distribution of the two potential outcomes
#' sim_potential_outcomes = function(m,...){
#'   muibar = rnorm(1, 0, 1)
#'   Y0 = rnorm(m, muibar, 10) 
#'   Y1 = Y0 + 1
#'   re = cbind(Y0, Y1)
#'   colnames(re) = list("Y0", "Y1")
#'   return(re)
#' }
#'
#' ## Test statistics: pooled difference in the two study arms
#' calc_teststat = function(W, data, ...){
#'   Y0 = data[,"Y0"]
#'   Y1 = data[,"Y1"]
#'   re = sum(Y1[W==1])/sum(W==1) - sum(Y0[W==0])/sum(W==0)
#'   return(re)
#' }
#' ## search for the optimal number of clusters 
#' simPower(nrep=1e4, nt=10, nc=10, FUN_clustersize=sim_cluster_size, FUN_Ys=sim_potential_outcomes, FUN_TestStat=calc_teststat)
#'}
#'
#' @export
simPower = function(nrep=1e4, nt, nc, alpha=0.05, FUN_TestStat, uppersided=NULL, ...){
  mcCores = parallel::detectCores()
  TH0THa = parallel::mclapply(1:nrep, TestStat_TH0THa, nt=nt, nc=nc, FUN_TestStat=FUN_TestStat, ..., mc.cores=mcCores-1)
  TH0THa = plyr::ldply(TH0THa, data.frame)
  power = simPowerTH0Ha(TH0=TH0THa$TH0, THa=TH0THa$THa, alpha=alpha, ...)
  return(power)
}  


#' estimate hypothesis test power for a CRT using simulation
#'
#' estimate hypothesis test power for a Cluster Randomization Trail (CRT) given the number of clusters in the treatment arm and the control arm for the Fisher's shrap null hypothesis and permutation test
#'
#' @param nrep the number of replications of the simulation procedure for generating the pseudo study data, default is \code{"1000"})
#' @param nt the number of cluster in the treatment arm to start the search
#' @param nc the number of cluster in the control arm to start the search
#' @param minpower the minimum power, default at 0.8
#' @param alpha type-1 errer, default at 0.05 
#' @param uppersided if a uppersided test is of interest, default at NULL (two-sided test), FALSE if a lowersided test if of interest
#' @param FUN_TestStat user-defined function for the test statistics, should take in two arguments, W=the treatment/control assignment indicator at the individual level; data=the dataset which includes the two potential outcomes at least
#' #' @param Npermutationtest the number of replications when estimating the permutation test power
#' @inheritParams simulate_CRT
#' 
#' @return a vector (length=1) for the estimated power
#'
#' @examples
#'\dontrun{
#' ## distribution of cluster sizes
#' sim_cluster_size = function(N, ...){
#'   size = round(100*rnorm(N, 0, 1), 0)
#'   size[size<=0] = 1
#'   return(size)
#' }
#' ## distribution of the two potential outcomes
#' sim_potential_outcomes = function(m,...){
#'   muibar = rnorm(1, 0, 1)
#'   Y0 = rnorm(m, muibar, 10) 
#'   Y1 = Y0 + 1
#'   re = cbind(Y0, Y1)
#'   colnames(re) = list("Y0", "Y1")
#'   return(re)
#' }
#'
#' ## Test statistics: pooled difference in the two study arms
#' calc_teststat = function(W, data, ...){
#'   Y0 = data[,"Y0"]
#'   Y1 = data[,"Y1"]
#'   re = sum(Y1[W==1])/sum(W==1) - sum(Y0[W==0])/sum(W==0)
#'   return(re)
#' }
#' ## search for the optimal number of clusters 
#' simPowerPT(nrep=1e4, nt=10, nc=10, FUN_clustersize=sim_cluster_size, FUN_Ys=sim_potential_outcomes, FUN_TestStat=calc_teststat)
#'}
#'
#' @export
simPowerPT = function (nrep = 10000, nt, nc, alpha = 0.05, FUN_TestStat, uppersided = NULL, Npermutationtest=100, ...){
  powersPT=NULL
  for(npt in 1:Npermutationtest){
    CRTsample = simulate_CRT(nt = nt, nc = nc, ...)
    mcCores = parallel::detectCores()
    assignment = parallel::mclapply(1:nrep, function(i, ...) {
      assignment_CRT(data = CRTsample, nt = nt, nc = nc, ...)
    }, mc.cores = mcCores - 1)
    assignment = as.data.frame(assignment)
    CRTsampleH0 = CRTsample
    CRTsampleH0$Y1 = CRTsampleH0$Y0
    TH0 = parallel::mclapply(1:nrep, function(i, ...) {
      FUN_TestStat(W = assignment[, i], data = CRTsampleH0, ...)
    }, mc.cores = mcCores - 1)
  THa = parallel::mclapply(1:nrep, function(i, ...) {
      FUN_TestStat(W = assignment[, i], data = CRTsample, ...) }
    , mc.cores = mcCores - 1)
    TH0 = unlist(TH0)
  THa = unlist(THa)
    power = simPowerTH0Ha(TH0 = TH0, THa = THa, alpha = alpha,...)
    powersPT = c(powersPT, power)
  }
  power = median(powersPT)
  return(power)
}

##powersPT = lapply(1:Npermutationtest, function(i, ...) {
##  CRTsample = simulate_CRT(nt = nt, nc = nc, ...)
##  mcCores = parallel::detectCores()
##  assignment = parallel::mclapply(1:nrep, function(i, 
##                                                   ...) {
##    assignment_CRT(data = CRTsample, nt = nt, nc = nc, 
##                   ...)
##  }, mc.cores = mcCores - 1)
##  assignment = as.data.frame(assignment)
##  CRTsampleH0 = CRTsample
##  CRTsampleH0$Y1 = CRTsampleH0$Y0
##  TH0 = parallel::mclapply(1:nrep, function(i, ...) {
##    FUN_TestStat(W = assignment[, i], data = CRTsampleH0, 
##                 ...)
##  }, mc.cores = mcCores - 1)
##  THa = parallel::mclapply(1:nrep, function(i, ...) {
##    FUN_TestStat(W = assignment[, i], data = CRTsample, 
##                 ...)
##  }, mc.cores = mcCores - 1)
##  TH0 = unlist(TH0)
##  THa = unlist(THa)
##  power = simPowerTH0Ha(TH0 = TH0, THa = THa, alpha = alpha, 
##                        ...)
##  return(power)
##}, ...)

##powersPT=NULL
##for(npt in 1:Npermutationtest){
##  CRTsample = simulate_CRT(nt = nt, nc = nc, ...)
##  mcCores = parallel::detectCores()
##  assignment = parallel::mclapply(1:nrep, function(i, ...) {
##    assignment_CRT(data = CRTsample, nt = nt, nc = nc, ...)
##  }, mc.cores = mcCores - 1)
##  assignment = as.data.frame(assignment)
##  CRTsampleH0 = CRTsample
##  CRTsampleH0$Y1 = CRTsampleH0$Y0
##  TH0 = parallel::mclapply(1:nrep, function(i, ...) {
##    FUN_TestStat(W = assignment[, i], data = CRTsampleH0, ...)
##  }, mc.cores = mcCores - 1)
##  THa = parallel::mclapply(1:nrep, function(i, ...) {
##    FUN_TestStat(W = assignment[, i], data = CRTsample, ...) }
##    , mc.cores = mcCores - 1)
##  TH0 = unlist(TH0)
##  THa = unlist(THa)
##  power = simPowerTH0Ha(TH0 = TH0, THa = THa, alpha = alpha,...)
##  powersPT = c(powersPT, power)
##}


##mcCores = parallel::detectCores()
##powersPT = parallel::mclapply(1:Npermutationtest, function(i, ...) {
##  CRTsample = simulate_CRT(nt = nt, nc = nc, ...)
##  mcCores = parallel::detectCores()
##  assignment = parallel::mclapply(1:nrep, function(i, 
##                                                   ...) {
##    assignment_CRT(data = CRTsample, nt = nt, nc = nc, 
##                   ...)
##  }, mc.cores = mcCores - 1)
##  assignment = as.data.frame(assignment)
##  CRTsampleH0 = CRTsample
##  CRTsampleH0$Y1 = CRTsampleH0$Y0
##  TH0 = parallel::mclapply(1:nrep, function(i, ...) {
##    FUN_TestStat(W = assignment[, i], data = CRTsampleH0, 
##                 ...)
##  }, mc.cores = mcCores - 1)
##  THa = parallel::mclapply(1:nrep, function(i, ...) {
##    FUN_TestStat(W = assignment[, i], data = CRTsample, 
##                 ...)
##  }, mc.cores = mcCores - 1)
##  TH0 = unlist(TH0)
##  THa = unlist(THa)
##  power = simPowerTH0Ha(TH0 = TH0, THa = THa, alpha = alpha, 
##                        ...)
##  return(power)
##}, ..., mc.cores = mcCores - 1)

#' simulate a smaple CRT study data
#' 
#' data generating mechanism
#'
#' @param nrep the number of replications of the simulation procedure for generating the pseudo study data, default is \code{"1000"})
#' @param nt the number of cluster in the treatment arm to start the search
#' @param nc the number of cluster in the control arm to start the search
#' @param FUN_clustersize user-defined function for the distribution of cluster sizes, should take in only one argument: n=the number of clsuters in the study
#' @param FUN_Ys user-defined function for the distribution of the two potential outcomes, should take in only one argument: m=the number of individual units in the a clsuter
#' @param dataset a dataset from which to generate the study data (no need to provide FUN_clustersize or FUN_Ys), should be individual level data, including columns for clutser ID, outcome of interest, and the variable for stratified randomization if desired.
#' @param outcome the name of the outcome variable in the dataset
#' @param clusterID the name of the cluster ID in the dataset
#' 
#' @return a dataframe
#'
#' @examples
#'\dontrun{
#' ## distribution of cluster sizes
#' sim_cluster_size = function(N, ...){
#'   size = round(100*rnorm(N, 0, 1), 0)
#'   size[size<=0] = 1
#'   return(size)
#' }
#' ## distribution of the two potential outcomes
#' sim_potential_outcomes = function(m,...){
#'   muibar = rnorm(1, 0, 1)
#'   Y0 = rnorm(m, muibar, 10) 
#'   Y1 = Y0 + 1
#'   re = cbind(Y0, Y1)
#'   colnames(re) = list("Y0", "Y1")
#'   return(re)
#' }
#' 
#' ## simulate a CRT study
#' simulate_CRT(nt=10, nc=10, FUN_clustersize=sim_cluster_size, FUN_Ys=sim_potential_outcomes, ...)
#'}
#'
#' @export
simulate_CRT = function(nt, nc, replacement=TRUE, FUN_clustersize=NULL, FUN_Ys=NULL, dataset=NULL, outcome=NULL, clusterID=NULL, ...){
  
  if(is.null(dataset)){
    if(is.null(FUN_clustersize) ){
      stop("please specify a function to simulate the size of clusters when no dataset is provided")
    }else{
      cluster_size = FUN_clustersize(nt+nc)
    }
    ## generate Yi(0) and Yi(1) 
    if(is.null(FUN_Ys)){
      stop("please specify a function to simulate the potential outcome under the control")
    }else{
      outcomes = lapply(cluster_size, simulate_Ys, FUN_Ys=FUN_Ys, ...)
      outcomes = plyr::ldply(outcomes, data.frame)
    }
    cID = rep(1:(nt+nc), cluster_size)
    re = as.data.frame(cbind(cID, outcomes))
  }else{
    if(is.null(FUN_Ys)){
      stop("please specify a function to simulate the potential outcomes under intervention")
    }
    colnames(dataset)[which(colnames(dataset)==clusterID)]="cID"
    colnames(dataset)[which(colnames(dataset)==outcome)]="Y0"
    cNdataset = length(unique(dataset$cID))
    if((!replacement)&(cNdataset<(nt+nc))){
      stop("required cluster size bigger than the one offered")
    }
    dataset$cID = as.numeric(dataset$cID)
    selected_cluster = sample(unique(dataset$cID), nt+nc, replace=replacement)
    real_baseline = function(i){
      baseline = dataset[dataset$cID==selected_cluster[i],]
      baseline$cID = i
      return(baseline)
    }
    baseline = lapply(1:(nt+nc), real_baseline)
    baseline = plyr::ldply(baseline, data.frame)
    ## Yi(1)
    outcomes = simulate_Ys(dataset=baseline, FUN_Ys=FUN_Ys, ...)
    baseline$Y0 = outcomes[,"Y0"]
    baseline$Y1 = outcomes[,"Y1"]
    re = baseline
  } 
  return(re)
}    


###################################################################
## other inner functions
###################################################################

## sample the two potential outcomes with/without a real dataset
simulate_Ys = function(ClusterSize=NULL, FUN_Ys=NULL, dataset=NULL, ...){
  if(is.null(FUN_Ys) & is.null(dataset)){
    stop("must specify FUN_Ys when no dataset is provided")
  }else if(is.null(dataset)){
    ## generate Ys when no dataset is provided
    Ys = FUN_Ys(ClusterSize)
  }else{
    ## generate Y1 when dataset is provided
    Ys = FUN_Ys(data=dataset, ...)
  }
  re = as.data.frame(Ys)
  return(re)
}

## approximate Power using the distribution of the test statistics under H0 and Ha
simPowerTH0Ha = function(TH0, THa, alpha=0.05, uppersided=NULL, ...){
  ## by default performs twosided test
  ## otherwise specify uppersided or not
  if(sum(is.na(TH0))>0){print(summary(TH0))}
  if(sum(is.na(TH0))>0){print(summary(THa))}
  if(is.null(uppersided)){
    Tstar = quantile(TH0, probs=c(alpha/2, 1-alpha/2))
    power = (sum(THa >= Tstar[2]) + sum(THa <= Tstar[1])) / length(THa)
  }else if(uppersided){
    Tstar = quantile(TH0, probs=1-alpha/2)
    power = sum(THa >= Tstar[1]) / length(THa)
  }else if(!uppersided){
    Tstar = quantile(TH0, probs=alpha/2)
    power = sum(THa <= Tstar[1]) / length(THa)
  }
  return(power)
}


## The value of the test statistics under H0 and Ha
TestStat_TH0THa = function(i, FUN_TestStat, nt, nc, ...){
  ## data generating mechanism
  CRTsample = simulate_CRT(nt=nt, nc=nc, ...)
  ## assignement mechanism
  assignment = assignment_CRT(data=CRTsample, nt=nt, nc=nc, ...)
  
  ## Test statistics under H0
  ##TH0 = calc_TestStat(assignment=assignment, Y0name="Y0", Y1name="Y0", data=CRTsample, ...)
  CRTsampleH0 = CRTsample
  CRTsampleH0$Y1 = CRTsampleH0$Y0
  TH0 = FUN_TestStat(assignment, CRTsampleH0, ...)
  ## Test statistics under Ha
  ##THa = calc_TestStat(assignment=assignment, Y0name="Y0", Y1name="Y1", data=CRTsample, ...)
  THa = FUN_TestStat(assignment, CRTsample, ...)

  re = cbind(TH0, THa)
  colnames(re) = list("TH0", "THa")
  return(re)
}


## randomization mechanism: assignment at the cluster level
assignment_CRT = function(nt, nc, data, stratifyBy=NULL, ...){
  ClusterAssignment = rep(0,nc+nt)
  if(is.null(stratifyBy)){
    TXTid=sample(1:(nc+nt), nt)
    ClusterAssignment[TXTid] = 1
  }else if(stratifyBy %in% colnames(data)){
    tmp = data[!duplicated(data$cID),c("cID", stratifyBy)]
    SizeStrata = as.vector(table(tmp[,stratifyBy]))
    Nstrata = length(SizeStrata)
    while(sum(ClusterAssignment)!=nt){
      ClusterAssignment = unlist(lapply(SizeStrata, function(n){ rbinom(n, 1, prob=1/2) }))
    }
    cIDorderByStrata = unique(data[,"cID"][order(data[,stratifyBy])])
    ClusterAssignment = ClusterAssignment[cIDorderByStrata]
  }else{
    stop(paste(stratifyBy, "is not found in the data"))
  }
  cluster_size = as.data.frame(table(data$cID))$Freq
  assignment = unlist(sapply(1:length(cluster_size), function(x){rep(ClusterAssignment[x],cluster_size[x])})) 
  return(assignment)
}

