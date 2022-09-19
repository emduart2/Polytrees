# notes:
# - to execute in parallel, uncomment corresponding code blocks and change
#   "%do%" into "%dopar%"
#   WARNING: parallel execution might run into problem on windows machines!

library(foreach)
library(utils)
library(progressr)

# parameter space to explore (sdatasets -1 is c(),
# interventionSize >= 1 is the absolute value of nodes to intervene)
df_params <- expand.grid(
  tsize = c(100),
  totalSamples = c(200),
  interventionSize = c(0.2,0.8),
  ndatasets = c(2,12),
  k = c(1:20),
  sdatasets = list(c())
)
comment <- ""

# setup parallel cluster
# n.cores <- parallel::detectCores() - 1
# my.cluster <- parallel::makeCluster(
#   n.cores, 
#   type = "FORK"
# )
# doParallel::registerDoParallel(cl = my.cluster)

# main (parallel) computation loop
# (if parallel not wished, just change %dopar% into %do%)
(vals <- foreach(
  p = df_params$tsize,
  totalSmpl = df_params$totalSamples,
  intervSize = df_params$interventionSize,
  nds = df_params$ndatasets,
  k = df_params$k,
  sd = df_params$sdatasets,
  .combine = 'rbind',
  .verbose=TRUE
) %do% {
  # parse parameters
  if(length(sd) == 1 && sd == -1) {
    sd <- c()
  } 
  if(intervSize <1){
    iss <- ceiling(p*intervSize)
  } else {
    iss <- intervSize
  }
  
  # create setting and sample interventional correlation matrices
  IS <- isetting(p,nds,iss,sd,totalSmpl)
  G <- graph_from_adjacency_matrix(IS$gTrued)
  ID <- interventionalData(G,IS$L,IS$targetsI)
  C_list <- ID$Rs
  
  
  # compute deviation from true skeleton for 3 diff aggregation functions
  lC<-Imatrix(C_list,ID$Ns)
  CL<-chowLiu(lC)
  E_e<-get.edgelist(CL)
  estimated_skeleton<-get.adjacency(CL)
  SHD_ltest_val<-sum(abs(estimated_skeleton-IS$gTrues))/(2*(p-1))
  
  meanC<-wmeanCorrels(C_list,ID$Ns)$Rmean
  CL<-chowLiu(meanC)
  E_e<-get.edgelist(CL)
  estimated_skeleton<-get.adjacency(CL)
  SHD_mean_val<-sum(abs(estimated_skeleton-IS$gTrues))/(2*(p-1))
  
  medianC<-wmedianCorrels(C_list,ID$Ns)$Rmedian
  CL<-chowLiu(medianC)
  E_e<-get.edgelist(CL)
  estimated_skeleton<-get.adjacency(CL)
  SHD_median_val<-sum(abs(estimated_skeleton-IS$gTrues))/(2*(p-1))
  
  # return results as vector
  res <- c(SHD_ltest_val,SHD_mean_val,SHD_median_val)
  return(res)
})

# create nice result data frame
df_median <- df_mean <- df_ltest <- df_params
df_ltest$method <- "ltest"
df_ltest$SHD <- vals[,1]
df_mean$method <- "mean"
df_mean$SHD <- vals[,2]
df_median$method <- "median"
df_median$SHD <- vals[,3]
df <- rbind(df_ltest,df_mean,df_median)

# save results in one variable
compRes$params <- df_params
compRes$df_vals <- df
compRes$comment <- comment
if(! exists("allResultsNew")) {allResultsNew = list()}
allResultsNew <- append(allResultsNew,list(compRes))

# close parallel cluster
# parallel::stopCluster(cl = my.cluster)
