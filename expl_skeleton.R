library(foreach)
library(utils)
library(progressr)

# computes the SHD of two skeletons given as (symmetric) adjacency
# matrices. If number of nodes larger three and the skeletons are both
# trees, the result lies in [0,1]
SHD_skeleton <- function(S1,S2){
  return(sum(abs(S1-S2))/2)
}

# estimates the skeleton from interventional data.
# Inputs: Rs=list of sampled correlation matrices, Ns=samplesizes,
#   method=computation method.
estimate_skeleton <- function(Rs, Ns, method="mean"){
  matr = switch(method,
               mean = wmeanCorrels(Rs,Ns)$Rmean,
               median = wmedianCorrels(Rs,Ns)$Rmedian,
               ltest = Imatrix(Rs,Ns)
         )
  Cl = chowLiu(matr)
  return(get.adjacency(Cl))
}



# example parameter space to explore
# (interventionSize >= 1 is the absolute value of nodes to intervene)
df_params <- expand.grid(
  tsize = c(100),
  totalSamples = c(300),
  interventionSize = c(1,10),
  ndatasets = c(10),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("random","imperfect","inhibitory","perfect"),
  conservative = FALSE
)


# notes:
# - to execute in parallel, uncomment corresponding code blocks and change
#   "%do%" into "%dopar%"
#   WARNING: parallel execution might run into problem on windows machines!
skeletonExploration <- function(df_params,allResults=list()){
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
    kindOfIntervention = df_params$kindOfIntervention,
    ensureDiff = df_params$ensureDiff,
    .packages = c("mpoly"),
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
    IS <- isetting(p,nds,iss,sd,totalSmpl, ensureDiff=ensureDiff)
    G <- graph_from_adjacency_matrix(IS$gTrued)
    ID <- interventionalData(G,IS$L,IS$targetsI,kindOfIntervention=kindOfIntervention)
    
    # compute deviation from true skeleton for 3 diff aggregation functions
    estimated_skeleton<-estimate_skeleton(ID$Rs,ID$Ns,method="ltest")
    SHD_ltest_val<-SHD_skeleton(estimated_skeleton, IS$gTrues)

    estimated_skeleton<-estimate_skeleton(ID$Rs,ID$Ns, method="mean")
    SHD_mean_val<-SHD_skeleton(estimated_skeleton, IS$gTrues)
    
    estimated_skeleton<-estimate_skeleton(ID$Rs,ID$Ns, method="median")
    SHD_median_val<-SHD_skeleton(estimated_skeleton, IS$gTrues)
    
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
  if(! exists("allResults")) {allResults = list()}
  allResults <- append(allResults,list(compRes))
  
  # close parallel cluster
  # parallel::stopCluster(cl = my.cluster)
  
  return(prepare_df_plot(df))
}


# prepare data for plotting
prepare_df_plot <- function(df){
  # convert to factors/strings
  df$totalSamples <- factor(df$totalSamples, levels=sort(unique(df$totalSamples),decreasing = FALSE))
  df$interventionSize <- factor(df$interventionSize, levels=sort(unique(df$interventionSize),decreasing = FALSE))
  df$ndatasets <- factor(df$ndatasets, levels=sort(unique(df$ndatasets),decreasing = FALSE))
  df$tsize <- factor(df$tsize, levels=sort(unique(df$tsize),decreasing = FALSE))
  df$sdatasets <- as.character(lapply(df$sdatasets, FUN=function(x){paste(x,collapse=", ")}))
  
  # create parameter strings
  str = ""
  if(length(unique(df$tsize))==1){
    str = append(str,paste(", nodes: ", unique(df$tsize)))
  }
  if(length(unique(df$totalSamples))==1){
    str = append(str,paste(", samples: ", unique(df$totalSamples)))
  }
  if(length(unique(df$interventionSize))==1){
    str = append(str,paste(", interventionSize: ", unique(df$interventionSize)))
  }
  if(length(unique(df$ndatasets))==1){
    str = append(str,paste(", ndatasets: ", unique(df$ndatasets)))
  }
  if(length(unique(df$kindOfIntervention))==1){
    str = append(str,paste(", intervention: ", unique(df$kindOfIntervention)))
  }
  if(length(unique(df$conservative))==1){
    str = append(str,paste(", conservative: ", unique(df$conservative)))
  }
  str <- paste(str,collapse="")
  str <- substr(str,3,1000000L)
  return(list(df=df,str=str))
}

