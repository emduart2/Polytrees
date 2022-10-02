library(foreach)
library(utils)
library(progressr)

# example parameter space to explore
# (interventionSize >= 1 is the absolute value of nodes to intervene)
# df_params <- expand.grid(
#   tsize = c(100),
#   totalSamples = c(300),
#   interventionSize = c(1,10),
#   ndatasets = c(10),
#   k = c(1:20),
#   sdatasets = list(c()),
#   kindOfIntervention = c("random","imperfect","inhibitory","perfect"),
#   conservative = FALSE,
#   estimateSkeleton = c("mean","FALSE")
# )

# setup parallel cluster
# n.cores <- parallel::detectCores() - 1
# my.cluster <- parallel::makeCluster(
#   n.cores,
#   type = "FORK"
# )
# doParallel::registerDoParallel(cl = my.cluster)

# close parallel cluster
# parallel::stopCluster(cl = my.cluster)


# function to estimate orientations (not used at the moment)
estimate_orientations <- function(p,Covlist,Ilist,Nlist,ESkel,lC,thres,procedure="1"){
  
  # procedure 1 simple
  if(procedure == "1simp"){
    dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,ESkel,lC,thres,method="simple")
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    i_cpdag_adj<-i_cpdag(Ilist,dag_adj)
    i_cpdag_graph<-as(i_cpdag_adj,"graphNEL")
    return(i_cpdag_graph)
  }
  
  # procedure 1 (with trycatch for failure of optimization)
  if(procedure == "1"){
    tryCatch({
      dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,ESkel,lC,thres,method="nothing");
      dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p);
      i_cpdag_adj<-i_cpdag(Ilist,dag_adj);
      i_cpdag_graph<-as(i_cpdag_adj,"graphNEL");
    },
    error = function(e){i_cpdag_graph <<- NULL})
    return(i_cpdag_graph) 
  }
  
  # procedure 2 simple
  if(procedure == "2simp"){
    dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,ESkel,p,method="simple")
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    i_cpdag_graph<-as(dag_adj,"graphNEL")
    return(i_cpdag_graph)
  } 
  
  # procedure 2 (with trycatch for failure of optimization)
  if(procedure == "2"){
    tryCatch({
      dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,ESkel,p,method="nothing");
      dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p);
      i_cpdag_graph<-as(dag_adj,"graphNEL");
    },
    error = function(e){i_cpdag_graph <<- NULL})
    return(i_cpdag_graph)
  }
  
  # procedure 3 simple
  if(procedure == "3simp"){
    dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,ESkel,p,method="simple")
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    i_cpdag_adj<-i_cpdag(Ilist,dag_adj)
    i_cpdag_graph<-as(i_cpdag_adj,"graphNEL")
    return(i_cpdag_graph)
  } 
  
  # procedure 3 (with catch for failure of optimization)
  if(procedure == "3"){
    tryCatch({
      dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,ESkel,p,method="nothing");
      dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p);
      i_cpdag_adj<-i_cpdag(Ilist,dag_adj);
      i_cpdag_graph<-as(i_cpdag_adj,"graphNEL");
    },
    error = function(e){i_cpdag_graph <<- NULL})
    return(i_cpdag_graph)
  }
  
  stop("Invalid procedure name.")
}


# notes:
# - to execute in parallel, uncomment corresponding code blocks and change
#   "%do%" into "%dopar%"
#   WARNING: parallel execution might run into problem on windows machines!
orientationExploration <- function(df_params,allResults=list(), allOneNodeIntv=FALSE){
  
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
    estimateSkeleton = df_params$estimateSkeleton,
    ensureDiff = df_params$ensureDiff,
    .combine = 'rbind',
    .verbose=TRUE,
    .errorhandling = "stop",
    .packages = c("mpoly")
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
    
    # create setting
    IS <- isetting(p,nds,iss,sd,totalSmpl, ensureDiff = ensureDiff)
    G<-graph_from_adjacency_matrix(IS$gTrued)
    ID<-interventionalData(G,IS$L,IS$targetsI)
    Covlist<-ID$Covs
    Xlist<-ID$Xs
    Ilist<-IS$targetsI
    Nlist<-ID$Ns      
    C_list<-ID$Rs
    for(i in c(1:length(Covlist))){
       Covlist[[i]]<-Covlist[[i]]*Nlist[i]
    }
    true_i_cpdag<-dag2essgraph(as_graphnel(G),Ilist)
    thres<-0.5*log(sum(Nlist))*length(Ilist)
    lC<-Imatrix(C_list,ID$Ns)
    if(estimateSkeleton != "FALSE"){
      skelEst = estimate_skeleton(ID$Rs,ID$Ns, method=estimateSkeleton)
      ESkel <- get.edgelist(graph_from_adjacency_matrix(skelEst))
    } else {
      ESkel <- get.edgelist(G)
    }
    
    
    ##Orientation with original skeleton strategy 1
    dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,ESkel,lC,thres,method="simple")
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    e_G<-graph_from_adjacency_matrix(dag_adj)
    i_cpdag_graph<-dag2essgraph(as_graphnel(e_G),Ilist)
    shd_1_simp <-shd(true_i_cpdag,i_cpdag_graph)
    
    tryCatch({
      dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,ESkel,lC,thres,method="nothing");
      dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p);
      e_G<-graph_from_adjacency_matrix(dag_adj)
      i_cpdag_graph<-dag2essgraph(as_graphnel(e_G),Ilist)
      shd_1 <-shd(true_i_cpdag,i_cpdag_graph)
    },
    error = function(e){shd_1 <<- NaN}
    )
    
    ##Orientation with original skeleton strategy 2
    dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,ESkel,p,method="simple")
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    i_cpdag_graph<-as(dag_adj,"graphNEL")
    shd_2_simp <- shd(true_i_cpdag,i_cpdag_graph)
    
    tryCatch({
      dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,ESkel,p,method="nothing");
      dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p);
      i_cpdag_graph<-as(dag_adj,"graphNEL");
      shd_2 <-shd(true_i_cpdag,i_cpdag_graph)
    },
    error = function(e){shd_2 <<- NaN}
    )
    
    
    ##Orientation with original skeleton strategy 3
    dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,ESkel,p,method="simple")
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    i_cpdag_adj<-i_cpdag(Ilist,dag_adj)
    i_cpdag_graph<-as(i_cpdag_adj,"graphNEL")
    shd_3_simp <-shd(true_i_cpdag,i_cpdag_graph)
    
    tryCatch({
      dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,ESkel,p,method="nothing");
      dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p);
      i_cpdag_adj<-i_cpdag(Ilist,dag_adj);
      i_cpdag_graph<-as(i_cpdag_adj,"graphNEL");
      shd_3 <-shd(true_i_cpdag,i_cpdag_graph)
    },
    error = function(e){shd_3 <<- NaN}
    )
    
    
    # return results as vector
    res <- c(shd_1_simp,shd_1, shd_2_simp,shd_2, shd_3_simp,shd_3, ID$edgesIntervened)
    return(res)
  })
  
  # create nice result data frame
  df_1simp <- df_1 <- df_2simp <- df_2 <- df_3simp <- df_3 <- df_params
  df_1simp$method <- "1simp"
  df_1simp$SHD <- vals[,1]
  df_1simp$edgIntv <- vals[,7]
  df_1$method <- "1"
  df_1$SHD <- vals[,2]
  df_1$edgIntv <- vals[,7]
  df_2simp$method <- "2simp"
  df_2simp$SHD <- vals[,3]
  df_2simp$edgIntv <- vals[,7]
  df_2$method <- "2"
  df_2$SHD <- vals[,4]
  df_2$edgIntv <- vals[,7]
  df_3simp$method <- "3simp"
  df_3simp$SHD <- vals[,5]
  df_3simp$edgIntv <- vals[,7]
  df_3$method <- "3"
  df_3$SHD <- vals[,6]
  df_3$edgIntv <- vals[,7]
  df <- rbind(df_1simp,df_1,df_2simp,df_2,df_3simp,df_3)
  
  # save results in one variable
  compRes$params <- df_params
  compRes$df_vals <- df
  if(! exists("allResults")) {allResults = list()}
  allResults <- append(allResults,list(compRes))
  l = prepare_df_plot(df)
  return(l)
}
