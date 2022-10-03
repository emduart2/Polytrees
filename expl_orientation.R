library(foreach)
library(utils)
library(pcalg)
library(igraph)

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
#   ensureDiff = TRUE,
#   alpha = 0.05
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
estimate_orientations <- function(
    p,Covlist,Ilist,Nlist,E,lC,thres,alpha,Xlist,procedure="1",pw_method="BIC"){
  
  # procedure 1 simple
  if(procedure == "1simp"){
    dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,E,lC,thres,method="simple"
                               ,pw_method=pw_method,alpha,Xlist)
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
    return(i_cpdag_graph)
  }
  
  # procedure 1 (with trycatch for failure of optimization)
  if(procedure == "1"){
    dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,E,lC,thres,method="nothing"
                               ,pw_method=pw_method,alpha,Xlist)
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
    return(i_cpdag_graph) 
  }
  
  # procedure 2 simple
  if(procedure == "2simp"){
    dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,E,p,method="simple",
                                   pw_method=pw_method,alpha,Xlist)
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
    return(i_cpdag_graph)
  } 
  
  # procedure 2 (with trycatch for failure of optimization)
  if(procedure == "2"){
    dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,E,p,method="nothing",
                                   pw_method=pw_method,alpha,Xlist)
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
    return(i_cpdag_graph)
  }
  
  # procedure 3 simple
  if(procedure == "3simp"){
    dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,E,p,method="simple"
                             ,pw_method=pw_method,alpha,Xlist)
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
    
    return(i_cpdag_graph)
  } 
  
  # procedure 3 (with catch for failure of optimization)
  if(procedure == "3"){
    dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,E,p,method="nothing"
                             ,pw_method=pw_method,alpha,Xlist)
    dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
    i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
    return(i_cpdag_graph)
  }
  
  stop("Invalid procedure name.")
}


# notes:
# - to execute in parallel, uncomment corresponding code blocks and change
#   "%do%" into "%dopar%"
#   WARNING: parallel execution might run into problem on windows machines!
orientationExploration <- function(
    df_params, 
    scoreFct_all = list(pcalg::shd, SID), 
    sFctNames_all=c("SHD","SID"),
    procedures_all = c("1","1simp","2","2simp","3","3simp"),
    pw_methods_all = c("BIC","TEST")){
  
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
    alpha = df_params$alpha,
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
    E <- get.edgelist(G)
    
    # computations
    res = array(NaN, length(procedures_all) * length(pw_methods_all) * length(scoreFct_all))
    i = 1
    for(proc in procedures_all){
      for(pw in pw_methods_all){
        
        # estimate orientation
        est = tryCatch({
          estimate_orientations(p,Covlist,Ilist,Nlist,E,lC,thres,alpha,Xlist,
                                       procedure=proc,pw_method=pw)},
          error = function(e){NULL})
        
        # score estimation
        j = 0
        for(sf in scoreFct_all){
          if( is.null(est)){
            out = NaN
          } else {
            out = sf(true_i_cpdag, est)
          }
          res[i+j] = out
          j = j + 1
        }
        i = i + j
      }
    }
    return(res)
  })
  
  # create nice result data frame
  dfs = vector("list", length(procedures_all) * length(pw_methods_all) * length(scoreFct_all))
  i = 1
  for(proc in procedures_all){
    for(pw in pw_methods_all){
      dfs[[i]] <- df_params
      dfs[[i]]$method <- proc
      dfs[[i]]$pw_method <- pw
      j = 0
      for(sf in scoreFct_all){
        dfs[[i]][,sFctNames_all[[j+1]]] <- vals[,i+j]
        j = j+1
      }
      i = i + j
    }
  }
  df <- Reduce(rbind, dfs, data.frame())
  
  # save results in one variable
  l = prepare_df_plot(df)
  return(l)
}
