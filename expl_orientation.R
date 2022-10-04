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
#   alpha = c(0.05)
# )
# methods_all = list(c("GIES"),c("mean","1"),c("true","2"))

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
  
  # procedure 1
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
  
  # procedure 2
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
  
  # procedure 3
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
# - to execute in parallel, initialize cluster and change "%do%" into "%dopar%"
#   WARNING: parallel execution might run into problem on windows machines!
orientationExploration <- function(
    df_params, 
    scoreFct_all = list(pcalg::shd, SID), 
    sFctNames_all = c("SHD","SID"),
    methods_all = list(c("GIES","GIES"),c("mean","1"),c("mean","3"),c("gtruth","3")),
    pw_methods_all = c("BIC")){
  
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
    use_dags = df_params$use_dags,
    dag_nbh = df_params$dag_nbh,
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
    IS <- isetting(p,nds,iss,sd,totalSmpl, ensureDiff = ensureDiff, 
                   use_dags=use_dags, dag_nbh=dag_nbh)
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
    
    # computations
    res = array(NaN, length(methods_all) * length(pw_methods_all) * (1+length(scoreFct_all)))
    i = 1
    for(m in methods_all){
      for(pw in pw_methods_all){
        
        # estimate essential graph
        if(m[1] == "GIES") {
          
          # parse into GIES setting
          allSamples = matrix(0, totalSmpl, p)
          targetindex = array(0, totalSmpl)
          targets = vector("list", length(IS$targetsI))
          a = b = 1
          for(j in 1:length(IS$targetsI)){
            b <- a + IS$targetsI[[j]][1]-1
            allSamples[a:b,] = ID$Xs[[j]]
            targetindex[a:b] = j
            targets[[j]] = IS$targetsI[[j]][-1]
            a <- b+1
          }
          setting_GIES = new("GaussL0penIntScore", 
                             data=allSamples, 
                             targets=targets, 
                             target.index=targetindex)
          
          # GIES estimate
          s = Sys.time()
          gies.fit = gies(setting_GIES)
          t = as.numeric(Sys.time() - s) * 1000
          est = as(gies.fit$essgraph,"graphNEL")
          
        } else {
          
          # estimate skeleton
          s = Sys.time()
          if( m[1] == "gtruth") {
            E <- get.edgelist(G)
          } else {
            skelEst = estimate_skeleton(ID$Rs,ID$Ns, method=m[1])
            E <- get.edgelist(graph_from_adjacency_matrix(skelEst))
          }
          
          # estimate orientations
          est = tryCatch({
            estimate_orientations(p,Covlist,Ilist,Nlist,E,lC,thres,alpha,Xlist,
                                  procedure=m[2],pw_method=pw)},
            error = function(e){NULL})
          t = as.numeric(Sys.time() - s) * 1000
        }


        # score estimation and add to results
        res[i] = t # add runtime to results
        j = 1
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
  dfs = vector("list", length(methods_all) * length(pw_methods_all))
  i = 1
  for(m in methods_all){
    for(pw in pw_methods_all){
      dfs[[i]] <- df_params
      dfs[[i]]$aggrFct <- m[1]
      dfs[[i]]$proc <- m[2]
      dfs[[i]]$pw_method <- pw
      
      dfs[[i]]$time <- vals[,i]
      j = 1
      for(sf in scoreFct_all){
        dfs[[i]][,sFctNames_all[[j]]] <- vals[,i+j]
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
