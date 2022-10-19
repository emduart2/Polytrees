library(foreach)
library(utils)
library(pcalg)
library(igraph)
library(doParallel)

# - when using GIES, we need ensureDiff=TRUE 


#---- wrapper for estimating skeleton and orientations ----
# estimates the skeleton from interventional data.
# Inputs: 
#   Rs: list of sampled correlation matrices, 
#   Ns: numeric vector of samplesizes,
#   method: String. Computation method. One of "mean","median","ltestl". 
#       Default: "mean".
# Output: returns skeleton as igraph
estimate_skeleton <- function(Rs, Ns, method="mean"){
  matr = switch(method,
                mean = wmeanCorrels(Rs,Ns)$Rmean,
                median = wmedianCorrels(Rs,Ns)$Rmedian,
                ltest = Imatrix(Rs,Ns)
  )
  Cl = chowLiu(matr)
  return(Cl)
}


# function to estimate orientations 
# Inputs:
#   p:        Integer. Number of nodes
#   Covlist:  List of covariance matrices.
#   Ilist:    List of intervention targets. Each entry is numeric vector with 
#             first elemnt being sample size, all other elemnts nodes on which 
#             interventions were performed.
#   E:        2-by-e matrix where e is the number of edges in the skeleton. 
#             Each undirected edge should only be present once.
#   lC:       Imatrix
#   thres:    TODO
#   alpha:    significance level
#   Xlist:    List of samples with same length as Ilist. Each entry is n-by-p where
#             n is the number of samples
#   procedure:String to denote which procedure to take. One of 
#             "1","1simp","2","2simp","3","3simp". Default: "3simp"
#   pw_method:String. Pairwise method to use. One of "BIC","IRC". Default: "IRC"
#
# Outpus:
#   graphNEL: i_cpdag_graph with the estimated orientations
estimate_orientations <- function(
    p,Covlist,Ilist,Nlist,E,lC,thres,alpha,Xlist,procedure="3simp",pw_method="BIC"){
  
  # parse pw_method for internal naming
  if(pw_method == "IRC"){
    pw_method = "TEST"
  }
  
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
  
  # Random Orientation
  if(procedure=="random"){
    for(i in c(1:(p-1))){
      if(rbinom(1,1,0.5)==0){
        E[i,]<-rev(E[i,])
      }
    }
    dag_adj<-cpdag_from_lists(E,c(),p)
    i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
    return(i_cpdag_graph)
  }
  
  stop("Invalid procedure name.")
}




#---- functions to explore performance on a large scale ----
# function to simulate skeleton recovery
# Inputs:
#   df_params: data.frame. Containing following columns:
#                 tsize: denoting number of nodes
#                 totalSamples: number of samples
#                 interventionSize: size of each intervention
#                 ndatasets: number datasets. See documentation of isetting.
#                 k: iteration index
#                 sdatasets: see documentation of isetting.
#                 kindOfIntervention: which kind of intervention to simulate.
#                   See documentation of interventionalData
#                 ensureDiff: see documentation of isetting.
#             Each row is one setting to simulate skeleton recovery on.
#   scoreFct: Function to evaluate the estimations. Takes two adjacency matrices
#             and returns a numeric.
# Outputs:
#   df: data.frame where each row corresponds to a computation. The columns are
#         self-explanatory and a superset of the columns of df_params. 
#         Additional columns are the result of the scoreFunction and a column
#         indicating the aggregation function used, where "baseObs" denotes 
#         chowLiu on only the observational data and "baseAll" denotes chowLiu
#         on all data treating observational as interventional.
#   str String. Describes some properties of df_params
explore_skeleton <- function(df_params, scoreFct = SHD_skeleton){
  
  # main (parallel) computation loop
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
    # .verbose=TRUE,
    .combine = 'rbind'
  ) %dopar% {
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
    est <-estimate_skeleton(ID$Rs,ID$Ns,method="ltest")
    SHD_ltest_val<-scoreFct(as_adjacency_matrix(est), IS$gTrues)
    
    est <- estimate_skeleton(ID$Rs,ID$Ns, method="mean")
    SHD_mean_val<-scoreFct(as_adjacency_matrix(est), IS$gTrues)
    
    est <- estimate_skeleton(ID$Rs,ID$Ns, method="median")
    SHD_median_val<-scoreFct(as_adjacency_matrix(est), IS$gTrues)
    
    # baseline from chowLiu on only observational data
    maxNSample = as.integer(floor(totalSmpl/nds))
    # if nsample not evenly distributed on all datasets, ensure that baseline fair
    X = ID$Xs[[1]][1:maxNSample,]
    Cov = sample.cov(X)
    R = cov2cor(Cov)
    CL<-chowLiu(R)
    SHD_baseObs <-scoreFct(as_adjacency_matrix(CL), IS$gTrues)
    
    # baseline from chowLiu on all data treating interventional as observational data
    Xall = Reduce(rbind, ID$Xs,c())
    CovAll = sample.cov(Xall)
    Rall = cov2cor(CovAll)
    CL<-chowLiu(Rall)
    SHD_baseAll <-scoreFct(as_adjacency_matrix(CL), IS$gTrues)
    
    # return results as vector
    res <- c(SHD_ltest_val,SHD_mean_val,SHD_median_val,SHD_baseObs,SHD_baseAll)
    return(res)
  })
  
  # create nice result data frame
  df_median <- df_mean <- df_ltest <- df_baseObs <- df_baseAll <- df_params
  df_ltest$method <- "ltest"
  df_ltest$SHD <- vals[,1]
  df_mean$method <- "mean"
  df_mean$SHD <- vals[,2]
  df_median$method <- "median"
  df_median$SHD <- vals[,3]
  df_baseObs$method = "baseObs"
  df_baseObs$SHD = vals[,4]
  df_baseAll$method = "baseAll"
  df_baseAll$SHD = vals[,5]
  df <- rbind(df_ltest,df_mean,df_median,df_baseObs,df_baseAll)
  
  return(prepare_df_plot(df))
}


# function to simulate orientation recovery and the full algorithm.
# Inputs:
#   df_params: data.frame. Containing following columns:
#                tsize: denoting number of nodes
#                totalSamples: number of samples
#                interventionSize: size of each intervention
#                ndatasets: number datasets. See documentation of isetting.
#                k: iteration index
#                sdatasets: see documentation of isetting.
#                kindOfIntervention: which kind of intervention to simulate.
#                  See documentation of interventionalData
#                ensureDiff: see documentation of isetting.
#              Each row is one setting to simulate the algorithm on.
#   scoreFct_all: list of functions that take two graphNELs, where the first one
#             is the grountruth and second estimations and returns a numeric
#             value to evaluate the estimation.
#   sFctNames_all: string vector containing names of the scoreFcts to put into
#             the resulting data frame.
#   methods_all: list of numeric vectors with length 2 denoting the methods
#             to explore. First string is method fo skeleton recovery or 
#             "gtruth" for ground truth. Second string is for orientation 
#             mrecovery, see documentation of estimate_orientations. Can also be
#             c("GIES","GIES") to use Greedy interventional equivalence search
#             from Hauser and Buehlmann 2012.
#   pw_methods_all: string vector for denoting which pairwise recovery method
#             to use. String can be one of "BIC", "IRC".
# Outputs:
#   df: data.frame where each row corresponds to a computation. The columns are
#       self-explanatory and a superset of the columns of df_params. 
#       Additional columns are 
#         "shd_rand" denoting SHD of random polytree,
#         "time_s" denoting the runtime of the method, 
#         "method" denoting a  string that summarizes the choices for the recovery and
#         "[scoreFct]" for each scoreFct denoting the result of that function.
#   str String. Describes some properties of df_params
explore <- function(
    df_params, 
    scoreFct_all = list(pcalg::shd, SID), 
    sFctNames_all = c("SHD","SID"),
    methods_all = list(c("GIES","GIES"),c("mean","1"),c("mean","3"),c("gtruth","3")),
    pw_methods_all = c("BIC")){
  
  
  # main (parallel) computation loop
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
    # .verbose=TRUE,
    .errorhandling = "stop",
    .packages = c("mpoly")
  ) %dopar% {
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
    ID<-interventionalData(G,IS$L,IS$targetsI, kindOfIntervention=kindOfIntervention)
    
    
    # computations
    res = array(NaN, length(methods_all) * length(pw_methods_all) * (4+length(scoreFct_all)))
    i = 1
    for(m in methods_all){
      for(pw in pw_methods_all){
        
        
        # estimate essential graph
        Ilist<-IS$targetsI
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
          t = as.numeric(Sys.time() - s, units="secs")
          est = as(gies.fit$essgraph,"graphNEL")
          
        } else {
          
          # estimate skeleton
          s = Sys.time()
          Covlist<-ID$Covs
          Xlist<-ID$Xs
          Nlist<-ID$Ns      
          C_list<-ID$Rs
          for(j in c(1:length(Covlist))){
            Covlist[[j]]<-Covlist[[j]]*Nlist[j]
          }
          thres<-0.5*log(sum(Nlist))*length(Ilist)
          lC<-Imatrix(C_list,ID$Ns)
          if( m[1] == "gtruth") {
            E <- get.edgelist(G)
          } else {
            skelEst = estimate_skeleton(ID$Rs,ID$Ns, method=m[1])
            E <- get.edgelist(skelEst)
          }
          
          # estimate orientations
          est = tryCatch({
            estimate_orientations(p,Covlist,Ilist,Nlist,E,lC,thres,alpha,Xlist,
                                  procedure=m[2],pw_method=pw)},
            error = function(e){NULL})
          # est = estimate_orientations(p,Covlist,Ilist,Nlist,E,lC,thres,alpha,Xlist,
          #                         procedure=m[2],pw_method=pw)
          t = as.numeric(Sys.time() - s, units="secs") + ID$timeCovCor
        }
        
        # add additional information about the computation
        res[i] = t  # add runtime to results
        if(use_dags){ 
          nedges = dim(as_edgelist(G))[1]
          c = count_components(G,mode="weak")
          res[i+1] = nedges - (p-1) + (c-1)   # add optimal SHD
          res[i+2] = nedges                  # add number edges
        }
        g_rand = graph_from_adjacency_matrix(pruferwithskeletonOpt(p)$Directed)
        icpdag_rand = dag2essgraph(as_graphnel(g_rand),Ilist)
        true_i_cpdag<-dag2essgraph(as_graphnel(G),Ilist)
        res[i+3] = pcalg::shd(icpdag_rand, true_i_cpdag) # add shd of random polytree
        
        # score estimation and add to results
        j = 4
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
  dfs_counter = 1
  i = 1
  for(m in methods_all){
    for(pw in pw_methods_all){
      
      # add methods used
      df <- df_params
      df$aggrFct <- m[1]
      df$proc <- m[2]
      df$pw_method <- pw
      df$method <- paste(m[1],m[2],pw,sep=',') # for easier grouping in plot
      
      # add additional metrics
      df$time_s <- vals[,i]
      if(unique(df_params$use_dags)){
        df$SHD_optimal = vals[,i+1]
        df$nEdgesDAG = vals[,i+2]
      }
      df$shd_rand = vals[,i+3]
      
      # add scores
      j = 1
      for(sf in scoreFct_all){
        df[,sFctNames_all[[j]]] <- vals[,i+3+j]
        j = j+1
      }
      i = i + 3 + j
      
      # add to list
      dfs[[dfs_counter]] = df
      dfs_counter = dfs_counter + 1
    }
  }
  dfAll <- Reduce(rbind, dfs, data.frame())
  
  # save results in one variable
  l = prepare_df_plot(dfAll)
  return(l)
}


# internal function to prepare the result data frame for plotting
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
  if((unique(df$use_dags))==TRUE){
    str = append(str,paste(", dags with nbh ", unique(df$dag_nbh)))
  }
  str <- paste(str,collapse="")
  str <- substr(str,3,1000000L)
  return(list(df=df,str=str))
}
