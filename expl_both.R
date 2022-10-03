library(foreach)
library(utils)

# script to explore the whole algorithm (i.e. skeleton and orientation estimation)


# setup parallel cluster
# n.cores <- parallel::detectCores() - 1
# my.cluster <- parallel::makeCluster(
#   n.cores,
#   type = "FORK"
# )
# doParallel::registerDoParallel(cl = my.cluster)

# close parallel cluster
# parallel::stopCluster(cl = my.cluster)

# methods = list(c("mean","p1), c("GIES"), c("mean","p3"))

# notes:
# - to execute in parallel, uncomment corresponding code blocks and change
#   "%do%" into "%dopar%"
#   WARNING: parallel execution might run into problem on windows machines!
bothExploration <- function(df_params,methods){
  
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
    for(i in c(1:length(Covlist))){
          Covlist[[i]]<-Covlist[[i]]*Nlist[i]
        }
    Xlist<-ID$Xs
    Ilist<-IS$targetsI
    Nlist<-ID$Ns      
    C_list<-ID$Rs
    true_i_cpdag<-dag2essgraph(as_graphnel(G),Ilist)
    thres<-0.5*log(sum(Nlist))*length(Ilist)
    lC<-Imatrix(C_list,ID$Ns)
    

    # calculations
    res = matrix(0,length(methods),5)
    for(i in 1:length(methods)){
      
      if(methods[[i]][1] == "GIES"){
        
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
        time_ort = as.numeric(Sys.time() - s) * 1000
        
        # assign outputs
        gies_est = as(gies.fit$essgraph,"graphNEL")
        shd_all = shd(as(as_adjacency_matrix(G),"graphNEL"), gies_est)
        shd_skel = NaN
        time_skel = 0
        
      } else {
        
        # estimate skeleton
        s = Sys.time()
        skelEst = estimate_skeleton(ID$Rs,ID$Ns, method=methods[[i]][1])
        time_skel = as.numeric(Sys.time() - s) * 1000
        ESkel = get.edgelist(graph_from_adjacency_matrix(skelEst))
        shd_skel = SHD_skeleton(skelEst,IS$gTrues)
        
        # estimate orientations
        s = Sys.time()
        est = estimate_orientations(p,Covlist,Ilist,Nlist,ESkel,lC,thres,procedure=methods[[i]][2])
        time_ort = as.numeric(Sys.time() - s) * 1000
        if(is.null(est)){
          sdh_all = NaN
        } else {
          shd_all = shd(est, true_i_cpdag)
        }
      }
      res[i,] = c(shd_all, shd_skel, time_ort, time_skel, ID$edgesIntervened)
    }

    return(res)
  })
  
  # create nice result data frame
  df = data.frame()
  for(i in 1:length(methods)){
    tmp = vals[i + (0:(dim(df_params)[1]-1))*length(methods),]
    df_tmp = df_params
    df_tmp$method = paste(methods[[i]],collapse = ", ")
    df_tmp$shd_all = tmp[,1]
    df_tmp$shd_skel = tmp[,2]
    df_tmp$time_ort = tmp[,3]
    df_tmp$time_skel = tmp[,4]
    df_tmp$time_all = df_tmp$time_ort + df_tmp$time_skel
    df_tmp$edgIntv = tmp[,5]
    df = rbind(df,df_tmp)
  }
  
  # save results in one variable
  l = prepare_df_plot(df)
  return(l)
}
