#-----------------------------
# Needed setup

# source("http://bioconductor.org/biocLite.R") 
# biocLite("RBGL")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("RBGL")
# BiocManager::install("graph")
# BiocManager::install("Rgraphviz")
# 
# install.packages("bnlearn")
# install.packages("igraph")
# install.packages("matrixcalc")
# install.packages("BNSL")
# install.packages("pcalg")
# install.packages("graph")
# install.packages("abind")
# install.packages("SID")
# install.packages("rmarkdown")
# install.packages('spatstat')
# install.packages("poolr")
# install.packages("harmonicmeanp")
# install.packages("DescTools")
# install.packages("igraph")
# install.packages("mpoly")
library(BNSL)
library(bnlearn)
library(igraph)
library(matrixcalc)
library(graphics)
library(pcalg)
library(abind)
# library(SID)
library(tidyverse)
# library(ggplot2)
# library(rmarkdown)
library(spatstat)
library(poolr)
library(harmonicmeanp)
library(DescTools)
library(Matrix)
library(graph)
library(mpoly)
library(sets)
library(gap)
#-------------------------

# Function for calculating sample covariance matrix
# https://github.com/sqyu/CorDiffViz/tree/master/R
sample.cov <- function(X, Y=NULL){ # Calculate sample covariance matrix
  n <- nrow(X)
  Xr <- X - t(matrix(colMeans(X), ncol(X), n))
  if (is.null(Y)){
    return(t(Xr)%*%Xr/n)
  } else {
    if (nrow(Y) != n) stop("X and Y must have the same number of rows.")
    Yr <- Y - t(matrix(colMeans(Y), ncol(Y), n))
    return(t(Xr)%*%Yr/n)
  }
}

#-----------------------------------------------------------
# Sampling n values independently 
# uniformly at random from two disjoint intervals (-maxx,-minn)U(minn,maxx)
ruunif<-function(n,minn,maxx){
  l<-c()
  for (i in c(1:n)){
    c<-rbinom(1,1,0.5)
    if(c==0){ 
      l<-append(l,runif(1,-maxx,-minn))
    }
    else{
      l<-append(l,runif(1,minn,maxx))
    }
  }
  return(l)
}

#----------
# Creating a random matrix of coefficients Lambda
# starting from a dag G and its nodes
# G = Dag,
coeffLambda<-function(G){
  p<-gorder(G)
  Lambda<-matrix(0,p,p)
  edg<-as_edgelist(G) # This is the list of directed edges
  lambdaCoeff<-ruunif(nrow(edg),0.3,1) #this is a random vector of coefficients
  # Upper triangular matrix with Lambda coefficients
  if(nrow(edg) > 0){
    for(j in (1:nrow(edg))){
      r<-edg[j,1]
      c<-edg[j,2]
      Lambda[r,c]<- lambdaCoeff[j]
    }
  }
  return(Lambda)
}
#-----------------
# Creating a sample from a dagModel specified by a
# matrix of coefficients Lambda
# Inputs:
#   n: Integer. Number of samples.
#   Lambda: Numeric Matrix. Matrix of coefficients
#   varvector: vector of error variances
# Outputs:
#   Matrix: n-by-dim(Lambda,1) matrix with samples.
samplingDAG<-function(n,Lambda,varvector){
  p<-nrow(Lambda)
  eps<-matrix(0,p,n)
  for(i in c(1:p)){
    eps[i,]<-rnorm(n,mean=0,varvector[i])
  }
  #eps <- matrix( rnorm(p*n,mean=0,sd=1), p, n) 
  lambdaInv<-solve(diag(rep(1,p))-Lambda)
  # XmatIncorrect <- t(lambdaInv %*% eps)
  XmatCorrect <- t(t(lambdaInv) %*% eps)
  
  return(XmatCorrect)
}
#------------------------
#---------------------------------------------------------------------
# This functions learns a skeleton from the absolute value
# of a weight matrix. We will usually input a correlation matrix
# to this function.
chowLiu<-function(R){
  edge.list<-kruskal(abs(R))
  g<- graph_from_edgelist(edge.list,directed = FALSE)
  return(g)
}

#-- interventionalData ---- #
#-- This function takes a dag G and a list
#INPUT: G = Dag 
#       L = coefficient matrix 
#       interventionTargets= list of intervention Targets, this
#                           is a list of vectors c(n_I,I) where n_I is the
#                           sample size of the experiment I and I is the set
#                           of intervened nodes.
#       kindOfIntervention = String. Denotes the kind of intervention as in 
#               Yang et al. 2018. Possible strings: "random"=random coefficient,
#               "perfect"= coefficient set to 0, "imperfect"=coefficient set to 
#               0 with probability 0.5, "inhibitory"=coefficient*0.1
#OUTPUT: A list with five sublist,
#       Covs= This is a list of covariance matrices
#              obtained from the interventional experiments.
#       Rs = This is a list of correlation matrices
#              obtained from the interventional experiments.
#       Ls = This is a list of coefficient matrices with the values
#              of the intervened coefficients
#       Ns = This is a list of the sample size of each intervention
#       Xs = The complete data set for all interventional experiments.      
# 

interventionalData<-function(G,L,interventionTargets,kindOfIntervention="random"){
  p<-nrow(L)
  varvector<-abs(rnorm(p,mean=0,sd=100))
  ndatasets<-length(interventionTargets)
  Covlist<-list()
  Rlist<-list()
  Llist<-list()
  Xlist<-list()
  nList<-c()
  edgesIntervened <- 0
  timeCovCor = 0
  for(I in interventionTargets){
    nI<-I[1] # first element is the sample size
    edgIntI <- 0
    if(length(I)==1){ # in this case there are no intervention targets.
      LI<-L
    } else{	
      tI<-I[-1] # remaining elements are the intervention targets
      LI<-matrix(0,p,p)
      LI<-L
      for (j in (1:length(tI))) { # this changes the coeff in the structural eqns
        parentsj<-neighbors(G,tI[j],mode = "in")
        edgIntI <- edgIntI + length(parentsj)
        if (length(parentsj)>0){
          for (k in parentsj){
            
            # parse kindOfIntervention
            if(kindOfIntervention == "random"){
              coeff <- ruunif(1,0.3,1)
            } else if(kindOfIntervention == "perfect"){
              coeff <- 0
            } else if(kindOfIntervention == "imperfect"){
              if(rbinom(1,1,0.5)==1) {
                coeff <- 0
              } else {
                coeff <- 0
              }
            } else if(kindOfIntervention == "inhibitory"){
              coeff <- LI[k,tI[j]] * 0.1
            } else {
              stop("Invalid value of kindOfIntervention.")
            }
            LI[k,tI[j]]<-coeff
          }
        }
      }
    }
    XI<-samplingDAG(nI,LI,varvector)
    s = Sys.time()
    CovI<-sample.cov(XI)
    RI<-cov2cor(CovI)
    timeCovCor = timeCovCor + as.numeric(Sys.time() - s, units="secs")
    Covlist<-append(Covlist,list(CovI))
    Rlist<-append(Rlist,list(RI))
    Llist<-append(Llist,list(LI))
    nList<-append(nList,nI)
    Xlist<-append(Xlist,list(XI))
    edgesIntervened <- edgesIntervened + nI * edgIntI
  }
  RLlist<-list(Rs=Rlist,Ls=Llist,Ns=nList,Xs=Xlist,Covs=Covlist,edgesIntervened=edgesIntervened,timeCovCor=timeCovCor)
  return(RLlist)
}
#Imatrix
# INPUT: corrsIs = Ordered List of observed correlation matrices
#                 associated to the interventional experiments
#        nIs = Ordeded List of sample sizes associated to each
#              interventional experiment.
# OUTPUT: C=a matrix, each entry of the matrix is the mutual information between the two variables
Imatrix<-function(C_list,Ns){
  p<-dim(C_list[[1]])[1]
  C<-matrix(0,nrow=p,ncol=p)
  for(i in c(1:(p-1))){
    for(j in c((i+1):p)){
      for(k in c(1:length(Ns))){
        C[i,j]<-C[i,j]-(Ns[k]/2)*log(1-C_list[[k]][i,j]^2)
      }
    }
  }
  C<-C+t(C)
  return(C)
}


#-wmedianCorrels
# INPUT: corrsIs = Ordered List of observed correlation matrices
#                 associated to the interventional experiments
#        nIs = Ordeded List of sample sizes associated to each
#              interventional experiment.
# OUTPUT: A list consisting of:
#       Rmedian = a matrix, each entry is the weighted median
#       probs = normalized vector of weights for each interventional setting.
#               
# Consolidate function
# The function computes the entrywise weighted median of the correlation 
# entries. The function weighted.median takes a vector of values
# and a vector of weights of the same length.
wmedianCorrels<-function(corrIs,nIs){
  p<-nrow(corrIs[[1]])
  k<-length(nIs)
  probs<- 1/sum(nIs)*nIs
  threewayT<-array(abs(unlist(corrIs)),c(p,p,k))
  Rmedian<-apply(threewayT,1:2,weighted.median,probs)
  return(list(Rmedian=Rmedian,probs=probs))
}



#---- Weighted Mean correlation matrix
#----
#-wmeanCorrels
# INPUT: corrsIs = Ordered List of observed correlation matrices
#                 associated to the interventional experiments
#        nIs = Ordeded List of sample sizes associated to each
#              interventional experiment.
# OUTPUT: A list consisting of:
#       Rmean = a matrix, each entry is the weighted median
#       probs = normalized vector of weights for each interventional setting.
#               
# Consolidate function
# The function computes the entrywise weighted median of the correlation 
# entries. The function weighted.median takes a vector of values
# and a vector of weights of the same length.
wmeanCorrels<-function(corrIs,nIs){
  p<-nrow(corrIs[[1]])
  k<-length(nIs)
  probs<- 1/sum(nIs)*nIs
  threewayT<-array(abs(unlist(corrIs)),c(p,p,k)) # this one had a mistake, now it is fixed
  Rmean<-apply(threewayT,1:2,weighted.mean,probs) # then we take the mean of that
  return(list(Rmean=Rmean,probs=probs))
}






# The next function does simulations with multiple node interventions
# This function generates a DAG, a list of intervention targets 
# with corresponding sample sizes and a matrix L of coefficients
# INPUT: p = number of nodes in the dag
#
#        ndatasets= number of settings
#
#        interventionsize= either an integer between 1 and p-1 or a vector of length ndatasets containing the 
#        the number of intervened nodes in that specific datasets. If interventionsize is an integer
#        1 observentional setting and (ndatasets-1) interventional settings are with interventionsize
#        randomly chosen intervention nodes each are created
#
#        sdatasets= either c() or a vector of length ndatasets containing the sample sizes for each setting.
#        If sdatasets==c(), the first setting gets (totalsample%/%ndatasets+totalsample%%ndatasets) samples
#        while the others get totalsample%/%ndatasets samples.
#
#        totalSample = sum of the sample sizes of all interventional and obsv. experiments, to be specified
#        only if sdatasets=c().
#
#        ensureDiff= Boolean. If TRUE, ensures that all interventional targets
#           are different. Default: TRUE
#
#        use_dags: Boolean. If TRUE, generate DAGs. If FALSE, generate Polytrees.
#           Default: TRUE.
#     
#         dag_nbh: Numeric. The expected number of neighbours per node if DAGs requested.
# If interventionsize==1 and sdatasets==c(), then you get the same output of "interventionalSetting" 
# where proprI=totalsample/ndatsets
isetting<-function(p,ndatasets,interventionsize,sdatasets,totalsample,
                   ensureDiff=TRUE,use_dags=FALSE, dag_nbh=p/10){
  
  if(use_dags) {
    
    # generate DAGs
    gOrig = randDAG(p, dag_nbh, method="er")
    adj = as(gOrig, "matrix")
    g = NULL
    g$Directed = 1*(abs(adj) > 0)
    g$Skeleton = 1* ((abs(adj) > 0) | t(abs(adj) > 0))
    rownames(g$Directed) = colnames(g$Directed) =
      rownames(g$Skeleton) = colnames(g$Skeleton) = NULL
    lambdaCoeffs<-coeffLambda(graph_from_adjacency_matrix(g$Directed))
  } else {
    
    # generate polytrees
    g<-pruferwithskeletonOpt(p)
    lambdaCoeffs<-coeffLambda(graph_from_adjacency_matrix(g$Directed))
  }


  if(length(sdatasets)==0){
    s<-sdatasets
    sdatasets<-0*c(1:ndatasets)+totalsample%/%ndatasets
    sdatasets[1]<-sdatasets[1]+totalsample%%ndatasets
  }
  if(length(interventionsize)==1){
    s<-interventionsize
    interventionsize<-0*c(1:ndatasets)+interventionsize
    interventionsize[1]<-0
  }
  interventionTargets<-list()
  if(! ensureDiff) {
    # just independently sample interventional sets
    for(i in c(1:ndatasets)){
      I<-sample(p,interventionsize[i],replace=FALSE)
      interventionTargets<-append(interventionTargets,list(c(sdatasets[i],I)))
    }
  } else {
  
    # check that generating different targets is possible
    warning("There might be an infinite loop. Please be sure that generating different intervention targets is possible.")
    
    # naive implementation with re-sampling (might be slow)
    intvTargetSets = vector(mode="list", ndatasets-1)
    for(i in (1:ndatasets)){
      repeat{
        I <- as.set(sample(p,interventionsize[i],replace=FALSE))
        j=1; found=FALSE;
        while(j<i && found==FALSE){
          found = intvTargetSets[[j]] == I
          j=j+1
        }
        if(found==FALSE){
          break
        }
      }
      intvTargetSets[[i]] = I
    }
    
    # parse into output list
    interventionTargets = vector(mode="list",ndatasets)
    interventionTargets[[1]] = c(ndatasets[1])
    for(i in 1:ndatasets){
      interventionTargets[[i]] = c(sdatasets[i],as.numeric(intvTargetSets[[i]]))
    }
  }

  return(list(gTrued= g$Directed, gTrues=g$Skeleton, L=lambdaCoeffs, targetsI=interventionTargets))
}


#-----------
#Functions to learn the CPDAG

#Functions to learn the CPDAG

#Performs procedure 1 of section 5.2
#INPUT:   Covlist= list of covariance matrices
#         Ilist= list of intervened nodes
#         Nlist= list of sample sizes
#         E= list of unoriented edges
#         C= matrix of dependece measures
#         thres= threshold for independence test
#         method= "simple" if you want to replace the collider test in 5.2 with the one in 5.1
#OUTPUT:  oriented=  list of oriented edges
#         unoriented= list of unoriented edges
complete_triplet<-function(p,Covlist,Ilist,Nlist,E,C,thres,method,pw_method,alpha,Xlist){
  O<-numeric()
  Trip_out<-triplets(Covlist,Ilist,Nlist,E,O,p,C,thres,method)
  U<-Trip_out$Ulist
  O<-Trip_out$Olist
  if(pw_method=="BIC"){
    repeat{
      if(length(U)==0){
        break
      }
      d<-dir_ident_edges(Ilist,U)
      v_list<-dir_ident_vert(U,d)
      if(sum(d)==0){
        break
      }
      Icomp<-i_component(matrix(U,ncol=2),v_list)
      TBO<-Icomp$tboriented
      U<-Icomp$unoriented
      
      d<-dir_ident_edges(Ilist,TBO)
      IDvertices<-dir_ident_vert(TBO,d)
      
      
      P_List<-aftercollider_test(Covlist,Ilist,Nlist,TBO,IDvertices)
      O<-rbind(O,P_List$o_edges)
      
    }
    return(list(oriented=O,unoriented=U))
  }else{
    if(length(U)==0){
      return(list(oriented=O,unoriented=U))
    }
    out<-test_aftercoll(Ilist,Nlist,U,p,alpha,Xlist)
    U<-out$unoriented
    O<-rbind(O,out$oriented)
    return(list(oriented=O,unoriented=U))
  }
}

test_aftercoll<-function(Ilist,Nlist,E,p,alpha,Xlist){
  TOUSE<-numeric()
  O<-numeric()
  not_orientable<-list()
  repeat{
    if(length(TOUSE)==0){
      E<-matrix(E,ncol=2)
      index<-first_dir_id_edge(Ilist,E,not_orientable)
      if(index!=0){
        o<-as(pairs(matrix(E[index,],ncol=2),Xlist,Ilist,alpha)$Olist,"numeric")
        if(length(o)>0){
          O<-rbind(O,o)
        }
        else{
          not_orientable<-rbind(not_orientable,E[index,])
        }
      }else{
        return(list(unoriented=E,oriented=O))
        break
      }
      if(length(o)>0){
        E<-E[-index,]
        E<-matrix(E,ncol=2)
        if(length(E)==0){
          return(list(unoriented=E,oriented=O))
          break
        }
      }
    }
    else{
      TOUSE<-matrix(TOUSE,ncol=2)
      o<-TOUSE[1,]
      TOUSE<-TOUSE[-1,]
    }
    if(length(o)>0){
      if(length(E)!=0){
        E<-matrix(E,ncol=2)
        l1<-which(E[,1]==o[2])
        l2<-which(E[,2]==o[2])
        for(i in l2){
          E[i,]<-rev(E[i,])
        }
        l<-append(l1,l2)
        E_o<-matrix(E[l,],ncol=2)
        if(length(E_o)==0){
        }
        else{
          E<-E[-l,]
          for(i in c(1:(length(E_o)/2))){
            O<-rbind(O,E_o[i,])
            TOUSE<-rbind(TOUSE,E_o[i,])
          }
        }
      }
      else{
        return(list(unoriented=E,oriented=O))
        break
      }
    }
  }
}

#Performs the first part of procedure 1 of section 5.2
#INPUT:   Covlist= list of covariance matrices
#         Ilist= list of intervened nodes
#         Nlist= list of sample sizes
#         E= list of unoriented edges
#         C= matrix of dependece measures
#         thres= threshold for independence test 
triplets<-function(Covlist,Ilist,Nlist,E,O,p,C,thres,method){
  V_v<-c(1:p)  #list of vertices to "check" i.e possible colliders
  UE<-E
  vcheck<-numeric()       #list of already checked vertices
  repeat{
    if(length(V_v)==0){
      LIST<-list(Olist=O,Ulist=UE)
      return(LIST)
      break
    }
    update<-FALSE               #needed to break the "repeat" loop when O can't be updated anymore
    for(i in c(1:length(V_v))){
      vupdate<-FALSE            #needed to check if there some edges have been oriented in this step of the for cycle
      vcheck<-numeric()          #list of vertices that can't be colliders
      v<-V_v[i]
      if(length(UE)>2){
        Lu1<-which(UE[,1]==v)    #edges in which v appears as 1st element
        Lu2<-which(UE[,2]==v)    #edges in which v appears as 2nd element
      }
      else{
        Lu1<-which(UE[1]==v)
        Lu2<-which(UE[2]==v)
      }
      L<-append(Lu1,Lu2)
      l<-length(L)
      if(l==1){
        vcheck<-append(vcheck,i)  #if there is only 1 edge containing it, v can't be a collider. 
      }
      if(length(O)>0){
        Li<-which(O[,2]==v)       #check if there is some already oriented edge that points toward v.
      }
      else{
        Li<-numeric()
      }
      if((length(Li)>0)&(l>0)){   #if there are edges pointing towards v, and unoriented edges containig v orient the unoriented edges.
        e<-O[Li[1],]  
        O<-withincoming(Covlist,Ilist,Nlist,e,v,Lu1,Lu2,O,UE,C,thres,method)   #updated list of oriented edges
        vcheck<-append(vcheck,i)
        update<-TRUE
        vupdate<-TRUE
      }
      else{
        if(l>1){                  #if there is more than 1 unoriented edge containing v, check if v is a collider.
          d<-onlyundirected(C,UE,Lu1,Lu2,thres)   #d[i,j]=1 if the two edges i,j forms a collider with v as center
          col<-which(d==1)       #position of the colliders 
          if(length(col)>0){
            a<-col[1]%%l         #position of the first collider
            if(a==0){
              a<-l
            }
            d<-d+t(d)
            d<-d+diag(dim(d)[1])
            
            d<-d[a,]             #edges that collide with the first collider 
            if(a<=length(Lu1)){
              e<-UE[Lu1[a],2]
              O<-withlist(Lu1,Lu2,O,UE,d)   #updated list of oriented edges
            }
            else{
              e<-UE[Lu2[a],1]
              O<-withlist(Lu1,Lu2,O,UE,d)   #updated list of oriented edges
            }
            vcheck<-append(vcheck,i)        
            update<-TRUE
            vupdate<-TRUE
          }
        }
      }
      if(vupdate){
        if(length(UE)>2){
          UE<-UE[-L,]
        }
        else{
          LIST<-list(Olist=O,Ulist=numeric())
          return(LIST)
          break
        }
      }
    }
    if(update==FALSE){
      LIST<-list(Olist=O,Ulist=UE)
      return(LIST)
      break
    }
    V_v<-V_v[-vcheck]        #remove the checked vertices from the list of possible colliders
  }
}








##AUXILIARY FUNCTIONS USED IN TRIPLETS
#return 1 if j is a collider and 0 otherwise
simpleitest<-function(C,i,j,k,thres){
  if(abs(C[i,k])<thres){
    r<-1
  }
  else{
    r<-0
  }
  return(r)
}

#Performs the likelihood test in 5.2 
#return 1 if u is a collider and 0 otherwise
tripletlikelihood<-function(Covlist,Ilist,Nlist,u,v,w){
  a<-tripletlikelihood_coll(Covlist,Ilist,Nlist,u,v,w)
  b<-tripletlikelihood_noncoll(Covlist,Ilist,Nlist,u,v,w)
  if(a>b){
    return(1)
  }else{
    return(0)
  }
}

#takes as input a list the correlation matrix, an oriented edge, a vertex v, the wo lists of edges containing v, 
#the lists of oriented/unoriented edges and a threshold
#gives as output the updated list of oriented edges
withincoming<-function(Covlist,Ilist,Nlist,e,v,Lu1,Lu2,O,UE,C,thres,method){
  d<-matrix(0,nrow = length(Lu1)+length(Lu2),ncol=1) 
  UE<-matrix(UE,ncol=2)
  if(length(Lu1)>0){
    for(j in c(1:length(Lu1))){
      if(method=="simple"){
        d[j]<-simpleitest(C,e[1],v,UE[Lu1[j],2],thres)
      }
      else{
        d[j]<-tripletlikelihood(Covlist,Ilist,Nlist,e[1],v,UE[Lu1[j],2])
      }
    }
  }
  if(length(Lu2)>0){
    for(j in c(1:length(Lu2))){
      if(method=="simple"){
        d[j]<-simpleitest(C,e[1],v,UE[Lu2[j],1],thres)
      }
      else{
        d[length(Lu1)+j]<-tripletlikelihood(Covlist,Ilist,Nlist,e[1],v,UE[Lu2[j],1])
      }
    }
  }
  O<-withlist(Lu1,Lu2,O,UE,d)
  return(O)
}

#takes as input a correlation matrix, a list of undirected egdes, the two lists of edges having v as 1st or 2nd vertex and a threshold
#gives as output a matrix d, with d[i,j]=1 if the two edges i,j forms a collider with v as center
onlyundirected<-function(C,UE,Lu1,Lu2,thres){
  L<-append(Lu1,Lu2)
  l<-length(L)
  d<-matrix(0,nrow =l,ncol=l)
  for(i in c(1:(l-1))){
    for(j in c((i+1):l)){
      if(i<=length(Lu1)){
        count_i<-2
      }
      else{
        count_i<-1
      }
      if(j<=length(Lu1)){
        count_j<-2
      }
      else{
        count_j<-1
      }
      d[i,j]<-simpleitest(C,UE[L[i],count_i],v,UE[L[j],count_j],thres)
    }
  }
  return(d)
}

#takes as input the two lists of edges containig a vertex v, and the lists of oriented/unoriented edges.
#gives as output the updated list of oriented edges.
withlist<-function(Lu1,Lu2,O,UE,d){
  for(j in c(1:length(d))){
    if(j<=length(Lu1)){
      if(d[j]==0){
        if(length(UE)>2){
          O<-rbind(O,UE[Lu1[j],])  
        }
        else{
          O<-rbind(O,UE)  
        }
      }
      else{
        if(length(UE)>2){
          O<-rbind(O,rev(UE[Lu1[j],]))
        }
        else{
          O<-rbind(O,rev(UE))
        }
      }
    }
    else{
      if(d[j]==0){
        if(length(UE)>2){
          O<-rbind(O,rev(UE[Lu2[j-length(Lu1)],]))
        }
        else{
          O<-rbind(O,rev(UE))
        }
      }
      else{
        if(length(UE)>2){
          O<-rbind(O,UE[Lu2[j-length(Lu1)],])
        }
        else{
          O<-rbind(O,UE)
        }
      }
    }
  }
  return(O)
}


# Function that implements strategy 3 in section 5.2
#INPUT:   Covlist= list of covariance matrices
#         Ilist= list of intervened nodes
#         Nlist= list of sample sizes
#         E= list of unoriented edges
#         C= matrix of dependece measures
#         thres= threshold for independence test
#         p= size of the tree
#         method= "simple" if you want to replace the collider test in 5.2 with the one in 5.1
#         pw_method= "BIC" to use 4.1 and "test" for 4.2
#         alpha= significance level for the "test" 4.2
#         Xlist= list of datasets
#OUTPUT:  oriented=  list of oriented edges
#         unoriented= list of unoriented edges


dir_i_or_first<-function(Covlist,Ilist,Nlist,C,thres,E,p,method,pw_method,alpha,Xlist){
  if(pw_method=="BIC"){
    d<-dir_ident_edges(Ilist,E)
    O<-c()
    for(i in c(1:length(d))){
      if(d[i]==1){
        a<-pairlikelihood(Covlist,Ilist,Nlist,E[i,1],E[i,2])
        b<-pairlikelihood(Covlist,Ilist,Nlist,E[i,2],E[i,1])
        if(a>b){
          O<-rbind(O,E[i,])
        }else{
          O<-rbind(O,rev(E[i,]))
        }
      }
    }
    E<-matrix(E,ncol=2)
    E<-E[-which(d==1),]
    E<-matrix(E,ncol=2)
  }else{
    out<-pairs(E,Xlist,Ilist,alpha)
    E<-out$Ulist
    O<-out$Olist
  }
  Trip_out<-triplets(Covlist,Ilist,Nlist,E,O,p,C,thres,method)
  E<-Trip_out$Ulist
  O<-Trip_out$Olist
  return(list(oriented=O,unoriented=E))
}


#Function that implements strategy 2 in section 5.2
#INPUT:   Covlist= list of covariance matrices
#         Ilist= list of intervened nodes
#         Nlist= list of sample sizes
#         E= list of unoriented edges
#         C= matrix of dependece measures
#         thres= threshold for independence test
#         p= size of the tree
#         method= "simple" if you want to do the collider test as in 5.1  
#         pw_method= "BIC" to use 4.1 and "test" for 4.2
#         alpha= significance levele for the "test" 4.2
#         Xlist= list of datasets
#OUTPUT:  oriented=  list of oriented edges
#         unoriented= list of unoriented edges

complete_alternating<-function(Covlist,Ilist,Nlist,C,thres,E,p,method,pw_method,alpha,Xlist){
  alt_out<-alternating(Covlist,Ilist,Nlist,C,E,p,thres,method,pw_method,alpha,Xlist)
  U<-alt_out$unoriented
  O<-alt_out$oriented
  if(length(U)!=0){
    trip_out<-triplets(Covlist,Ilist,Nlist,U,c(),p,C,thres,method)
    U<-trip_out$unoriented
    O<-rbind(O,trip_out$oriented)
  }
  return(list(unoriented=U,oriented=O))
}


alternating<-function(Covlist,Ilist,Nlist,C,E,p,thres,method,pw_method,alpha,Xlist){
  TOUSE<-numeric()
  O<-numeric()
  not_orientable<-list()
  repeat{
    if(length(TOUSE)==0){
      E<-matrix(E,ncol=2)
      index<-first_dir_id_edge(Ilist,E,not_orientable)
      if(index!=0){
        if(pw_method=="BIC"){
          a<-pairlikelihood(Covlist,Ilist,Nlist,E[index,1],E[index,2])
          b<-pairlikelihood(Covlist,Ilist,Nlist,E[index,2],E[index,1])
          if(a>b){
            o<-E[index,]
            O<-rbind(O,o)
          }else{
            o<-rev(E[index,])
            O<-rbind(O,o)
          }
        }else{
          o<-as(pairs(matrix(E[index,],ncol=2),Xlist,Ilist,alpha)$Olist,"numeric")
          if(length(o)>0){
            O<-rbind(O,o)
          }
          else{
            not_orientable<-rbind(not_orientable,E[index,])
          }
        }
      }else{
        return(list(unoriented=E,oriented=O))
        break
      }
      if(length(o)>0){
        E<-E[-index,]
        E<-matrix(E,ncol=2)
        if(length(E)==0){
          return(list(unoriented=E,oriented=O))
          break
        }
      }
    }
    else{
      TOUSE<-matrix(TOUSE,ncol=2)
      o<-TOUSE[1,]
      TOUSE<-TOUSE[-1,]
    }
    if(length(o)>0){
      if(length(E)!=0){
        E<-matrix(E,ncol=2)
        l1<-which(E[,1]==o[2])
        l2<-which(E[,2]==o[2])
        for(i in l2){
          E[i,]<-rev(E[i,])
        }
        l<-append(l1,l2)
        E_o<-matrix(E[l,],ncol=2)
        if(length(E_o)==0){
        }
        else{
          E<-E[-l,]
          for(i in c(1:(length(E_o)/2))){
            if(method=="simple"){
              outcome<-simpleitest(C,o[1],o[2],E_o[i,2],thres)
            }else{
              outcome<-tripletlikelihood(Covlist,Ilist,Nlist,o[1],o[2],E_o[i,2])
            }
            if(outcome==1){
              O<-rbind(O,rev(E_o[i,]))
            }
            else{
              O<-rbind(O,E_o[i,])
              TOUSE<-rbind(TOUSE,E_o[i,])
            }
          }
        }
      }
      else{
        return(list(unoriented=E,oriented=O))
        break
      }
    }
  }
}


first_dir_id_edge<-function(Ilist,E,not_orientable){
  E<-matrix(E,ncol=2)
  if(length(E)==0){
    return(0)
  }
  if(length(not_orientable)==0){
    for(i in c(1:length(E[,1]))){
      for(j in c(1:length(Ilist))){
        if(length(intersect(Ilist[[j]][-1],E[i,]))==1){
          return(i)
          break
        }
      }
    }
    return(0)
  }else{
      for(i in c(1:length(E[,1]))){
        for(j in c(1:length(Ilist))){
          if(length(intersect(Ilist[[j]][-1],E[i,]))==1){
            check<-TRUE
            for(k in c(1:length(not_orientable[,1]))){
              if((E[i,1]==not_orientable[k,1])&&(E[i,2]==not_orientable[k,2])){
                check<-FALSE
              }
            }
            if(check==TRUE){
              return(i)
              break
            }
          }
        }
      }
    return(0)
  }
}
 
 
##Auxiliary Functions

#Weight Matrix
pruferwithskeletonOpt <- function(k){
  undirected = sample_tree(k, directed=FALSE, method="prufer")
  directed = as.directed(undirected, mode="random")
  return(list(Directed=as_adjacency_matrix(directed), 
              Skeleton=as_adjacency_matrix(undirected)))
}


# get edgelist from adjacecncy matrix
edgelist_toadjmatrix<-function(L){
  n<-dim(L)[1]
  I<-Matrix(nrow=(n+1),ncol=(n+1),data=0,sparse=TRUE)
  for(i in c(1:n)){
    I[L[i,1],L[i,2]]<-1
  }
  return(I)
}
#--- Functions to orient edges using regression coefficients ------#
# E = list of edges
# interventionTargets = list of interventional settings
# Environment matrix, Env_mat[i,j]=1 if the intervention j only affect E[i,1], 
# Env_mat[i,j]=-1 if the intervention j only affect E[i,2], Env[i,j]=0 if none of
# the nodes in edge i are targeted in intervention j, Env[i,j]=2 if both E[i,1]
# and E[i,2] are targets in intervention j. In this case we can not use 
# intervention j to test invariance of edge i.
# Env_list is the list of edges orientable using the given interventions
I_env<-function(E,interventionTargets){
  n<-nrow(E)
  l<-length(interventionTargets)
  Env_mat<-matrix(0,n,l) 
  for(i in c(1:n)){
    for(j in c(1:l)){
      if(is.element(E[i,1],interventionTargets[[j]])&&is.element(E[i,2],interventionTargets[[j]])){
        Env_mat[i,j]<-2
      }
      if(is.element(E[i,1],interventionTargets[[j]])&&!is.element(E[i,2],interventionTargets[[j]])){
        Env_mat[i,j]<-1
      }
      if(!is.element(E[i,1],interventionTargets[[j]])&&is.element(E[i,2],interventionTargets[[j]])){
        Env_mat[i,j]<--1
      }
    }

  }
  r_list<-list(Env_matrix=Env_mat)
  return(r_list)
}
#----
#----- Edgewise orientations with pairs
# New version of pairs
# E = matrix. each row is an edge and each column is a vertex of the edge. these 
#     are unoriented edges
# Xs = list of data sets
# Ilist= list of interventions settings with sample sized
# alpha = a significance level
# meth= "min", reject the one with lowest p-value, or accept
#       # right now only min is implemented, there are perhaps other possible ways to do this
pairs<-function(E,Xs,Ilist,alpha,meth="min"){
  O<-c()
  U<-c()
  ndatasets<-length(Ilist)
  if(nrow(E)!=0){
    IE<-I_env(E,Ilist)
    for(i in c(1:nrow(E))){
      e1<-E[i,1]
      e2<-E[i,2] 
      Ie1noe2<-which(IE$Env_matrix[i,]==1) # indices where e1 is a target but e2 is not
      Inoe1e2<-which(IE$Env_matrix[i,]==-1) # indices where e2 is a target but e1 is not
      Xobsv<-cbind(Xs[[1]][,e1],Xs[[1]][,e2]) # This is assuming the first sample is always the observed one
      Xe1noe2 <- vector(mode='list', length=length(Ie1noe2)) # list to save the data sets relevant to test the direction e1->e2
      Xnoe1e2<- vector(mode='list', length=length(Inoe1e2))  # list to save the data sets relevant to test the direction e2->e1
      p_vals_t1<-c()
      p_vals_t2<-c()
      
      if(length(Ie1noe2)>0 & length(Inoe1e2)>0){ # In this case we can test both directions
        for(k in Ie1noe2) {Xe1noe2<-append(Xe1noe2,list(cbind(Xs[[k]][,e1],Xs[[k]][,e2])))} # matters to perform the test e1 is indp var, e2 is dep
        for(k in Inoe1e2) {Xnoe1e2<-append(Xnoe1e2,list(cbind(Xs[[k]][,e2],Xs[[k]][,e1])))} # here e2 is ind, e1 is the dep
        Xe1noe2<-Xe1noe2[-which(sapply(Xe1noe2, is.null))]
        Xnoe1e2<-Xnoe1e2[-which(sapply(Xnoe1e2, is.null))]
        for(i in c(1:length(Xe1noe2))){ # These are the p_values to test e1->e2
          p_val<-F_test(Xobsv,Xe1noe2[[i]])
          p_vals_t1<-cbind(p_vals_t1,p_val)
        }
        for(i in c(1:length(Xnoe1e2))){ # These are the p_values to test e2->e1
          p_val<-F_test(cbind(Xobsv[,2],Xobsv[,1]),Xnoe1e2[[i]]) # The order of the columns is swapped here
          p_vals_t2<-cbind(p_vals_t2,p_val)
        }
        mint1<-min(p_vals_t1)
        mint2<-min(p_vals_t2)
        if(mint1<mint2 & mint1<alpha/length(Ie1noe2)){O<-rbind(O,c(e2,e1))}
        else if (mint2<mint1 & mint2<alpha/length(Inoe1e2) ){O<-rbind(O,c(e1,e2))}
        else{U<-rbind(U,c(e1,e2))}
      }
      else if(length(Ie1noe2)>0 & length(Inoe1e2)==0){ # Can only use the test for e1->e2
        # the order of the two columns matters
        for(k in Ie1noe2) {Xe1noe2<-append(Xe1noe2,list(cbind(Xs[[k]][,e1],Xs[[k]][,e2])))} 
        Xe1noe2<-Xe1noe2[-which(sapply(Xe1noe2, is.null))]
        p_vals_t1<-c()
        for(i in c(1:length(Xe1noe2))){ # These are the p_values to test e1->e2
          p_val<-F_test(Xobsv,Xe1noe2[[i]])
          p_vals_t1<-cbind(p_vals_t1,p_val)
        }
        c<-which(p_vals_t1<  alpha/length(Ie1noe2))
        if (length(c)==0){O<-rbind(O,c(e1,e2))}
        else {O<-rbind(O,c(e2,e1))}
      }
      else if(length(Ie1noe2)==0 & length(Inoe1e2)>0){ # Can only use the test for e2->e1
        for(k in Inoe1e2) {Xnoe1e2<-append(Xnoe1e2,list(cbind(Xs[[k]][,e2],Xs[[k]][,e1])))}
        Xnoe1e2<-Xnoe1e2[-which(sapply(Xnoe1e2, is.null))]
        p_vals_t2<-c()
        for(i in c(1:length(Xnoe1e2))){ # These are the p_values to test e2->e1
          p_val<-F_test(cbind(Xobsv[,2],Xobsv[,1]),Xnoe1e2[[i]])
          p_vals_t2<-cbind(p_vals_t2,p_val)
        }
        c<-which(p_vals_t2<alpha/length(Inoe1e2))
        if (length(c)==0){O<-rbind(O,c(e2,e1))}
        else {O<-rbind(O,c(e1,e2))}

      }
      else{
        U<-rbind(U,c(e1,e2))
        }
    }
  }
  return(list(Olist=O,Ulist=U))
}
#---------------------------------------------------------------------
# Auxiliary functions for the test of equality of regression coefficients
#------- Drops the first element of a list
dropFirst<-function(x){
  x<-x[-1]
  return(x)
}

# ---- Modifying regCoeff
regCoeff<-function(x,y){
  x<-as.numeric(x)
  y<-as.numeric(y)
  xbar<-mean(x)
  ybar<-mean(y)
  sxx<-sum((x-xbar)^2)
  syy<-sum((y-ybar)^2)
  sxy<-sum((x-xbar)*(y-ybar))
  b<-(1/sxx)*sxy
  a<-ybar-b*xbar
  return(c(a,b))
}
#----------

#---- F_test -- Based on Chow 1960- Test of regression coefficients with modified dfs
#------
# If F_test receives Xe1e2 it means we test e1->e2 i.e if e1 is a causal predictor for e2
#
# INPUT: This function receives in the first entry an observational data set, 
# in the second entry a list of data sets in which e1 is intervened on To test the
# hypothesis that the first column c1 is a parent of the second column c2 using comparison
# of regression coefficients. The output is a list of pvalues, one for each data set
# in the secon list
# Test c1->c2. 
# OUTPUT: A list of p-values

F_test<-function(Xobsv,Xint){
  y1<-Xobsv[,2]
  x1<-Xobsv[,1]
  y2<-Xint[,2]
  x2<-Xint[,1]
  x<-append(x1,x2)
  y<-append(y1,y2)
  n1<-nrow(Xobsv)
  n2<-nrow(Xint)
  cx<-regCoeff(x,y)
  cx1<-regCoeff(x1,y1)
  cx2<-regCoeff(x2,y2)
  S1<-sum((y-cx[1]-cx[2]*x)^2)
  S2<-sum((y1-cx1[1]-cx1[2]*x1)^2)
  S3<-sum((y2-cx2[1]-cx2[2]*x2)^2)
  s1hatsq<-(1/(n1-2))*S2
  s2hatsq<-(1/(n2-2))*S3
  f<-((n1-2)*s1hatsq+(n2-2)*s2hatsq)^2/((n1-2)*s1hatsq^2+(n2-2)*s2hatsq^2)
  den<-(S1-S2-S3)/2
  num<-(S2+S3)/f
  FFmod<-den/num
  df1<-2
  p_val<-pf(FFmod,df1,f,lower.tail = FALSE)
  return(p_val)
}

##Function that computes the i-cpdag
#INPUT:   Ilist= list of interventional settings
#         A= adjacency matrix
#OUPUT:   adjacency matrix of the i-cpdag

i_cpdag<-function(Ilist,A){
  check<-FALSE
  l_int<-length(Ilist)
  p<-length(A[,1])
  
  for(i in c(1:l_int)){
    if(length(Ilist[[i]][-1])==0){
      check<-TRUE
    }
  }
  
  if(check==TRUE){
    A_aug<-matrix(0,nrow=(p+l_int),ncol = (p+l_int))
    for(i in c(1:p)){
      for(j in c(1:p)){
        A_aug[i,j]<-A[i,j]
      }
    }
    for(i in c(1:l_int)){
      A_aug[p+i,Ilist[[i]][-1]]<-1
    }
    G_aug<-graph_from_adjacency_matrix(A_aug)
    V(G_aug)$name<-as.character(c(1:(p+l_int)))
    BN_aug<-cpdag(as.bn(G_aug))
    CPDAG_aug<-as.igraph(BN_aug)
    CPDAG<-induced_subgraph(CPDAG_aug,c(1:p))
    return(get.adjacency(CPDAG))
  }
  else{
    G<-graph_from_adjacency_matrix(A)
    V(G)$name<-as.character(c(1:p))
    BN<-cpdag(as.bn(G))
    CPDAG<-as.igraph(BN)
    A_aug<-get.adjacency(CPDAG)
    
    for(i in c(1:l_int)){
      Jlist<-list()
      for(j in c(1:l_int)){
        if(j!=i){
          Jlist[[j]]<-union(Ilist[[i]][-1],Ilist[[j]][-1])
        }else{
          Jlist[[j]]<-c()
        }
      }
      A_i_aug<-matrix(0,nrow=(p+l_int),ncol = (p+l_int))
      for(j in c(1:p)){
        for(k in c(1:p)){
          A_i_aug[j,k]<-A[j,k]
        }
      }
      for(j in c(1:l_int)){
        if(j!=i){
          A_i_aug[p+j,Jlist[[j]]]<-1
        }
      }
      
      G_i_aug<-graph_from_adjacency_matrix(A_i_aug)
      V(G_i_aug)$name<-as.character(c(1:(p+l_int)))
      BN_i_aug<-as.igraph(cpdag(as.bn(G_i_aug)))
      
      CPDAG_i_aug<-induced_subgraph(BN_i_aug,c(1:p))
      
      A_i_aug<-get.adjacency(CPDAG_i_aug)
      A_aug<-pmin(A_i_aug,A_aug)
    }
    return(A_aug)
  }
}

#----------------------------------
##<<<<<<< HEAD

##=======
##Takes as input the the lists of oriented/unoriented edges and the size of the tree
#gives as output the adjacency matrix of the cpdag
cpdag_from_lists<-function(Olist,Ulist,p){
  A<-Matrix(matrix(0,p,p),sparse=TRUE)
  o_n<-length(Olist)/2
  u_n<-length(Ulist)/2
  if(o_n>0){
    Olist<-matrix(Olist,nrow=o_n)
    for(i in c(1:o_n)){
      A[Olist[i,1],Olist[i,2]]<-1
    }
  }
  if(u_n>0){
    Ulist<-matrix(Ulist,nrow=u_n)
    for(i in c(1:u_n)){
      A[Ulist[i,1],Ulist[i,2]]<-A[Ulist[i,2],Ulist[i,1]]<-1
    }
  }
  return(A)
}

