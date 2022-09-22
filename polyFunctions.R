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
#install.packages('spatstat')
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
library(SID)
library(tidyverse)  
library(ggplot2)
library(rmarkdown)
library(spatstat)
library(poolr)
library(harmonicmeanp)
library(DescTools)
library(Matrix)
library(graph)
library(mpoly)
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

# Function for calculating the sample correlation matrix
# https://github.com/sqyu/CorDiffViz/tree/master/R
sample.cor <- function(X, Y=NULL){ # Calculate sample correlation matrix
  if (is.null(Y)){
    covhat <- sample.cov(X,Y)
    corhat <- covhat/sqrt(diag(covhat)); corhat <- t(corhat)/sqrt(diag(covhat))
  } else {
    n <- nrow(X)
    if (nrow(Y) != n) stop("X and Y must have the same number of rows.")
    Xr <- X - t(matrix(colMeans(X), ncol(X), n))
    Yr <- Y - t(matrix(colMeans(Y), ncol(Y), n))
    covhat <- t(Xr)%*%Yr/n
    corhat <- t(covhat)/sqrt(colMeans(Yr^2))
    corhat <- t(corhat)/sqrt(colMeans(Xr^2))
  }
  return(corhat)
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
  for(j in (1:nrow(edg))){
    r<-edg[j,1]
    c<-edg[j,2]
    Lambda[r,c]<- lambdaCoeff[j]
  }
  return(Lambda)
}
#-----------------
# Creating a sample from a dagModel specified by a
# matrix of coefficients Lambda
# Inputs:
#   n: Integer. Number of samples.
#   Lambda: Numeric Matrix. Matrix of coefficients
# Outputs:
#   Matrix: n-by-dim(Lambda,1) matrix with samples.
samplingDAG<-function(n,Lambda){
  p<-nrow(Lambda)
  eps <- matrix( rnorm(p*n,mean=0,sd=1), p, n) 
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
  ndatasets<-length(interventionTargets)
  Covlist<-list()
  Rlist<-list()
  Llist<-list()
  Xlist<-list()
  nList<-c()
  for(I in interventionTargets){
    nI<-I[1] # first element is the sample size
    if(length(I)==1){ # in this case there are no intervention targets.
      LI<-L
    } else{
      tI<-I[-1] # remaining elements are the intervention targets
      LI<-matrix(0,p,p)
      LI<-L
      for (j in (1:length(tI))) { # this changes the coeff in the structural eqns
        parentsj<-neighbors(G,tI[j],mode = "in")
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
    XI<-samplingDAG(nI,LI)
    CovI<-sample.cov(XI)
    RI<-cov2cor(CovI)
    Covlist<-append(Covlist,list(CovI))
    Rlist<-append(Rlist,RI)
    Llist<-append(Llist,list(LI))
    nList<-append(nList,nI)
    Xlist<-append(Xlist,list(XI))
  }
  RLlist<-list(Rs=Rlist,Ls=Llist,Ns=nList,Xs=Xlist,Covs=Covlist)
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


#-----
# This function takes and id=identifier, a graph G and a
# set of intervention targets. 
# The set of intervention targets contains sample sizes for each experiment
# and the targeted nodes.
# The output of this function is a collection
# of undirected graphs together with the SHD between the 
# true graph and the one learned by Chow-Liu. 
# Three learning algorithms are considered here using
# Chow-Liu with three different weight matrices, using
# an observed correlation matrix, a weighted median and a weighted mean
# the ONE in the name refers to the fact that this function does the
# learning for one instance of a dag a coefficient matrix and a set of
# intervention targets.
testLearningONE<-function(id,G,L,interventionTargets){
  intervExps<-interventionalData(G,L,interventionTargets)
  R1<-intervExps$Rs[[1]] # Observed correlations
  R2<-wmeanCorrels(intervExps$Rs,intervExps$Ns) # weighted mean Correls
  R3<-wmedianCorrels(intervExps$Rs,intervExps$Ns) # weighted median Correls
  G1<-chowLiu(R1)
  G2<-chowLiu(R2$Rmean)
  G3<-chowLiu(R3$Rmedian)
  I1<- as_adjacency_matrix(G1)
  I2<- as_adjacency_matrix(G2)
  I3<- as_adjacency_matrix(G3)
  gTrue <- graph_from_edgelist( as_edgelist(G), directed = FALSE)
  Itrue <- as_adjacency_matrix(gTrue)
  p <- nrow(R1)
  SHDG1 <- sum(abs((Itrue-I1)))/(2*(p-1)) # compute structural Hamming distance for each learned graph
  SHDG2 <- sum(abs((Itrue-I2)))/(2*(p-1))
  SHDG3 <- sum(abs((Itrue-I3)))/(2*(p-1))
  return(list(Gobsv=G1, Gmean=G2, Gmedian=G3,shd1=SHDG1,shd2=SHDG2,shd3=SHDG3))  
}



# The next function does simulations with single node interventions
# This function generates a DAG, a list of intervention targets 
# with corresponding sample sizes and a matrix L of coefficients
# INPUT:  p = number of nodes in the dag
#        proprI = percentage of nodes that will get a single target intervention
#        propObsSample= proportion of the samples in the Observed experiment.
#        totalSample = sum of the sample sizes of all interventional and obsv. experiments
interventionalSetting<-function(p,proprI,propObsSample,totalSample){
  g<-pruferwithskeletonOpt(p)
  numInterv<- ceiling(p*proprI) # proportion of nodes to intervene on
  tI <- sample(1:p, numInterv, replace=F)
  lambdaCoeffs<-coeffLambda(graph_from_adjacency_matrix(g$Directed))
  nObsv<- propObsSample*totalSample
  nInterv<-(1-propObsSample)*totalSample/numInterv
  interventionTargets<- list(nObsv)
  for(i in tI){
    interventionTargets<-append(interventionTargets,list(c(nInterv,i)))
  }
  return(list(gTrued= g$Directed, gTrues=g$Skeleton, L=lambdaCoeffs, targetsI=interventionTargets))
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
#        checkConservativity= Boolean. If TRUE, ensures that resulting family
#           of interventional targets is conservative. Default: FALSE
#
# If interventionsize==1 and sdatasets==c(), then you get the same output of "interventionalSetting" 
# where proprI=totalsample/ndatsets
isetting<-function(p,ndatasets,interventionsize,sdatasets,totalsample,checkConservativity=FALSE){
  g<-pruferwithskeletonOpt(p)
  lambdaCoeffs<-coeffLambda(graph_from_adjacency_matrix(g$Directed))
  
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
  if(! checkConservativity) {
    # just independently sample interventional sets
    for(i in c(1:ndatasets)){
      I<-sample(p,interventionsize[i],replace=FALSE)
      interventionTargets<-append(interventionTargets,list(c(sdatasets[i],I)))
    }
  } else {

    interventionTargets <- vector("list", ndatasets)
    
    # first sample 
    I <- sample(p,interventionsize[i],replace=FALSE)
    interventionTargets[[1]] <- c(sdatasets[1],I)
    
    # sample where to ensure that all nodes that were intervened in the first sample
    # will not be intervened to guarantee conservativity
    stillNeedsToBeNotIntervened <- I
    notToIntervene <- vector("list", ndatasets)
    maxNumNotToSample <- p-interventionsize
    ind <- 2:ndatasets
    if(sum(maxNumNotToSample[ind]) < length(stillNeedsToBeNotIntervened))
      stop("Conservative family of intervention targets not possible with given parameters.")
    for(node in stillNeedsToBeNotIntervened){
      
      # sample intervention index where node ensured not to be intervened
      j <- sample(x = ind, 1, prob = maxNumNotToSample[ind]/sum(maxNumNotToSample[ind]))
      notToIntervene[[j]] <- append(notToIntervene[[j]], node)
      maxNumNotToSample[j] <- maxNumNotToSample[j] - 1
    }
    
    # sample remaining interventions
    for(i in 2:ndatasets){
      I <- sample(setdiff(1:p,notToIntervene[[i]]), interventionsize[i], replace=FALSE)
      interventionTargets[[i]] <- c(sdatasets[i],I)
    }
  }

  return(list(gTrued= g$Directed, gTrues=g$Skeleton, L=lambdaCoeffs, targetsI=interventionTargets))
}


#-----------
#Functions to learn the CPDAG

#Performs procedure 1 of section 5.2
#INPUT:   Covlist= list of covariance matrices
#         Ilist= list of intervened nodes
#         Nlist= list of sample sizes
#         E= list of unoriented edges
#         C= matrix of dependece measures
#         thres= threshold for independence test 
#OUTPUT:  oriented=  list of oriented edges
#         unoriented= list of unoriented edges
complete_triplet<-function(p,Covlist,Ilist,Nlist,E,C,thres){
  Trip_out<-triplets(Covlist,Ilist,Nlist,E,C,thres)
  U<-Trip_out$Ulist
  O<-Trip_out$Olist
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
}

#Performs the first part of procedure 1 of section 5.2
#INPUT:   Covlist= list of covariance matrices
#         Ilist= list of intervened nodes
#         Nlist= list of sample sizes
#         E= list of unoriented edges
#         C= matrix of dependece measures
#         thres= threshold for independence test 
triplets<-function(Covlist,Ilist,Nlist,E,C,thres){
  m<-(nrow(E)+1)
  V_v<-c(1:m)  #list of vertices to "check" i.e possible colliders
  print(V_v)
  UE<-E
  O<-numeric()       
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
        O<-withincoming(Covlist,Ilist,Nlist,e,v,Lu1,Lu2,O,UE)   #updated list of oriented edges
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
withincoming<-function(Covlist,Ilist,Nlist,e,v,Lu1,Lu2,O,UE){
  d<-matrix(0,nrow = length(Lu1)+length(Lu2),ncol=1) 
  if(length(Lu1)>0){
    for(j in c(1:length(Lu1))){
      if(length(UE)>2){
        d[j]<-tripletlikelihood(Covlist,Ilist,Nlist,e[1],v,UE[Lu1[j],2])
      }
      else{
        d[j]<-tripletlikelihood(Covlist,Ilist,Nlist,e[1],v,UE[2])
      }
    }
  }
  if(length(Lu2)>0){
    for(j in c(1:length(Lu2))){
      if(length(UE)>2){
        d[length(Lu1)+j]<-tripletlikelihood(Covlist,Ilist,Nlist,e[1],v,UE[Lu2[j],1])
      }
      else{
        d[length(Lu1)+j]<-tripletlikelihood(Covlist,Ilist,Nlist,e[1],v,UE[1])
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

##Auxiliary Functions

#Weight Matrix
## 
pruferwithskeleton <- function(k){
  if(k>2){
    P_orig <- ceiling(runif(k-2, min = 0, max = k))
    P <- P_orig
    V_orig <- 1:k
    V <- V_orig
    adj_matrix <- Matrix(matrix(numeric(k * k), ncol = k),sparse=TRUE)
    adj_matrixs <-Matrix(matrix(numeric(k * k), ncol = k),sparse=TRUE) 
    
    for(i in 1:(k-2)){
      complement <- setdiff(V, P)
      v_0 <- min(complement)
      V <- setdiff(V, v_0)
      adj_matrixs[which(V_orig == v_0), which(V_orig == P_orig[i])] <- adj_matrixs[which(V_orig == P_orig[i]), which(V_orig == v_0)]<-1
      
      m<-rbinom(1,1,0.5)
      if(m==0){
        adj_matrix[which(V_orig == v_0), which(V_orig == P_orig[i])] <- 1
      }
      else{
        adj_matrix[which(V_orig == P_orig[i]), which(V_orig == v_0)] <- 1
      }
      P <- P[2:length(P)]
    }
    adj_matrixs[which(V_orig == V[1]), which(V_orig == V[2])]<-adj_matrixs[which(V_orig == V[2]), which(V_orig == V[1])]<-1
    m<-rbinom(1,1,0.5)
    if(m==0){
      adj_matrix[which(V_orig == V[1]), which(V_orig == V[2])] <- 1
    }
    else{
      adj_matrix[which(V_orig == V[2]), which(V_orig == V[1])] <- 1
    }
  }
  else{
    adj_matrixs<-matrix(c(0,1,1,0),nrow=2,byrow=TRUE)
    m<-rbinom(1,1,0.5)
    if(m==0){
      adj_matrix<-matrix(c(0,1,0,0),nrow=2,byrow=TRUE)
    }
    else{
      adj_matrix<-matrix(c(0,0,1,0),nrow=2,byrow=TRUE)
    }
  }
  ret<-list(Directed=adj_matrix,Skeleton=adj_matrixs)
  return(ret)
}

# same as pruferwithskeleton, but with optimized built-in method
pruferwithskeletonOpt <- function(k){
  undirected = sample_tree(k, directed=FALSE, method="prufer")
  directed = as.directed(undirected, mode="random")
  return(list(Directed=as_adjacency_matrix(directed), 
              Skeleton=as_adjacency_matrix(undirected)))
}


#get edgelist from adjacecncy matrix
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
        c<-TRUE
      }
      if(!is.element(E[i,1],interventionTargets[[j]])&&is.element(E[i,2],interventionTargets[[j]])){
        Env_mat[i,j]<--1
        c<-TRUE
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
pairs<-function(E,Xs,Ilist,alpha,meth="min"){
  O<-c()
  U<-c()
  ndatasets<-length(Ilist)
  if(length(E)!=0){
    l<-length(Xs)
    IE<-I_env(E,Ilist)
    for(i in c(1:nrow(E))){
      e1<-E[i,1]
      e2<-E[i,2] 
      Ie1e2<-append(which(IE$Env_matrix[i,]==0),which(IE$Env_matrix[i,]==1)) # indices suitable to test e1 e2
      Ie2e1<-append(which(IE$Env_matrix[i,]==0),which(IE$Env_matrix[i,]==-1)) # indices suitable to test e2 e1
      if(length(Ie1e2)>1 & length(Ie2e1)>1){ # In this case we can test both directions
        Xe1e2 <- vector(mode='list', length=length(Ie1e2)) # list to save the data sets relevant to test the direction e1->e2
        Xe2e1<- vector(mode='list', length=length(Ie2e1))  # list to save the data sets relevant to test the direction e2->e1
        for(k in Ie1e2) {Xe1e2<-append(Xe1e2,list(cbind(Xs[[k]][,e1],Xs[[k]][,e2])))} # the order in of this two columns
        for(k in Ie2e1) {Xe2e1<-append(Xe2e1,list(cbind(Xs[[k]][,e2],Xs[[k]][,e1])))} # matters to perform the test
        Xe1e2<-Xe1e2[-which(sapply(Xe1e2, is.null))]
        Xe2e1<-Xe2e1[-which(sapply(Xe2e1, is.null))]
        o1<- F_test(Xe1e2) # Test e1->e2---
        o2<- F_test(Xe2e1)  # Test e2->e1
        if (meth=="min"){
          if(min(o1)< min(o2)){O<-rbind(O,c(e2,e1))}   # If o1 has the smallest p-value we reject  e1->e1
          else if(min(o2)<min(o1)){O<-rbind(O,c(e1,e2))}
        }
        else if (meth=="max"){
          if(max(o1)< max(o2)){O<-rbind(O,c(e2,e1))}   # If o2 has the hightes p-value we accept  e2->e1
          else if(max(o2)<max(o1)){O<-rbind(O,c(e1,e2))}
        }
      }
      else if(length(Ie1e2)<2 & length(Ie2e1)>1){ # Can only use the test for e2->e1
        Xe2e1<- vector(mode='list', length=length(Ie2e1))
        for(k in Ie2e1) {Xe2e1<-append(Xe2e1,list(cbind(Xs[[k]][,e2],Xs[[k]][,e1])))} # order matters to perform the test
        Xe2e1<-Xe2e1[-which(sapply(Xe2e1, is.null))]
        o2<-F_test(Xe2e1)
        if(length(which(o2<alpha/length(Ie2e1)))==0) O<-rbind(O,c(e2,e1)) else O<-rbind(O,c(e1,e2))
      }
      else if(length(Ie1e2)>1 & length(Ie2e1)<1){
        Xe1e2<- vector(mode='list', length=length(Ie1e2))
        for(k in Ie1e2) {Xe1e2<-append(Xe1e2,list(cbind(Xs[[k]][,e1],Xs[[k]][,e1])))} # order matters to perform the test
        Xe1e2<-Xe1e2[-which(sapply(Xe1e2, is.null))]
        o1<-F_test(Xe1e2)
        if(length(which(o1<alpha/length(Ie1e2)))==0) O<-rbind(O,c(e1,e2)) else O<-rbind(O,c(e2,e1))
      }
      else{U<-rbind(U,c(e1,e2))}
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
#----------
#---- GroupData
# X = a list of data sets
# This makes one single data set from all the X's
groupD<-function(X){
  Xregroup<-c()
  for (i in X){
    Xregroup<-rbind(Xregroup,i)
  }
  return(Xregroup)
}
#--- regCoeff
#-- Quick function to compute the regression coeff. Only works for vectors
# for matrices we need to change to matrix multiplication and inverse of a matrix.
# x= predictor variable y= response variable
regCoeff<-function(x,y){
  xbar<-mean(x)
  ybar<-mean(y)
  sxx<-sum((x-xbar)^2)
  syy<-sum((y-ybar)^2)
  sxy<-sum((x-xbar)*(y-ybar))
  b<-(1/sxx)*sxy
  a<-ybar-b*xbar
  return(c(a,b))
}
#---------------

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

#---- F_test -- Based on Chow's test 1960
#------
# If F_test receives Xe1e2 it means we test e1->e2 i.e if e1 is a causal predictor for e2
# If F_test receives Xe2e1 it means we test e2->e1
# INPUT: This function receives a list of data sets with two columns. To test the
# hypothesis that the first column c1 is a parent of the second column c2 using comparison
# of regression coefficients.
# Test c1->c2. 
# OUTPUT: A list of p-values
F_test<-function(Xlist){
  n<-length(Xlist)
  pvals<-c()
  for (i in c(1:length(Xlist))){
    XnotI<-groupD(Xlist[-i]) # data set that doesnt have the i-th sample
    XI<-Xlist[[i]] # data set for the sample I
    Xe<-matrix(cbind(rep(1,nrow(XI)),XI[,1]),nrow=nrow(XI),ncol=2)
    Xenot<-matrix(cbind(rep(1,nrow(XnotI)),XnotI[,1]),nrow=nrow(XnotI),ncol=2)
    Ye<-matrix(cbind(rep(1,nrow(XI)),XI[,2]),nrow=nrow(XI),ncol=2)
    Yenot<-matrix(cbind(rep(1,nrow(XnotI)),XnotI[,2]),nrow=nrow(XnotI),ncol=2)
    regCoeffNotI<-regCoeff(Xenot[,2],Yenot[,2]) # OLS estimator computed on Ie
    #print(regCoeffNotI)
    D<-matrix(Ye-(regCoeffNotI[1]+regCoeffNotI[2]*Xe),nrow = nrow(Ye),ncol=1)
    ne<-nrow(Xe)
    nenot<-nrow(Xenot)
    SigmaD<- Xe%*%solve(t(Xenot)%*%Xenot)%*%t(Xe) +diag(1,nrow=ne) #+  # Xe%*%solve(t(Xenot)%*%Xenot)%*%t(Xe)
    sigmahatYnote<-var(Yenot[,2]-(regCoeffNotI[1]+regCoeffNotI[2]*Xenot[,2]))
    Q<-t(D)%*%SigmaD%*%D/(sigmahatYnote*nrow(Xe))
    p_val<-pf(Q,ne,nenot-2, lower.tail=FALSE,log.p = FALSE)
    #print("This is a pvalue")
    #print(p_val)
    pvals<-cbind(pvals,c(p_val))
    #print(c(nrow(D),nrow(SigmaD),nrow(Xe),nrow(XI),Q,sigmahatXnote,p_val))
    #if (p_val< bfcorrection){
    #  return(0) # reject c1->c2
    #}
  }
  return(pvals) # do not reject c2->c1
}
#--------- end of F_test
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

