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
# Example.
X<-cbind(rnorm(100),rnorm(100),rnorm(100))
sample.cor(X)
#----------
# Code for creating a random matrix of coefficients Lambda
# starting from a dag G and its nodes
# G = Dag,
coeffLambda<-function(G){
  p<-gorder(G)
  Lambda<-matrix(0,p,p)
  edg<-as_edgelist(G) # This is the list of directed edges
  lambdaCoeff<-runif(nrow(edg),min=-1,max = 1) #this is a random vector of coefficients
  # Upper triangular matrix with Lambda coefficients
  for(j in (1:nrow(edg))){
    r<-edg[j,1]
    c<-edg[j,2]
    Lambda[r,c]<- lambdaCoeff[j]
  }
  return(Lambda)
}
# Example:
el<-matrix(c(1,2,2,3,3,4),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g) # matrix of coefficients according to G
#-----------------
# Code for creating a sample from a dagModel specified by a
# matrix of coefficients Lambda
samplingDAG<-function(n,Lambda){
  p<-nrow(Lambda)
  eps<-c() # sample from the errors
  for( i in (1:p)){
    eps<-cbind(eps,rnorm(n))
  }
  Xmat<-matrix(0,n,p)
  Id<-diag(rep(1,p))
  Lambda <- Lambda+Id
  for( j in (1:p)){
    for(i in (1:p)){
      Xmat[,j]<- Xmat[,j]+Lambda[i,j]*eps[,i] 
    }
  }
  return(Xmat)
}
# Example
el<-matrix(c(1,2,2,3,3,4),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g)
X<-samplingDAG(100,L) 
sample.cor(X)
sample.cov(X)
X
# the output is a data set sampled from the structural
#equations implied by the coefficient matrix L
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
# Example 
el<-matrix(c(1,2,2,3,3,4),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g)
X<-samplingDAG(1000,L)
g1<-chowLiu(sample.cor(X))
plot(g1)
# The problem with the ChowLiu function from the bnlearn package
# is that its input is a dataframe, we want the input to be the
# The correlation matrix, the absolute value of this matrix
# is used as the weight matrix to perform the maximum weight
# spanning tree calculation.
#-------------------------------------------------------------
#-- interventionalData ---- #
#-- This function takes a dag G and a list
#INPUT: G = Dag 
#       L = coefficient matrix 
#       interventionTargets= list of intervention Targets, the 
#                            first entry of this list is the sample size
#                            of the observational setting
#OUTPUT: A list with three sublist,
#       List1= This is a list of correlation matrices
#              obtained from the interventional experiments.
#       List2= This is a list of coefficient matrices with the values
#              of the intervened coeffocients
#       List3= This is a list of the sample size of each intervention
# 
interventionalData<-function(G,L,interventionTargets){
  p<-nrow(L)
  n<-interventionTargets[[1]]
  nList<-c(n)
  interventionTargets<-interventionTargets[-1]
  X<-samplingDAG(n,L)
  Rlist<-list(sample.cor(X)) # Initialize with an observational correlation matrix
  Llist<-list(L)
  for(I in interventionTargets){
    nI<-I[1] # first element is the sample size
    tI<-I[-1] # remaining elements are the intervention targets
    LI<-matrix(0,p,p)
    LI<-L
    if (length(I)>0){
    for (j in (1:length(tI))) {
      parentsj<-neighbors(G,tI[j],mode = "in")
      if (length(parentsj)>0){
      for (k in parentsj){
        LI[k,tI[j]]<-runif(1,-1,1)
      }
      }
    }
    }
    XI<-samplingDAG(nI,LI)
    RI<-list(sample.cor(XI))
    Rlist<-append(Rlist,RI)
    Llist<-append(Llist,list(LI))
    nList<-append(nList,nI)
  }
  RLlist<-list(Rs=Rlist,Ls=Llist,Ns=nList)
  return(RLlist)
}
# Example
el<-matrix(c(1,3,2,3,3,4,4,5),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g)
intervExps<-interventionalData(g,L,list(100,c(100,2),c(300,2,3)))
intervExps$Rs # list of sample correlations
intervExps$Ls # list of actual intervened coefficients
intervExps$Ns # list of sample sizes
sum(intervExps$Ns)

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
  threewayT<-array(unlist(corrIs),c(p,p,k))
  Rmedian<-apply(threewayT,1:2,weighted.median,probs)
  return(list(Rmedian=Rmedian,probs=probs))
}

## Tests for the median correlation matrix with three matrices

m1<-matrix(c(1,1,0,1),nrow=2, byrow = TRUE)
m2<-matrix(c(1,1,0.2,1),nrow = 2, byrow = TRUE)
m3<-matrix(c(1.5,1,0.5,1),nrow = 2, byrow = TRUE)
wmedianCorrels(list(m1,m2,m3),c(50,50,50))
wmeanCorrels(list(m1,m2,m3),c(50,50,10))
## Tests for the median correlation matrix with 
# intrventional experiments.

el<-matrix(c(1,3,2,3,3,4,4,5),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g)
intervExps<-interventionalData(g,L,list(100,c(100,2),c(300,2,3)))
intervExps$Rs # list of sample correlations
intervExps$Ls # list of actual intervened coefficients
intervExps$Ns # list of sample sizes
wmedianCorrels(intervExps$Rs,intervExps$Ns)
##-------------------------------------------

#---- Weighted Mean correlation matrix
#----
#-wmedianCorrels
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
  threewayT<-log(abs(array(unlist(corrIs),c(p,p,k)))) # for this one we take the  logs of the absolute values
  Rmean<-apply(threewayT,1:2,weighted.mean,probs) # then we take the mean of that
  return(list(Rmean=Rmean,probs=probs))
}

# Tests with wmeanCorrels
M1<-wmeanCorrels(intervExps$Rs,intervExps$Ns)
M2<-wmedianCorrels(intervExps$Rs,intervExps$Ns)
M1$Rmean
M2$Rmedian
plot(chowLiu(M1$Rmean))

M2
p<-nrow(intervExps$Rs[[1]])
k<-length(intervExps$Ns)
thway<-array(unlist(intervExps$Rs),c(p,p,k))
thway[1,2,]
intervExps$Rs
thway[,1,]

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
# an observed correlation matrix, an weighted median and a weighted mean
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
plot(g)
gs<-testLearningONE(1,g,L,list(100,c(200,2),c(200,2,3)))
gs$shd1
gs$shd2
gs$shd3
par(mfrow=c(1,1))
plot(g)
plot(gs$Gobsv)
plot(gs$Gmean)
plot(gs$Gmedian)


# The next function does simulations with single node interventions
# This function generates a DAG, a list of intervention targets 
# with corresponding sample sizes and a matrix L of coefficients
# INPUT:  p = number of nodes in the dag
#        proprI = percentage of nodes that will get a single target intervention
#        propObsSample= proportion of the samples in the Observed experiment.
#        totalSample = sum of the sample sizes of all interventional and obsv. experiments
interventionalSetting<-function(p,proprI,propObsSample,totalSample){
  g<-pruferwithskeleton(p)
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

settingI <-interventionalSetting(6,0.3,0.4, 100)
plot(graph_from_adjacency_matrix( settingI$gTrued))
settingI$gTrues
settingI$L
settingI$targetsI
#------------
#----
# Setting up low-dimensional learning experiments
# 25 <= p <= 300
# 1 <= n/p <= 300
pstart<- 5 # start of number of nodes
pend<- 6 # maximum number of nodes
propI<- 0.2 # percentage of nodes to do single interventions on
propObsvSample<- 0.1 # percentage of samples to take from the observational distribution.
x<-0
y<-0
results<-c()
for(p in c(pstart:pend)){
  y<-y+1
  for (n in seq(pstart+50,pstart*pend*50,by=50)){
    x<-x+1
    print(c(p,n))
    oneSetting<-interventionalSetting(p,propI,propObsvSample,n)
    g<-graph_from_adjacency_matrix(oneSetting$gTrued)
    result<- testLearningONE(1,g,oneSetting$L,oneSetting$targetsI)
    print(c(result$shd1,result$shd2,result$shd3))
   
    results<-rbind(results,c(x,n,p,n/p,result$shd1,result$shd2,result$shd3))
  }
}

colnames(results)<-c("id","n","p","n/p","SHDobsv", "SHDmean","SHDmedian")
summary(results)

x<-data.frame(results)
par(mfrow=c(1,1))
boxplot(x$SHDobsv~x$n,
        data=x,
        main="ChowLiu Robsv",
        xlab="Sample Size",
        ylab="SHD",
        col="orange",
        border="brown"
)
boxplot(x$Rmedian~x$n,
        data=x,
        main="ChowLiu Rmedian",
        xlab="Sample Size",
        ylab="Number of Wrong Edges",
        col="orange",
        border="brown"
)
boxplot(x$Rmean~x$n,
        data=x,
        main="ChowLiu Rmean",
        xlab="Sample Size",
        ylab="Number of Wrong Edges",
        col="orange",
        border="brown"
)


g<-graph_from_adjacency_matrix(oneSetting$gTrued)
plot(g)
oneSetting$targetsI
result<- testLearningONE(1,g,oneSetting$L,oneSetting$targetsI)
plot(result$Gmean)
plot(result$Gmedian)
plot(result$Gobsv)
plot()
result
#-----
# This function takes and id=identifier, a graph G and a
# set of intervention targets. 
# for each of the sample sizes 50,100,500,1000,500
# it creates interventional data sampled
# from the same structural equation model.
# Using the list of correlation matrices for the interventional
# experiments, it then does a ChowLiu maximum-weight spanning tree
# to learn a skeleton. It then records the number of edges
# that the algorithm got wrong for each learning method
testLearning<-function(id,G,L,II){
  results<-c()
  for (n in c(50,100,500,700,1000,2000,3000,4000,5000)){
   intervExps<-interventionalData(n,G,L,II)
   k<-length(intervExps[[1]])
   p<-nrow(L)
   threewayT<-array(unlist(intervExps[[1]]),c(p,p,k))
   Rmedian<-apply(threewayT,1:2,median)
   Rmean<-apply(threewayT,1:2,prod)
   Gobsv<-chowLiu(intervExps[[1]][[1]])  # Learn a skeleton from an obsv correlation,
   Gmedian<-chowLiu(Rmedian)           # Learn skeleton from median correlation
   Gmean<-chowLiu(Rmean)               # Learn skeleton from mean correlation
   Gu<-graph_from_edgelist(as_edgelist(G),directed = FALSE)
   dif1<-Gu%m%Gmedian
   dif2<-Gmedian%m%Gu
   difMedian<-max(length(as_edgelist(dif1)),length(as_edgelist(dif2)))
   dif1<-Gu%m%Gmean
   dif2<-Gmean%m%Gu
   difMean<-max(length(as_edgelist(dif1)),length(as_edgelist(dif2)))
   dif1<-Gu%m%Gobsv
   dif2<-Gobsv%m%Gu
   difObsv<-max(length(as_edgelist(dif1)),length(as_edgelist(dif2)))/ (p-1)
   results<-rbind(results,c(id,n,p,difObsv,difMedian,difMean))
  }
  return(results)
}
#------------------


#--------------------
#--- Data set of polytrees
el1<-matrix(c(1,2,2,3),nc=2,byrow = TRUE)
g1<-graph_from_edgelist(el1,directed = TRUE)
el2<-matrix(c(1,2,2,3,3,4),nc=2,byrow = TRUE)
g2<-graph_from_edgelist(el2,directed = TRUE)
el3<-matrix(c(1,3,2,3,3,4),nc=2,byrow = TRUE)
g3<-graph_from_edgelist(el3,directed = TRUE)
el4<-matrix(c(1,4,2,4,3,4),nc=2,byrow = TRUE)
g4<-graph_from_edgelist(el4,directed = TRUE)
el5<-matrix(c(1,2,2,3,3,4,4,5),nc=2,byrow = TRUE) 
g5<-graph_from_edgelist(el5,directed = TRUE)
el6<-matrix(c(1,5,2,3,3,4,3,5),nc=2,byrow = TRUE)
g6<-graph_from_edgelist(el6,directed = TRUE)
el7<-matrix(c(1,2,2,3,3,4,4,5,4,6),nc=2,byrow = TRUE)
g7<-graph_from_edgelist(el7,directed = TRUE)
el8<-matrix(c(1,2,2,3,3,4,4,5,5,6),nc=2,byrow = TRUE)
g8<-graph_from_edgelist(el8,directed = TRUE)
polytreeList<-list(g1,g2,g3,g4,g5,g6,g7,g8)

#---
# List of random polytrees on fixed number of nodes

#--- Data set of polytrees
results<-c()
for (i in (1:length(polytreeList))){
  newTest<-testLearning(i,polytreeList[[i]],list(c(2)))
  results<-rbind(results,newTest)
}




