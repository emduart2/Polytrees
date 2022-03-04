#-----------------------------
# Needed setup

source("http://bioconductor.org/biocLite.R") 
biocLite("RBGL")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RBGL")
BiocManager::install("graph")
BiocManager::install("Rgraphviz")

install.packages("bnlearn")
install.packages("igraph")
install.packages("matrixcalc")
install.packages("BNSL")
install.packages("pcalg")
install.packages("graph")
install.packages("abind")
install.packages("SID")
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
X<-samplingDAG(100,L) # the output is a data set sampled from the structural
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

#-- This function takes a sample size, a dag G and a list
# of intervention targets and returns a list,
# The first entry of the list is a list of correlation matrices
# obtained from the interventional experiments.
# The second entry of the list, is a list of the coefficient matrices
# of the intervention, this is to check which coefficients are being used
# to run the sampling from the DAG.
interventionalData<-function(n,G,L,interventionTargets){
  X<-samplingDAG(n,L)
  Rlist<-list(sample.cor(X)) # Initialize with an observational correlation matrix
  Llist<-list(L)
  for(I in interventionTargets){
    LI<-matrix(0,p,p)
    LI<-L
    if (length(I)>0){
    for (j in (1:length(I))) {
      parentsj<-neighbors(G,I[j],mode = "in")
      if (length(parentsj)>0){
      for (k in parentsj){
        LI[k,I[j]]<-runif(1,-1,1)
      }
      }
    }
    }
    XI<-samplingDAG(n,LI)
    RI<-list(sample.cor(XI))
    Rlist<-append(Rlist,RI)
    Llist<-append(Llist,list(LI))
  }
  return(list(Rlist,Llist))
}
# Example
el<-matrix(c(1,3,2,3,3,4,4,5),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g)
intervExps<-interventionalData(500,g,L,list(c(2),c(3)))
intervExps[[1]][1]
#---

testLearning<-function(id,G,II){
  results<-c()
  L<-coeffLambda(G)
  for (n in c(50,100,500,1000,5000)){
   intervExps<-interventionalData(n,G,L,II)
   k<-length(intervExps[[1]])
   p<-nrow(L)
   threewayT<-array(unlist(intervExps[[1]]),c(p,p,k))
   Rmedian<-apply(threewayT,1:2,median)
   Rmean<-apply(threewayT,1:2,prod)
   Gobsv<-chowLiu(intervExps[[1]][[1]])  # Learn a skeleton from an obsv correlation,
   Gmedian<-chowLiu(Rmedian)           # median correlation
   Gmean<-chowLiu(Rmean)               # mean correlation
   Gu<-graph_from_edgelist(as_edgelist(G),directed = FALSE)
   dif1<-Gu%m%Gmedian
   dif2<-Gmedian%m%Gu
   difMedian<-max(length(as_edgelist(dif1)),length(as_edgelist(dif2)))
   dif1<-Gu%m%Gmean
   dif2<-Gmean%m%Gu
   difMean<-max(length(as_edgelist(dif1)),length(as_edgelist(dif2)))
   dif1<-Gu%m%Gobsv
   dif2<-Gobsv%m%Gu
   difObsv<-max(length(as_edgelist(dif1)),length(as_edgelist(dif2)))
   results<-rbind(results,c(id,n,p,difObsv,difMedian,difMean))
  }
  return(results)
}


#-----
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

#--- Data set of polytrees
results<-c()
for (i in (1:length(polytreeList))){
  newTest<-testLearning(i,polytreeList[[i]],list(c(2)))
  results<-rbind(results,newTest)
}

mydf<-data.frame(results)
names(mydf)<-c("id","n","p","Robsv", "Rmedian","Rmean")
mydf


names(mydf)

mydf %>%
  ggplot(aes(x=n, y=Rmedian))+
         geom_point()
mydf
ggplot(mydf,mapping = aes(x=n))+
  geom_histogram(bins=3)

df <- gather(mydf, event, total, Rmedian:Rmean)
df
plot <- ggplot(df, aes(n, total, fill=event))
plot <- plot + geom_bar(stat = "identity", position = 'dodge')
plot


gl<-testLearning(1,g,list(c(3)))
gl
plot(gl)
plot(g)

list(diag(c(1,1,1)),diag(c(2,3,4)),diag(c(1,1,1)))
intervExps[[1]]
threewayT <- do.call(abind, c(list(diag(c(1,1,-1)),diag(c(2,3,4)),diag(c(-1,1,1))), list(along=3)))
Rmedian<-apply(threewayT,1:2,median)
Rmean<-apply(threewayT,1:2,prod)

Rmedian
Rmean
prod(c(1,2))
?prod
Rmean<-hadamard.prod( intervExps[[1]][[1]],intervExps[[1]][[2]])
gc<-chowLiu(Rmean)
gm<-chowLiu(Rmedian)
plot(gm)
gc<-chowLiu(intervExps[[1]][[2]])
gc<-chowLiu(intervExps[[1]][[3]])
plot(gc)
dif1<-gu%m%gc
dif2<-gc%m%gu
plot(dif1)
plot(dif2)
# Example
el<-matrix(c(1,3,2,3,3,4,4,5,4,6),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
gu<-graph_from_edgelist(el,directed = FALSE)
plot(g)
plot(chowLiu(interventionalData(1000,g,0)))

#-----


plot(g1)
?randomDAG
admat<-randomDAG(4,0.5,c(1:4))
g<-graph_from_adjacency_matrix(admat)
plot(g)
typeof(sample(p,p))
c(1,2,3)
