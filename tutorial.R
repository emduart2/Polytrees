#
# Examples 1
X<-cbind(rnorm(100),rnorm(100),rnorm(100))
sample.cor(X)
# Example 2
el<-matrix(c(1,2,2,3,3,4),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g) # matrix of coefficients according to G

# Example 3
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

# Example: Getting a matrix of coefficients for a DAG
# this samples the coefficients from the uniform(-1,1)
el<-matrix(c(1,2,2,3,3,4),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g) # matrix of coefficients according to G

# Example: Sampling from the joint distribution
# of a DAG using a fixed matrix L of coefficients.
el<-matrix(c(1,2,2,3,3,4),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g)
X<-samplingDAG(100,L)

# Example: Learning the skeleton with a ChowLiu
# algorithm that uses the abs value of the 
# correlation as a matrix of weights.
el<-matrix(c(1,2,2,3,3,4),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g)
X<-samplingDAG(1000,L)
g1<-chowLiu(sample.cor(X))
plot(g1)

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

# Example: Create a list of interventional correlation
# matrices, starting from a dag, a list of coefficients, and
# a list of observations.
el<-matrix(c(1,3,2,3,3,4,4,5),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g)
intervExps<-interventionalData(500,g,L,list(c(2),c(3)))
intervExps[[1]] # List of correlation matrices
intervExps[[2]] # List of coefficient matrices from which the sample correlation matrices
# were sampled. In particular these matrices have equal entries except for the entries
# of the edges that point to a node that is an intervention target.

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

#----
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

# Example:-------------------
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
#----------------------------

# Example:
#Setting an interventional sample
settingI <-interventionalSetting(6,0.3,0.4, 100)
plot(graph_from_adjacency_matrix( settingI$gTrued))
settingI$gTrues
settingI$L
settingI$targetsI

# Example: 
#See the results of the learning algorithm using different
# sample sizes in a single DAG
el<-matrix(c(1,3,2,3,3,4,4,5),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g)
result<-testLearning(1,g,L,list(c(3),c(4),c(5)))
colnames(result)<-c("id","n","p","Robsv", "Rmedian","Rmean")
result
# Example: Test the learning algorithm on the different correlation matrices
# Robs, Rmean, Rmedian
# We fix a DAG, sample from its structural equations with two interventions
# 100 times and record the number of edges that the chow Liu got wrong
# 04.03.2022
el<-matrix(c(1,2,2,6,2,3,3,5,5,7,7,8,4,5),nc=2,byrow = TRUE)
el
gtrue<-graph_from_edgelist(el,directed = TRUE)
L<-coeffLambda(gtrue)
L
plot(gtrue)
results<-c()
for (i in (1:1000)){
  result<-testLearning(i,gtrue,L,list(c(2),c(4),c(6)))
  results<-rbind(results,result)
}
colnames(results)<-c("id","n","p","Robsv", "Rmedian","Rmean")
x<-data.frame(results)
par(mfrow=c(1,1))
boxplot(x$Robsv~x$n,
        data=x,
        main="ChowLiu Robsv",
        xlab="Sample Size",
        ylab="Number of Wrong Edges",
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

#--- Trying Daniele's code

G<-pruferwithskeleton(3)
G$Directed
G$Skeleton
gdag<-graph_from_adjacency_matrix(G$Directed, mode = "directed")
plot(gdag)
#---------------------------------------------------------------
# E = list of edges
E= matrix(c(1,2),nc=2,byrow=TRUE)
# S =
S=matrix(c(1,0.3,0.3,1),nc=2, byrow = TRUE)
thres=0.01
triplets(E,S,thres)

# Structural Hamming distance example
g1 <- randomDAG(10, prob = 0.2)
g2 <- randomDAG(10, prob = 0.2)
shd(G1,G2)
plot(G2)
G1<-graph_from_adjacency_matrix(g1,mode="directed")
G2<-graph_from_adjacency_matrix(g2,mode="directed")
#--- Checking out Daniele's code
E<- matrix(c(1,2,2,3,3,4), ncol = 2, byrow=TRUE)
it<-list(list(1),list(2),list(3,4))
envs<-I_env(E,it)
envs$Env_matrix
envs$Env_list
#----
edgeList<-matrix(c(1,2,2,3,3,4),nc=2,byrow = TRUE)
g<-graph_from_edgelist(edgeList,directed = TRUE)
plot(g)
L<-coeffLambda(g) 
XObsv<-samplingDAG(100,L)
sample.cor(XObsv)
typeof(XObsv)
XObsv
nrow(XObsv)
#-- A list of Polytrees -----#
#--------------------
#--- Some examples of polytrees given by list of edges.
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
#--------------
#--- Data set of polytrees
results<-c()
for (i in (1:length(polytreeList))){
  newTest<-testLearning(i,polytreeList[[i]],list(c(2)))
  results<-rbind(results,newTest)
}

gg<-pruferwithskeleton(3)
plotgg
gg$Skeleton
g<-graph_from_adjacency_matrix(gg$Skeleton, directed=FALSE)
plot(g)
source("Intervention_functions.R")
source("polyFunctions.R")
setwd("/Users/ElianaDuarte/Documents/gitHubRepos/Polytrees")
