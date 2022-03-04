# Learning polytrees from interventions

# Underlying true DAG
# Setting up the number of nodes and sample size
n<-100
p<-4
el<-matrix(c(1,2,2,3,3,4),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
# This part samples from the errors of the SEM
eps<-c()
for( i in (1:p)){
  eps<-cbind(eps,rnorm(n))
}
# Using the errors, we create a data set according to the structural equations
Xmat<-matrix(0,n,p)
Lambda<-matrix(0,p,p)
edg<-as_edgelist(g) # This is the list of directed edges
lambdaCoeff<-runif(nrow(edg),min=-1,max = 1) #this is a random vector of coefficients
# Upper triangular matrix with Lambda coefficients
for(j in (1:nrow(edg))){
  r<-edg[j,1]
  c<-edg[j,2]
  Lambda[r,c]<- lambdaCoeff[j]
}
#----- Write the data set according to the structural equations
Id<-diag(rep(1,p))
Lambda <- Lambda+Id
for( j in (1:p)){
   for(i in (1:p)){
     Xmat[,j]<- Xmat[,j]+Lambda[i,j]*eps[,i] 
   }
}
# The matrix Xmat is the data matrix that samples from the SEM defined by G.
#-----
# Changing an entry of the Lambda matrix above and a column
# for the errors, we can simulate an intervention.
#-----
Rhat<-sample.cor(Xmat) #-- We compute the correlation matrix
# The function kruskal takes a matrix of weights and applies the maximum
# weight spanning tree for the matrix. Our matrix of weights
# is the entrywise absolute value of the correlation matrix
# The output of Kruskal is a list of edges, we can use that list to
# create a graph and plot it.
el1<-kruskal(abs(Rhat))
g1<-graph_from_edgelist(el1,directed=FALSE)
plot(g1)
plot(g)
#--- Comparison with the Chow-Liu
cw<-chow.liu(data.frame(Xmat))
plot(cw)
#---------------
Rhat<-sample.cor(X)

#----
#-- Simulations with the function
n<-500
p<-6
el<-matrix(c(1,3,2,3,3,4,4,5,5,6),nc=2,byrow = TRUE)
g<-graph_from_edgelist(el,directed = TRUE)
plot(g)
L<-coeffLambda(g,p)
X0<-samplingDAG(n,p,L)
L1<-L
L1[2,3]<-runif(1,-1,1)
X1<-samplingDAG(n,p,L1)
Rhat<-sample.cor(X0)
R1hat<-sample.cor(X1)
el1<-kruskal(abs(Rhat))
g1<-graph_from_edgelist(el1,directed=FALSE)
plot(g1)
el2<-kruskal(abs(R1hat))
g2<-graph_from_edgelist(el2,directed=FALSE)
Rmean<-hadamard.prod(Rhat,R1hat)
elr<-kruskal(abs(Rmean))
gr<-graph_from_edgelist(elr,directed = FALSE)
plot(gr)
plot(g2)
plot(g)
Rhat
R1hat
#-----
length(interventionTargets[[3]])
p4<-neighbors(g,4,mode = "in")
p4[2]
#
gu<-graph_from_edgelist(el,directed = FALSE)
dif<-gu%m%g1
plot(dif)
length(as_edgelist(dif))
?difference
c<-chow.liu(data.frame(X))
plot(c)
