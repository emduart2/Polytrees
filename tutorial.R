#
# Examples
X<-cbind(rnorm(100),rnorm(100),rnorm(100))
sample.cor(X)

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

# Example: Test the learning algorithm on the different correlation matrices
# Robs, Rmean, Rmeadian
el<-matrix(c(1,2,2,6,3,2,3,5,5,7,7,8,4,5),nc=2,byrow = TRUE)
gtrue<-graph_from_edgelist(el,directed = TRUE)
plot(gtrue)
results<-c()
for (i in (1:100)){
  result<-testLearning(i,gtrue,list(c(2),c(4)))
  results<-rbind(results,result)
}
colnames(results)<-c("id","n","p","Robsv", "Rmedian","Rmean")
summary(results)
results[,50,]
rest50<- data.frame(results) %>% filter(n==50)
rest100<- data.frame(results) %>% filter(n==100)
rest500<- data.frame(results) %>% filter(n==500)
rest1000<-data.frame(results) %>% filter(n==1000)
summary(rest50)
summary(rest100)
summary(rest500)
summary(rest1000)
x<-data.frame(results)
par(mfrow=c(2,2))
plot(x$n,x$Robsv)
plot(x$n,x$Rmedian)
plot(x$n,x$Rmean)
