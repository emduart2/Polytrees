# This script contains the code to generate Table 2 of the paper.


# Load Data
library(readxl)
sachs_dir = ""  # needs to point to a folder where the sachsdata is located in a format as described in the appendix of Sachs et al. 2005
file.list <- list.files(sachs_dir, pattern='*.xls', recursive = TRUE)
df.list <- lapply(paste(sachs_dir,file.list,sep=""), read_excel)

# Data preparation: the experimental setting is taken from Wang, Solus, Yang, Uhler neurips 2017
Xlist<-list()
Ilist<-list()
Nlist<-c()
Xlist[[1]]<-rbind(as.matrix(df.list[[1]]),as.matrix(df.list[[7]]))
Nlist[1]<-1755
for(i in c(1:5)){
  Xlist[[i+1]]<-as.matrix(df.list[[i+7]])
  Xlist[[i+1]]<-Xlist[[i+1]]-mean(Xlist[[i+1]])
  Nlist[i+1]<-dim(df.list[[i+7]])[1]
}
Ilist <- list()
Ilist[[1]]<-c(Nlist[1])     # baseline
Ilist[[2]]<-c(Nlist[2],7)   # akt -> akt
Ilist[[3]]<-c(Nlist[3],9)   # g0076 -> PKC
Ilist[[4]]<-c(Nlist[4],4)   # psitect -> PIP2
Ilist[[5]]<-c(Nlist[5],2) # u0126 -> MEK
Ilist[[6]]<-c(Nlist[6],5)   # ly -> PIP3?
p<-11


# Create ground truth
edg<-list()
edg[[1]]<-c("PKA")
edg[[2]]<-c("praf","PKA")
edg[[3]]<-c()
edg[[4]]<-c("PIP3","plcg")
edg[[5]]<-c("plcg")
edg[[6]]<-c("pmek","PKA")
edg[[7]]<-c("p44/42","PKA","PIP3")
edg[[8]]<-c("PKC")
edg[[9]]<-c("PIP2","plcg")
edg[[10]]<-c("PKA","pjnk")
edg[[11]]<-c("PKC","PKA")
trueAD<-matrix(0,nrow=11,ncol=11)
rownames(trueAD)<-colnames(trueAD)<-colnames(Xlist[[1]])
for(i in c(1:11)){
  for(j in c(edg[[i]])){
    trueAD[j,i]<-1
  }
}
trueG<-as(trueAD,"graphNEL")



# define important variables
Clist<-list()
Covlist<-list()
for(i in c(1:6)){
  Covlist[[i]]<-cov(Xlist[[i]])
  Clist[[i]]<-cov2cor(Covlist[[i]])
}
lC<-Imatrix(Clist,Nlist)
thres<-3*log(sum(Nlist))


# our estimations
skel = estimate_skeleton(Clist,Nlist,method="mean")
E <- get.edgelist(skel)
est_p1BIC = estimate_orientations(p,Covlist,Ilist,Nlist,E,lC,thres,0.05,Xlist,
                                procedure="1simp",pw_method="TEST")
est_p3BIC = estimate_orientations(p,Covlist,Ilist,Nlist,E,lC,thres,0.05,Xlist,
                                procedure="3simp",pw_method="TEST")

# Rebane and Pearl only observational
skel = estimate_skeleton(Clist[1],Nlist[1],method="mean")
E <- get.edgelist(skel)
lC_Pearl = Imatrix(Clist[1],Nlist[1])
est_PearlObs = estimate_orientations(p,Covlist[1],Ilist[1],Nlist[1],E,lC_Pearl,thres,0.05,Xlist[1],
                                procedure="1",pw_method="BIC")

# Rebane and Pearl treat all data as observational
totalSmpl = sum(Nlist)
Xall = Reduce(rbind, Xlist,c())
CovAll = sample.cov(Xall)
Covlist_Pearl = list(); Covlist_Pearl[[1]] = CovAll * totalSmpl
Clist_Pearl = list(); Clist_Pearl[[1]] = cov2cor(CovAll)
lC_Pearl = Imatrix(Clist_Pearl,c(totalSmpl))
x = list(); x[[1]] = Xall;
skel = estimate_skeleton(Clist_Pearl,c(totalSamples),method="mean")
E <- get.edgelist(skel)
est_PearlAll = estimate_orientations(p,Covlist_Pearl,list(c(totalSmpl)),c(totalSmpl),E,lC_Pearl,thres,0.05,x,
                            procedure="1",pw_method="BIC")

# GIES
totalSmpl = sum(Nlist)
allSamples = matrix(0, totalSmpl, p)
targetindex = array(0, totalSmpl)
targets = vector("list", length(Ilist))
a = b = 1
for(j in 1:length(Ilist)){
  b <- a + Ilist[[j]][1]-1
  allSamples[a:b,] = Xlist[[j]]
  targetindex[a:b] = j
  targets[[j]] = Ilist[[j]][-1]
  a <- b+1
}
setting_GIES = new("GaussL0penIntScore", 
                   data=allSamples, 
                   targets=targets, 
                   target.index=targetindex)
gies.fit = gies(setting_GIES)
est_gies = as(gies.fit$essgraph,"graphNEL")



# table 2: compare results on the sachsdatasets
scoreFcts = list(pcalg::shd,SID)
scoreFctNames = c("SHD","SID")
estimates = list(est_p1BIC,est_p3BIC,est_gies)
estimates_names = c("p1,IRC","p2,IRC","GIES")
mat = matrix(0, length(estimates), length(scoreFcts))
for(i in 1:length(estimates)){
  for(j in 1:length(scoreFcts)){
    mat[i,j] = scoreFcts[[j]](dag2essgraph(trueG,Ilist), estimates[[i]])
  }
}
rownames(mat) = estimates_names
colnames(mat) = scoreFctNames
mat



