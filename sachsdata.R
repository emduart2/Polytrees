library(readxl)
#Download Data
sachs_dir = "/home/lenny/Documents/Mathematische Statistik Lehrstuhl/sachs.som.datasets/"
file.list <- list.files(sachs_dir, pattern='*.xls', recursive = TRUE)
df.list <- lapply(paste(sachs_dir,file.list,sep=""), read_excel)



#Data preparation: the experimental setting is taken from Wang, Solus, Yang, Uhler neurips 2017
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


##Ground truth
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




#DAG computation

Clist<-list()
Covlist<-list()
for(i in c(1:6)){
  Covlist[[i]]<-cov(Xlist[[i]])
  Clist[[i]]<-cov2cor(Covlist[[i]])
}

lC<-Imatrix(Clist,Nlist)
meanC<-wmeanCorrels(Clist,Nlist)$Rmean

CL<-chowLiu(meanC)
E_e<-get.edgelist(CL)

thres<-3*log(sum(Nlist))
# e_s_dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,E_e,p,method="simple")
# e_s_dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,E_e,lC,thres,method="simple")
e_s_dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,E_e,p,method="simple",pw_method = "BIC")

e_s_dag_adj<-cpdag_from_lists(e_s_dag_list$oriented,e_s_dag_list$unotiented,p)
colnames(e_s_dag_adj)<-colnames(Xlist[[1]])
G<-graph_from_adjacency_matrix(e_s_dag_adj)
G<-as_graphnel(G)

#Comparison
shd(trueG,G)
SID::structIntervDist(trueG,G)$sidUpperBound

par(mfrow=c(1,2))
plot(G)
plot(trueG)




# estimate with GIES
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


# GIES estimate
gies.fit = gies(setting_GIES)
gies_est = as(gies.fit$essgraph,"graphNEL")

shd(trueG, gies_est)
SID::structIntervDist(trueG,gies_est)$sidUpperBound


# check true and false positives
adj_true = as(trueG,"matrix")
adj_est = as(gies_est,"matrix")
(fp = sum((adj_est==1) & (adj_true==0)))
(tp = sum((adj_est==1) & (adj_true==1)))

par(mfrow=c(1,2))
plot(gies_est)
plot(trueG)

