library(readxl)
#Download Data
file.list <- list.files(pattern='*.xls', recursive = TRUE)
df.list <- lapply(file.list, read_excel)



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

Ilist<-list()
Ilist[[1]]<-c(Nlist[1])
Ilist[[2]]<-c(Nlist[2],7)
Ilist[[3]]<-c(Nlist[3],9)
Ilist[[4]]<-c(Nlist[4],4)
Ilist[[5]]<-c(Nlist[5],2)
Ilist[[6]]<-c(Nlist[6],5)
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
#e_s_dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,E_e,p,method="simple")
#e_s_dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,E_e,lC,thres,method="simple")
e_s_dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,E_e,p,method="simple")

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
