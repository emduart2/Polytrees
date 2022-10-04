nndatasets<-c(2,5,10)
iss<-c(1)
sdatasets<-c()
totalSamples<-c(100,1000,10000)
p<-2
k<-200
SHD_3_simp<-SHD_3_simp_test<-array(0,c(length(nndatasets),length(iss),length(totalSamples),k))
alpha<-0.05

for(ss in c(1:length(totalSamples))){
  for(ii in c(1:length(iss))){
    for(nn in c(1:length(nndatasets))){
      for(kk in c(1:k)){
        
        IS<-isetting(p,nndatasets[nn],iss[ii],sdatasets,totalSamples[ss],ensureDiff=FALSE)
        G<-graph_from_adjacency_matrix(IS$gTrued)
        E<-get.edgelist(G)
        ID<-interventionalData(G,IS$L,IS$targetsI)
        Covlist<-ID$Covs
        Xlist<-ID$Xs
        Ilist<-IS$targetsI
        Nlist<-ID$Ns
        for(i in c(1:length(Covlist))){
          Covlist[[i]]<-Covlist[[i]]*Nlist[i]
        }
        C_list<-ID$Rs
        true_i_cpdag<-dag2essgraph(as_graphnel(G),Ilist)
        
        
        
        thres<-0.5*log(sum(Nlist))*length(Ilist)
        
        lC<-Imatrix(C_list,ID$Ns)
        
        ##Orientation with original skeleton strategy 3
        dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,E,p,method="simple"
                                 ,pw_method="BIC",alpha,Xlist)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_3_simp[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        ##Orientation with original skeleton strategy 3+TEST
        dag_list<-pairs(E,Xlist,Ilist,alpha)
        dag_adj<-cpdag_from_lists(dag_list$Olist,dag_list$Ulist,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        SHD_3_simp_test[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
      }
    }
  }
}

for(ss in c(1:length(totalSamples))){
  for(ii in c(1:length(iss))){
    for(nn in c(1:length(nndatasets))){
      print(paste0("sample size=",totalSamples[ss]," n.of datasets=",nndatasets[nn]))
      print(c(mean(SHD_3_simp[nn,ii,ss,]),var(SHD_3_simp[nn,ii,ss,])))
      print(c(mean(SHD_3_simp_test[nn,ii,ss,]),var(SHD_3_simp_test[nn,ii,ss,])))
      print("--")
    }
  }
}

