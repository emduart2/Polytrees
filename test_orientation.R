##Orientation testing 
nndatasets<-c(5)

iss<-c(3)
sdatasets<-c()
totalSamples<-c(1000)
p<-15
k<-20
SHD_1_simp<-SHD_1<-SHD_2<-SHD_2_simp<-SHD_3<-SHD_3_simp<-
  SHD_1_simp_test<-SHD_1_test<-SHD_2_test<-SHD_2_simp_test<-SHD_3_test<-SHD_3_simp_test<-array(0,c(length(nndatasets),
                                           length(iss),length(totalSamples),k))
alpha<-0.05
for(ss in c(1:length(totalSamples))){
  for(ii in c(1:length(iss))){
    for(nn in c(1:length(nndatasets))){
      for(kk in c(1:k)){
        
        IS<-isetting(p,nndatasets[nn],iss[ii],sdatasets,totalSamples[ss])
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
        
        #meanC<-wmeanCorrels(C_list,ID$Ns)$Rmean
        #CL<-chowLiu(meanC)
        #E_e<-get.edgelist(CL)
        
        
        ##Orientation with original skeleton strategy 1+BIC
        dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,E,lC,thres,method="simple"
                                   ,pw_method = "BIC",alpha,Xlist)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_1_simp[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,E,lC,thres,method="nothing"
                                   ,pw_method = "BIC",alpha,Xlist)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_1[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        ##Orientation with original skeleton strategy 1+TEST
        dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,E,lC,thres,method="simple"
                                   ,pw_method = "TEST",alpha,Xlist)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        SHD_1_simp_test[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,E,lC,thres,method="nothing",
                                   pw_method = "TEST",alpha,Xlist)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_1_test[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        ##Orientation with original skeleton strategy 2+BIC
        dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,E,p,method="simple",
                                       pw_method = "BIC",alpha,Xlist)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_2_simp[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,E,p,method="nothing",
                                       pw_method = "BIC",alpha,Xlist)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_2[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        ##Orientation with original skeleton strategy 2+TEST
        dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,E,p,
                                       method="simple",pw_method="Test",alpha,Xlist)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_2_simp_test[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,E,p,method="nothing"
                                       ,pw_method="Test",alpha,Xlist)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_2_test[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        ##Orientation with original skeleton strategy 3
        dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,E,p,method="simple"
                                 ,pw_method="BIC",alpha,Xlist)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_3_simp[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,E,p,method="nothing"
                                 ,pw_method="BIC",alpha,Xlist)
        
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_3[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        ##Orientation with original skeleton strategy 3+TEST
        dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,E,p,method="simple"
                                 ,pw_method="Test",alpha,Xlist)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_3_simp_test[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
        
        dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,E,p,method="nothing"
                                 ,pw_method="Test",alpha,Xlist)
        
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-dag2essgraph(as_graphnel(graph_from_adjacency_matrix(dag_adj)),Ilist)
        
        
        SHD_3_test[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)
      }
    }
  }
}
