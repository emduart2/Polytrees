##Orientation testing 

##Orientation testing 

nndatasets<-c(3:5)
iss<-c(3:6)
sdatasets<-c()
totalSamples<-c(100,1000,10000)
p<-20
k<-20
SHD_1<-SHD_2<-SHD_2_simp<-SHD_3<-array(0,c(length(nndatasets),
                                                         length(iss),length(totalSamples),k))
for(ss in c(1:length(totalSamples))){
  for(ii in c(1:length(iss))){
    for(nn in c(1:length(nndatasets))){
      for(kk in c(1:k)){
        IS<-isetting(p,nndatasets[nn],iss[ii],sdatasets,totalSamples[ss])
        G<-graph_from_adjacency_matrix(IS$gTrued)
        ID<-interventionalData(G,IS$L,IS$targetsI)
        Covlist<-ID$Covs
        Xlist<-ID$Xs
        Ilist<-IS$targetsI
        Nlist<-ID$Ns      
        C_list<-ID$Rs
        true_i_cpdag_adj<-i_cpdag(Ilist,IS$gTrued)
        true_i_cpdag<-as(true_i_cpdag_adj,"graphNEL")
        
        E<-get.edgelist(G)

        
        thres<-0.5*log(sum(Nlist))*length(Ilist)
        
        lC<-Imatrix(C_list,ID$Ns)

        meanC<-wmeanCorrels(C_list,ID$Ns)$Rmean
        CL<-chowLiu(meanC)
        E_e<-get.edgelist(CL)

       
        ##Orientation with original skeleton strategy 1
        dag_list<-complete_triplet(p,Covlist,Ilist,Nlist,E_e,lC,thres)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_adj<-i_cpdag(Ilist,dag_adj)
        i_cpdag_graph<-as(i_cpdag_adj,"graphNEL")
        
        SHD_1[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)/(sum(true_i_cpdag_adj)+sum(i_cpdag_adj))
        
        ##Orientation with original skeleton strategy 2
        dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,E_e,p,method="simple")
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-as(dag_adj,"graphNEL")
        
        SHD_2_simp[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)/(sum(true_i_cpdag_adj)+sum(i_cpdag_adj))
        
        dag_list<-complete_alternating(Covlist,Ilist,Nlist,lC,thres,E,p,method="nothing")
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_graph<-as(dag_adj,"graphNEL")
        
        SHD_2[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)/(sum(true_i_cpdag_adj)+sum(i_cpdag_adj))
        
        ##Orientation with original skeleton strategy 3
        dag_list<-dir_i_or_first(Covlist,Ilist,Nlist,lC,thres,E_e,p)
        dag_adj<-cpdag_from_lists(dag_list$oriented,dag_list$unoriented,p)
        i_cpdag_adj<-i_cpdag(Ilist,dag_adj)
        i_cpdag_graph<-as(i_cpdag_adj,"graphNEL")
        
        SHD_3[nn,ii,ss,kk]<-shd(true_i_cpdag,i_cpdag_graph)/(sum(true_i_cpdag_adj)+sum(i_cpdag_adj))
      }
    }
  }
}
