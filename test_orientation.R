##Orientation testing 

nndatasets<-c(3:5)
iss<-c(3)
sdatasets<-c()
totalSamples<-c(10000)
p<-100
k<-10
SHD_true_skeleton<-SHD_estimated_skeleton<-array(0,c(length(nndatasets),
                                                         length(iss),length(totalSamples),k))
for(ss in c(1:length(totalSamples))){
  for(ii in c(1:length(iss))){
    for(nn in c(1:length(nndatasets))){
      for(kk in c(1:k)){
        IS<-isetting(p,nndatasets[nn],iss[ii],sdatasets,totalSamples[ss])
        G<-graph_from_adjacency_matrix(IS$gTrued)
        ID<-interventionalData(G,IS$L,IS$targetsI)
        Xlist<-ID$Xs
        Ilist<-IS$targetsI
        Nlist<-ID$Ns      
        C_list<-ID$Rs
        true_i_cpdag_adj<-i_cpdag(Ilist,IS$gTrued)
        true_i_cpdag<-as(true_i_cpdag_adj,"graphNEL")
        
        E<-get.edgelist(G)

        
        thres<-0.5*log(sum(Nlist))
        
        lC<-Imatrix(C_list,ID$Ns)

        meanC<-wmeanCorrels(C_list,ID$Ns)$Rmean
        CL<-chowLiu(meanC)
        E_e<-get.edgelist(CL)
        
        ##Orientation with estimated skeleton
        e_s_dag_list<-complete_triplet(p,Xlist,Ilist,Nlist,E_e,lC,thres)
        e_s_dag_adj<-cpdag_from_lists(e_s_dag_list$oriented,e_s_dag_list$unotiented,p)
        e_s_i_cpdag_adj<-i_cpdag(Ilist,e_s_dag_adj)
        e_s_i_cpdag<-as(e_s_i_cpdag_adj,"graphNEL")
        
        SHD_estimated_skeleton[nn,ii,ss,kk]<-shd(true_i_cpdag,e_s_i_cpdag)/(sum(true_i_cpdag_adj)+sum(e_s_i_cpdag_adj))
        
        
        ##Orientation with original skeleton
        o_s_dag_list<-complete_triplet(p,Xlist,Ilist,Nlist,E,lC,thres)
        o_s_dag_adj<-cpdag_from_lists(o_s_dag_list$oriented,o_s_dag_list$unoriented,p)
        o_s_i_cpdag_adj<-i_cpdag(Ilist,o_s_dag_adj)
        o_s_i_cpdag<-as(e_s_i_cpdag_adj,"graphNEL")
        
        SHD_true_skeleton[nn,ii,ss,kk]<-shd(true_i_cpdag,o_s_i_cpdag)/(sum(true_i_cpdag_adj)+sum(o_s_i_cpdag_adj))
      }
    }
  }
}
