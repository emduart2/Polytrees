p<-10
proprI<-0.5
propObsSample<-0.2


k<-100               #number of repetitions
totalSamples<-c(100,200)
distributions<-list("beta","chi_2_p","chi_2_test")
methods<-list("Stouffer","Fisher","Harmonic")
thresholds<-c(0.01,0.05,0.1,0.2)


##These 3 arrays contains the SHD computed in all the experiments for all the combinations of 
#samplesizes/distribution/method/threshold/skeleton

SHD_cpdag_estimated_skeleton<-array(0,c(k,length(thresholds),
    length(distributions),length(methods),length(totalSamples),2)) #The last component indicated if the skeleton 
                                                                  #is computed with wmean or wmedian
SHD_cpdag_original_skeleton<-array(0,c(k,length(thresholds),
    length(distributions),length(methods),length(totalSamples)))  #This vector contains the SHD when triplets is fed with 
                                                                  #the correct skeleton

SHD_skeleton_only<-array(0,c(length(totalSamples),prod(dim(SHD_cpdag_estimated_skeleton)[1:4]),2)) #This vector contains the SHD only for 
                                                                              #the skeleton


#Simluation
for(ss in c(1:length(totalSamples))){
  count_1<-1
  count_2<-1
  for(d in c(1:length(distributions))){
    for(m in c(1:length(methods))){
      for(t in c(1:length(thresholds))){
        for(i in c(1:k)){
          for(j in c(1,2)){
            IS<-interventionalSetting(p,proprI,propObsSample,totalSamples[ss])
            G<-graph_from_adjacency_matrix(IS$gTrued)
            true_cpdag<-dag2cpdag(as(IS$gTrued,"graphNEL"))
            ID<-interventionalData(G,IS$L,IS$targetsI)
            C_list<-ID$Rs
            
            if(j==1){
              medianC<-wmedianCorrels(C_list,ID$Ns)$Rmedian
              CL<-chowLiu(medianC)
              E_e<-get.edgelist(CL)
              estimated_skeleton<-get.adjacency(CL)
              SHD_skeleton_only[ss,count_1,j]<-sum(abs(estimated_skeleton-IS$gTrues))/(4*(p-1))
              count_1<-count_1+1
            }
            if(j==2){
              meanC<-wmeanCorrels(C_list,ID$Ns)$Rmean
              CL<-chowLiu(meanC)
              E_e<-get.edgelist(CL)
              estimated_skeleton<-get.adjacency(CL)
              SHD_skeleton_only[ss,count_2,j]<-sum(abs(estimated_skeleton-IS$gTrues))/(4*(p-1))
              count_2<-count_2+1
            }
            e_s_cpdag_list<-triplets(E_e,C_list,thresholds[t],distributions[d],methods[m],ID$Ns)
            e_s_cpdag<-as(cpdag_from_lists(e_s_cpdag_list$Olist,e_s_cpdag_list$Ulist,p),"graphNEL")
            
            E<-get.edgelist(graph_from_adjacency_matrix(IS$gTrued))
            o_s_cpdag_list<-triplets(E,C_list,thresholds[t],distributions[d],methods[m],ID$Ns)
            o_s_cpdag<-as(cpdag_from_lists(o_s_cpdag_list$Olist,o_s_cpdag_list$Ulist,p),"graphNEL")
              
            
            SHD_cpdag_estimated_skeleton[i,t,d,m,ss,j]<-shd(true_cpdag,e_s_cpdag)/(dim(edgeMatrix(e_s_cpdag))[2]+dim(edgeMatrix(true_cpdag))[2])
            SHD_cpdag_original_skeleton[i,t,d,m,ss]<-0.5*shd(true_cpdag,o_s_cpdag)/(dim(edgeMatrix(o_s_cpdag))[2]+dim(edgeMatrix(true_cpdag))[2])
          }
        }
      }
    }
  }
}


##Boxplots

#SHD for the skeleton only
df_skeleton_only<-data.frame( #Needs to be adapted if the length of "totalSamples" changes
  SHD_skeleton_only[1,,1],
  SHD_skeleton_only[2,,1],
  
  SHD_skeleton_only[1,,2],
  SHD_skeleton_only[2,,2]
)

skeleton_names<-rep(totalSamples,2)
skeleton_col<-c(rep("red3",6),rep("royalblue3",6))

{
  dev.new(width=5,height=4)
  par(mar=c(4.1,5.1,4.1,3.1))
  boxplot(df_skeleton_only,names=skeleton_names,col=skeleton_col,
        main=c("SHD for the skeleton"), ylab="SHD",xlab="Total Sample size")
  legend("bottomright", legend = c("median","mean"),col=c("red3","royalblue3"),bty = "o",pch = 15)

  abline(v=6.5,lty=1, col="black")
}

#Triplets fed with the right skeleton
for(dd in c(1:3)){  
  dev.new(width=5,height=4)
  par(mar=c(4.1,5.1,4.1,3.1))
  par(mfrow=c(2,2))
  for(t in c(1:4)){
    df_original_skeleton<-data.frame(  #Needs to be adapted if the length of "totalSamples" changes
      SHD_cpdag_original_skeleton[,t,dd,1,1],
      SHD_cpdag_original_skeleton[,t,dd,1,2],
      
      
      SHD_cpdag_original_skeleton[,t,dd,2,1],
      SHD_cpdag_original_skeleton[,t,dd,2,2],

      
      
      SHD_cpdag_original_skeleton[,t,dd,3,1],
      SHD_cpdag_original_skeleton[,t,dd,3,2]
      
    )
    
    {

      
      original_names<-rep(totalSamples,3)
      original_col<-c(rep("red2",length(totalSamples)),rep("royalblue2",length(totalSamples)),rep("seagreen2",length(totalSamples)))
    
      boxplot(df_original_skeleton,names=original_names,col=original_col,
            main=paste("Threshold =",thresholds[t],",","Distribution=",distributions[dd]),ylim=c(0,0.35),
            xlabel="Total Sample Size",ylab="SHD")
    
      legend("topright", legend = methods,col=c("red2","royalblue2","seagreen2"),bty = "o",pch = 15,cex=0.5)
      
      for(i in c(0:3)){
        abline(v=0.5+i*6,lty=1, col="black")
      } 
    }
  }
}

#Triplets fed with the estimated skeleton
skeletons=c("Wmedians","Wmean")
for(m in c(1,2)){
  for(dd in c(1:3)){
    dev.new(width=5,height=4)
    par(mar=c(4.1,5.1,4.1,3.1))
    par(mfrow=c(2,2))
    for(t in c(1:4)){
      df<-data.frame(   #Needs to be adapted if the length of "totalSamples" changes
        SHD_cpdag_estimated_skeleton[,t,dd,1,1,m],
        SHD_cpdag_estimated_skeleton[,t,dd,1,2,m],
        
        
        SHD_cpdag_estimated_skeleton[,t,dd,2,1,m],
        SHD_cpdag_estimated_skeleton[,t,dd,2,2,m],
        
        
        SHD_cpdag_estimated_skeleton[,t,dd,3,1,m],
        SHD_cpdag_estimated_skeleton[,t,dd,3,2,m]
        
      )
      
      {
        e_names<-rep(totalSamples,3)
        cols<-c(rep("red2",2),rep("royalblue2",2),rep("seagreen2",2))
        
        boxplot(df,names=e_names,col=cols,
                main=paste("Threshold =",thresholds[t],",","Distribution=",distributions[dd],",","Skeleton=",skeletons[m]),ylim=c(0,0.1),
                xlabel="Total Sample Size",ylab="SHD")
        
        legend("bottomright", legend = methods,col=c("red2","royalblue2","seagreen2"),bty = "o",pch = 15,cex=0.5)
        
        for(i in c(0:3)){
          abline(v=0.5+i*6,lty=1, col="black")
        } 
      }
    }
  }
}
