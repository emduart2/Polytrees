##LIBRARIES
library(igraph)
library(Matrix)

##CPDAG

#takes as input a list of unorieted edges, a correlation matrix and a threshold 
#gives as output two lists of orietned/unoriented  edges  
triplets<-function(E,S,thres){
  m<-(nrow(E)+1)
  V_v<-c(1:(m+1))  #list of vertices to "check" i.e possible colliders
  UE<-E
  O<-numeric()       
  #c<-matrix(0,nrow = m,ncol=1)
  vcheck<-numeric()       #list of already checked vertices
  repeat{
    if(length(V_v)==0){
      LIST<-list(Olist=O,Ulist=UE)
      return(LIST)
      break
    }
    update<-FALSE               #needed to break the "repeat" loop when O can't be updated anymore
    for(i in c(1:length(V_v))){
      vupdate<-FALSE            #needed to check if there some edges have been oriented in this step of the for cycle
      vcheck<-numeric()          #list of vertices that can't be colliders
      v<-V_v[i]
      if(length(UE)>2){
        Lu1<-which(UE[,1]==v)    #edges in which v appears as 1st element
        Lu2<-which(UE[,2]==v)    #edges in which v appears as 2nd element
      }
      else{
        Lu1<-which(UE[1]==v)
        Lu2<-which(UE[2]==v)
      }
      L<-append(Lu1,Lu2)
      l<-length(L)
      if(l==1){
        vcheck<-append(vcheck,i)  #if there is only 1 edge containing it, v can't be a collider. 
      }
      if(length(O)>0){
        Li<-which(O[,2]==v)       #check if there is some already oriented edge that points toward v.
      }
      else{
        Li<-numeric()
      }
      if((length(Li)>0)&(l>0)){   #if there are edges pointing towards v, and unoriented edges containig v orient the unoriented edges.
        e<-O[Li[1],]  
        O<-withincoming(S,e,v,Lu1,Lu2,O,UE,thres)   #updated list of oriented edges
        vcheck<-append(vcheck,i)
        update<-TRUE
        vupdate<-TRUE
      }
      else{
        if(l>1){                  #if there is more than 1 unoriented edge containing v, check if v is a collider.
          d<-onlyundirected(S,UE,Lu1,Lu2,thres)   #d[i,j]=1 if the two edges i,j forms a collider with v as center
          col<-which(d==1)       #position of the colliders 
          if(length(col)>0){
            a<-col[1]%%l         #position of the first collider
            if(a==0){
              a<-l
            }
            d<-d+t(d)
            d<-d+diag(dim(d)[1])
            
            d<-d[a,]             #edges that collide with the first collider 
            if(a<=length(Lu1)){
              e<-UE[Lu1[a],2]
              O<-withlist(Lu1,Lu2,O,UE,d)   #updated list of oriented edges
            }
            else{
              e<-UE[Lu2[a],1]
              O<-withlist(Lu1,Lu2,O,UE,d)   #updated list of oriented edges
            }
            vcheck<-append(vcheck,i)        
            update<-TRUE
            vupdate<-TRUE
          }
        }
      }
      if(vupdate){
        if(length(UE)>2){
          UE<-UE[-L,]
        }
        else{
          LIST<-list(Olist=O,Ulist=numeric())
          return(LIST)
          break
        }
      }
    }
    if(update==FALSE){
      LIST<-list(Olist=O,Ulist=UE)
      return(LIST)
      break
    }
    V_v<-V_v[-vcheck]        #remove the checked vertices from the list of possible colliders
  }
}

##AUXILIARY FUNCTIONS USED IN TRIPLETS
#return 1 if j is a collider and 0 otherwise
tripletcomp<-function(S,i,j,k,thres){
  if(abs(S[i,k])<thres){
    r<-1
  }
  else{
    r<-0
  }
  return(r)
}

#takes as input a list the correlation matrix, an orieted edge, a vertex v, the two lists of edges containing v, 
#the lists of oriented/unoriented edges and a threshold
#gives as output the updated list of oriented edges
withincoming<-function(S,e,v,Lu1,Lu2,O,UE,thres){
  d<-matrix(0,nrow = length(Lu1)+length(Lu2),ncol=1) 
  if(length(Lu1)>0){
    for(j in c(1:length(Lu1))){
      if(length(UE)>2){
        d[j]<-tripletcomp(S,e[1],v,UE[Lu1[j],2],thres)
      }
      else{
        d[j]<-tripletcomp(S,e[1],v,UE[2],thres)
      }
    }
  }
  if(length(Lu2)>0){
    for(j in c(1:length(Lu2))){
      if(length(UE)>2){
        d[length(Lu1)+j]<-tripletcomp(S,e[1],v,UE[Lu2[j],1],thres)
      }
      else{
        d[length(Lu1)+j]<-tripletcomp(S,e[1],v,UE[1],thres)
      }
    }
  }
  O<-withlist(Lu1,Lu2,O,UE,d)
  return(O)
}

#takes as input a correlation matrix, a list of undirected egdes, the two lists of edges having v as 1st or 2nd vertex and a threshold
#gives as output a matrix d, with d[i,j]=1 if the two edges i,j forms a collider with v as center
onlyundirected<-function(S,UE,Lu1,Lu2,thres){
  L<-append(Lu1,Lu2)
  l<-length(L)
  d<-matrix(0,nrow =l,ncol=l)
  for(i in c(1:(l-1))){
    for(j in c((i+1):l)){
      if(i<=length(Lu1)){
        count_i<-2
      }
      else{
        count_i<-1
      }
      if(j<=length(Lu1)){
        count_j<-2
      }
      else{
        count_j<-1
      }
      d[i,j]<-tripletcomp(S,UE[L[i],count_i],v,UE[L[j],count_j],thres)
    }
  }
  return(d)
}

#takes as input the two lists of edges containig a vertex v, and the lists of oriented/unoriented edges.
#gives as output the updated list of oriented edges.
withlist<-function(Lu1,Lu2,O,UE,d){
  for(j in c(1:length(d))){
    if(j<=length(Lu1)){
      if(d[j]==0){
        if(length(UE)>2){
          O<-rbind(O,UE[Lu1[j],])  
        }
        else{
          O<-rbind(O,UE)  
        }
      }
      else{
        if(length(UE)>2){
          O<-rbind(O,rev(UE[Lu1[j],]))
        }
        else{
          O<-rbind(O,rev(UE))
        }
      }
    }
    else{
      if(d[j]==0){
        if(length(UE)>2){
          O<-rbind(O,rev(UE[Lu2[j-length(Lu1)],]))
        }
        else{
          O<-rbind(O,rev(UE))
        }
      }
      else{
        if(length(UE)>2){
          O<-rbind(O,UE[Lu2[j-length(Lu1)],])
        }
        else{
          O<-rbind(O,UE)
        }
      }
    }
  }
  return(O)
}



##Orientation with interventions

#Environment matrix, Env:mat[i,j]=1 if the intervention j only affect E[i,1], 
#Env_mat[i,j]=-1 if the intervention j only affect E[i,2], Env[i,j]=0 otherwise
#Env_list is the list of edges orientable using the given interventions

I_env<-function(E,interventionTargets){
  n<-nrow(E)
  l<-length(interventionTargets)
  Env_mat<-matrix(0,n,l)  
  Env_list<-c()
  for(i in c(1:n)){
    c<-FALSE
    for(j in c(1:l)){
      if(is.element(E[i,1],interventionTargets[[j]])&&is.element(E[i,2],interventionTargets[[j]])){
        Env_mat[i,j]<-2
        c<-TRUE
      }
      if(is.element(E[i,1],interventionTargets[[j]])&&!is.element(E[i,2],interventionTargets[[j]])){
        Env_mat[i,j]<-1
        c<-TRUE
      }
      if(!is.element(E[i,1],interventionTargets[[j]])&&is.element(E[i,2],interventionTargets[[j]])){
        Env_mat[i,j]<--1
        c<-TRUE
      }
    }
    if(c==TRUE){
      Env_list<-cbind(Env_list,i)
    }
  }
  r_list<-list(Env_matrix=Env_mat,Env_list=Env_list)
  return(r_list)
}


#Takes as input the two lists of oriented/unoriented edges, a list of correlation matrices, and a list 
#of intervention targets
#gives as output two lists of oriented/unoriented edges
#need to implement an I_test function that takes as input an edge (u,v), the list of interventions 
#in which only u is affected and returns 1 if u->v is the correct orientation, 0 otherwise.
I_orient<-function(E,O,Rlist,interventionTargets){
  Env<-I_env(E,interventionTargets)
  Env_mat<-Env$Env_matrix
  Env_list<-Env$Env_list
  TOUSE<-numeric()
  n<-dim(S)[1]
  if(length(Env_list)==0){
    ret_list(unoriented=E,oriented=0)
    return(ret_list)
  }
  repeat{
    if(length(Env_list)!=0){
      edge_index<-Env_list[1]
      edge<-E[edge_index,]
      edge_env<-which(Env_mat[edge_index,]==1)
      if(length(edge_env)!=0){
        if(I_test[edge,Rlist,edge_env]==1){
          o_edge<-edge
          O<-rbind(O,o_edge)
          TOUSE<-rbind(TOUSE,o_edge)
        }
        else{
          o_edge<-rev(edge)
          O<-rbind(O,O_edge)
          TOUSE<-rbind(TOUSE,o_edge)
        }
      }
      else {
        edge<-rev(edge)
        edge_env<-which(Env_mat[edge_index,]==-1)
        if(I_test[edge,Rlist,edge_env]==1){
          o_edge<-edge
          O<-rbind(O,o_edge)
          TOUSE<-rbind(TOUSE,o_edge)
        }
        else{
          o_edge<-rev(edge)
          O<-rbind(O,o_edge)
          TOUSE<-rbind(TOUSE,o_edge)
        }
      }
      Env_list<-Env_list[-1]
    }
  }
  E<-E[-Env_list,]
  repeat{
    if(length(TOUSE)!=0){
      TOUSE<-matrix(TOUSE,nrow=length(TOUSE)/2)
      o<-TOUSE[1,]
      TOUSE<-TOUSE[-1,]
    
      if(length(E)!=0){
        E<-matrix(E,nrow=length(E)/2)
        l1<-which(E[,1]==o[2])
        l2<-which(E[,2]==o[2])
        for(i in l2){
          E[i,]<-rev(E[i,])
        }
        l<-append(l1,l2)
        E_o<-matrix(E[l,],nrow=length(E[l,])/2)
        if(length(E_o)==0){
        }
        else{
          E<-E[-l,]
          for(i in c(1:(length(E_o)/2))){
            O<-rbind(O,E_o[i,])
            TOUSE<-rbind(TOUSE,E_o[i,])
          }
        }
      }
      else{
        ret_list(unoriented=E,oriented=0)
        return(ret_list)
      }
    }
    else{
      ret_list(unoriented=E,oriented=0)
      return(ret_list)
    }
  }
}