#install.packages("mpoly")
library("mpoly")

#Computes the maximal marginal log-likelihood of the model u->v
#as in section 4.1 of the paper
#INPUT:   Xlist= list of samples
#         Ilist= list of intervened nodes
#         Nlist= list of sample sizes
#         u,v= indeces of the nodes
#OUTPUT:  Maximal log-likelihood 
pairlikelihood<-function(Xlist,Ilist,Nlist,u,v){
  
  l<-length(Xlist)
  notvlist<-c()
  vlist<-c()
  for(i in c(1:l)){
    if(length(intersect(Ilist[[i]][-1],v))==0){
      notvlist<-append(notvlist,i)
    }
    else{
      vlist<-append(vlist,i)
    }
  }
  
  
  polys<-list()
  derpolys<-list()
  for(i in c(1:length(notvlist))){
    a<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,u])
    b<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,v])
    c<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,v])
    
    polys[[i]]<-mp(paste0(a," l l-",2*b," l+",c))
    derpolys[[i]]<-deriv(polys[[i]],"l")
  }


  llike<-function(l){
    ll<-0
    for(i in c(1:length(notvlist))){
      ll<-ll+(Nlist[notvlist[i]]/2)*log(as.function(polys[[i]])(l))
    }
    return(ll)
  }


  derllike<-function(l){
    derll<-0
    for(i in c(1:length(notvlist))){
      derll<-derll+(Nlist[notvlist[i]]/2)*(as.function(derpolys[[i]])(l)/as.function((polys[[i]])(l)))
    }
    return(derll)
  }

  SOL<-optim(0,method="Brent",llike,derllike,lower=-1,upper=1)
  
  maxloglik<-SOL$value
  
  if(length(vlist)>0){
    for(i in c(1:length(vlist))){
      
      a<-(Xlist[[vlist[i]]][,u]%*%Xlist[[vlist[i]]][,u])
      b<-(Xlist[[vlist[i]]][,u]%*%Xlist[[vlist[i]]][,v])
      c<-(Xlist[[vlist[i]]][,v]%*%Xlist[[vlist[i]]][,v])
      
      maxloglik<-maxloglik+(Nlist[vlist[i]]/2)*log(as.function(mp(paste0(a," l l-",2*b," l+",c)))(b/a))
    }
  }
  for(i in c(1:l)){
    maxloglik<-maxloglik+(Nlist[i]/2)*log((Xlist[[i]][,u]%*%Xlist[[i]][,u]))
  }
  return(-maxloglik)
}


#Computes the maximal marginal log-likelihood of the model u->v<-w
#as in section 5.2 of the paper
#INPUT:   Xlist= list of samples
#         Ilist= list of intervened nodes
#         Nlist= list of sample sizes
#         u,v,w  indeces of the nodes
#OUTPUT:  Maximal log-likelihood 


tripletlikelihood_coll<-function(Xlist,Ilist,Nlist,u,v,w){
  l<-length(Xlist)
  notvlist<-c()
  vlist<-c()
  for(i in c(1:l)){
    if(length(intersect(Ilist[[i]][-1],v))==0){
      notvlist<-append(notvlist,i)
    }
    else{
      vlist<-append(vlist,i)
    }
  }
  

  polys<-list()
  derpolys_u<-list()
  derpolys_w<-list()
  for(i in c(1:length(notvlist))){
    a<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,u])
    b<-(Xlist[[notvlist[i]]][,w]%*%Xlist[[notvlist[i]]][,w])
    c<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,v])
    d<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,u])
    e<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,w])
    f<-(Xlist[[notvlist[i]]][,w]%*%Xlist[[notvlist[i]]][,u])  
    
    polys[[i]]<-mp(paste0(a," l_u l_u+",b," l_w l_w+",c,"-",2*d," l_u-",2*e," l_w+",2*f," l_u l_w"))
    derpolys_u[[i]]<-deriv(polys[[i]],"l_u")
    derpolys_w[[i]]<-deriv(polys[[i]],"l_w")
  }
  
  llike<-function(l){
    ll<-0
    for(i in c(1:length(notvlist))){
      ll<-ll+(Nlist[notvlist[i]]/2)*log(as.function(polys[[i]])(l))
    }
    return(ll)
  }
  
  
  derllike<-function(l){
    derll<-c(0,0)
    for(i in c(1:length(notvlist))){
      derll<-derll+(Nlist[notvlist[i]]/(2*as.function(polys[[i]])(l)))*c(as.function(derpolys_u[[i]])(l),
                     as.function(derpolys_w[[i]])(l))
    }
    return(derll)
  }
  
  SOL<-optim(c(0,0),method="L-BFGS-B",llike,derllike,lower=c(-1,-1),upper=c(1,1))
  maxloglik<-SOL$value
  
  if(length(vlist)>0){
    for(i in c(1:length(vlist))){
      a<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,u])
      b<-(Xlist[[notvlist[i]]][,w]%*%Xlist[[notvlist[i]]][,w])
      c<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,v])
      d<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,u])
      e<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,w])
      f<-(Xlist[[notvlist[i]]][,w]%*%Xlist[[notvlist[i]]][,u])  
      
      pol<-mp(paste0(a," l_u l_u+",b," l_w l_w+",c,"-",2*d," l_u-",2*e," l_w+",2*f," l_u l_w"))
      lambda_u<-(d*b-f*e/(a*b-f^2))
      lambda_w<-(e*a-f*d/(a*b-f^2))
      
      maxloglik<-maxloglik+(Nlist[vlist[i]]/2)*log(as.function(pol)(c(lambda_u,lambda_w)))
    }
  }
  for(i in c(1:l)){
    maxloglik<-maxloglik+(Nlist[i]/2)*log((Xlist[[i]][,w]%*%Xlist[[i]][,w]))
  }
  return(-maxloglik)
}


#Computes the maximal marginal log-likelihood of the model u->v->w
#as in section 5.2 of the paper
#INPUT:   Xlist=list of samples
#         Ilist=list of intervened nodes
#         Nlist=list of sample sizes
#         u,v,w indeces of the nodes
#OUTPUT:  Maximal log-likelihood 
tripletlikelihood_noncoll<-function(Xlist,Ilist,Nlist,u,v,z){
  
  l<-length(Xlist)
  notvlist<-c()
  vlist<-c()
  for(i in c(1:l)){
    if(length(intersect(Ilist[[i]][-1],v))==0){
      notvlist<-append(notvlist,i)
    }
    else{
      vlist<-append(vlist,i)
    }
  }
  
  
  
  polys<-list()
  derpolys<-list()
  for(i in c(1:length(notvlist))){
    a<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,u])
    b<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,v])
    c<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,v])
    
    polys[[i]]<-mp(paste0(a," l l-",2*b," l+",c))
    derpolys[[i]]<-deriv(polys[[i]],"l")
  }
  
  
  llike<-function(l){
    ll<-0
    for(i in c(1:length(notvlist))){
      ll<-ll+(Nlist[notvlist[i]]/2)*log(as.function(polys[[i]])(l))
    }
    return(ll)
  }
  
  
  derllike<-function(l){
    derll<-0
    for(i in c(1:length(notvlist))){
      derll<-derll+(Nlist[notvlist[i]]/2)*(as.function(derpolys[[i]])(l)/as.function((polys[[i]])(l)))
    }
    return(derll)
  }
  
  SOL<-optim(0,method="Brent",llike,derllike,lower=-1,upper=1)
  
  maxloglik<-SOL$value

  if(length(vlist)>0){
    for(i in c(1:length(vlist))){
      
      a<-(Xlist[[vlist[i]]][,u]%*%Xlist[[vlist[i]]][,u])
      b<-(Xlist[[vlist[i]]][,u]%*%Xlist[[vlist[i]]][,v])
      c<-(Xlist[[vlist[i]]][,v]%*%Xlist[[vlist[i]]][,v])
      
      maxloglik<-maxloglik+(Nlist[vlist[i]]/2)*log(as.function(mp(paste0(a," l l-",2*b," l+",c)))(b/a))
    }
  }
  
  notwlist<-c()
  wlist<-c()
  for(i in c(1:l)){
    if(length(intersect(Ilist[[i]][-1],w))==0){
      notwlist<-append(notwlist,i)
    }
    else{
      wlist<-append(wlist,i)
    }
  }
  
  

  
  polys<-list()
  derpolys<-list()
  for(i in c(1:length(notwlist))){
    a<-(Xlist[[notwlist[i]]][,v]%*%Xlist[[notwlist[i]]][,v])
    b<-(Xlist[[notwlist[i]]][,v]%*%Xlist[[notwlist[i]]][,w])
    c<-(Xlist[[notwlist[i]]][,w]%*%Xlist[[notwlist[i]]][,w])
    
    polys[[i]]<-mp(paste0(a," l l-",2*b," l+",c))
    derpolys[[i]]<-deriv(polys[[i]],"l")
  }
  
  
  llike<-function(l){
    ll<-0
    for(i in c(1:length(notwlist))){
      ll<-ll+(Nlist[notwlist[i]]/2)*log(as.function(polys[[i]])(l))
    }
    return(ll)
  }
  
  
  derllike<-function(l){
    derll<-0
    for(i in c(1:length(notwlist))){
      derll<-derll+(Nlist[notwlist[i]]/2)*(as.function(derpolys[[i]])(l)/as.function((polys[[i]])(l)))
    }
    return(derll)
  }
  
  SOL<-optim(0,method="Brent",llike,derllike,lower=-1,upper=1)
  
  maxloglik<-maxloglik+SOL$value
  
  if(length(wlist)>0){
    for(i in c(1:length(wlist))){
      
      a<-(Xlist[[wlist[i]]][,v]%*%Xlist[[wlist[i]]][,v])
      b<-(Xlist[[wlist[i]]][,v]%*%Xlist[[wlist[i]]][,w])
      c<-(Xlist[[wlist[i]]][,w]%*%Xlist[[wlist[i]]][,w])
      
      maxloglik<-maxloglik+(Nlist[wlist[i]]/2)*log(as.function(mp(paste0(a," l l-",2*b," l+",c)))(b/a))
    }
  }
  return(-maxloglik)
}



#Computes a list of directly I-identifiable edges (using the I-settings)
#definition 2.2 in the paper
#INPUT:   Ilist= list of intervened nodes
#         E= list of unoriented edges
#OUTPUT:  d: vector of the same length of E. d[i]=0 is E[i,] is not directly I-identifiable
#         and it's equal to 1 otherwise.
dir_ident_edges<-function(Ilist,E){
  d<-c()
  for(i in c(1:length(E[,1]))){
    d[i]<-0
    for(j in c(1:length(Ilist))){
      if(length(intersect(Ilist[[j]][-1],E[i,]))==1){
        d[i]<-1
      }
    }
  }
  return(d)
}


#Computes a list of vertices inside the directly identifiable edges
#definition 2.2 in the paper
#INPUT:   d= 0/1 vector compute with dir_ident_edges
#         E= list of unoriented edges
#OUTPUT:  v_list= list of vertices
dir_ident_vert<-function(E,d){
  v_list<-c()
  for(i in c(1:length(d))){
    if(d[i]==1){
      if(length(intersect(v_list,E[i,1]))==0){
        v_list<-append(v_list,E[i,1])
      }
      if(length(intersect(v_list,E[i,2]))==0){
        v_list<-append(v_list,E[i,2])
      }
    }
  }
  return(v_list)
}

#Function that computes an identifiable connected component after the collider search
#INPUT:   E= list of unoriented edges
#         IDvertices= list of possible root vertices
#OUTPUT:  unoriented= list of still unoriented edges
#         tboriented= list of edges in the connected component that has to be oriented


i_component<-function(E,IDvertices){
  
  G<-graph_from_edgelist(E,directed=FALSE)
  v_list<-c(1:max(E))
  V(G)$name<-v_list
  clusters_G<-components(G)
  members<-clusters_G$membership
  c_size<-clusters_G$csize
  
  check<-FALSE
  for(i in c(1:length(IDvertices))){
    ind_int<-which(v_list==IDvertices[i])   #MAYBE USELESS, IDvertices only contains the right ones
    if(check==FALSE){
      if(length(ind_int)>0){
        cluster_ind<-members[[ind_int]]
        if(c_size[cluster_ind]>1){
          id_component<-which(members==cluster_ind)
          check<-TRUE
        }
      }
    }
  }
  
  G_sub<-induced_subgraph(G,id_component)
  E_sub<-get.edgelist(G_sub,names=TRUE)
  
  G_comp<-induced_subgraph(G,v_list[-id_component])
  E<-get.edgelist(G_comp)
  
  return(list(unoriented=E,tboriented=E_sub))
}



#Finds the maximum likelihood root as in 4.2
#INPUT:   Xlist=list of samples
#         Ilist=list of intervened nodes
#         Nlist=list of sample sizes
#         E= list of edges in the connected component
#         IDvertices= list of possible roots in the connected component
aftercollider_test<-function(Xlist,Ilist,Nlist,E,IDvertices){
  G<-graph_from_edgelist(E,directed=FALSE)
  v_list<-V(G)
  likelihood_matrix<-matrix(0,nrow=length(v_list),ncol=length(v_list))
  
  for(i in c(1:length(E[,1]))){
    u<-E[i,1]
    v<-E[i,2]
    likelihood_matrix[u,v]<-aftercollider_condlikelihood(Xlist,Ilist,Nlist,u,v)
    likelihood_matrix[v,u]<-aftercollider_condlikelihood(Xlist,Ilist,Nlist,v,u)
    if(likelihood_matrix[u,u]==0){
      likelihood_matrix[u,u]<-aftercollider_marginallikelihood(Xlist,Ilist,Nlist,u)
    }
    if(likelihood_matrix[v,v]==0){
      likelihood_matrix[v,v]<-aftercollider_marginallikelihood(Xlist,Ilist,Nlist,v)
    }
  }
  
  ll<-c()
  for(i in c(1:length(IDvertices))){
    ll[i]<-root_likelihood(0,likelihood_matrix,IDvertices[i],c(),G)
    ll[i]<-ll[i]+likelihood_matrix[IDvertices[i],IDvertices[i]]
  }
  return(IDvertices[which.max(ll)])
}


#Computes the marginal maximum likelihood of v after the collider search 
#as in 4.2 of the paper
#INPUT:   Xlist=list of samples
#         Ilist=list of intervened nodes
#         Nlist=list of sample sizes
#         v= index of the node
#OUTPUT:  marginal maximum likelihood of v
aftercollider_marginallikelihood<-function(Xlist,Ilist,Nlist,v){
  l<-length(Ilist)
  notvlist<-c()
  vlist<-c()
  
  for(i in c(1:l)){
    if(length(intersect(Ilist[[i]][-1],v))==0){
      notvlist<-append(notvlist,i)
    }
    else{
      vlist<-append(vlist,i)
    }
  }
  
  NN<-sum(Nlist[notvlist[]])
  
  i_variance<-0
  for(i in c(1:length(notvlist))){
    i_variance<-i_variance+Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,v]
  }
  maxloglik<- (-1/NN)*log(i_variance)
  
  if(length(vlist)>0){
    for(i in c(1:length(vlist))){
      maxloglik<-maxloglik-(Nlist[i]/2)*log((Xlist[[i]][,v]%*%Xlist[[i]][,v]))
    }
  }
  return(maxloglik)
}


#Computes the conditional maximum likelihood of v given u after the collider search 
#as in 4.2 of the paper
#INPUT:   Xlist=list of samples
#         Ilist=list of intervened nodes
#         Nlist=list of sample sizes
#         u,v= indices of the nodes
#OUTPUT:  conditional maximum likelihood of v given u
aftercollider_condlikelihood<-function(Xlist,Ilist,Nlist,u,v){
  
  l<-length(Ilist)
  notvlist<-c()
  vlist<-c()
  
  for(i in c(1:l)){
    if(length(intersect(Ilist[[i]][-1],v))==0){
      notvlist<-append(notvlist,i)
    }
    else{
      vlist<-append(vlist,i)
    }
  }
  
  NN<-sum(Nlist[notvlist[]])

  polys<-list()
  i_covariance<-0
  i_length<-0

  for(i in c(1:length(notvlist))){
    a<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,u])
    b<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,v])
    c<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,v])
    i_covariance<-i_covariance+b
    i_length<-i_length+a
    
    polys[[i]]<-mp(paste0(a," l l-",2*b," l+",c))
  }
  
  i_reg_coef<-i_covariance/i_length
  
  i_residual<-0
  
  for(i in c(1:length(notvlist))){
    i_residual<-i_residual+as.function(polys[[i]])(i_reg_coef)
  }
  
  maxloglik<- (NN/2)*log(i_residual)
  
  
  if(length(vlist)>0){
    for(i in c(1:length(vlist))){
      
      a<-(Xlist[[vlist[i]]][,u]%*%Xlist[[vlist[i]]][,u])
      b<-(Xlist[[vlist[i]]][,u]%*%Xlist[[vlist[i]]][,v])
      c<-(Xlist[[vlist[i]]][,v]%*%Xlist[[vlist[i]]][,v])
      
      maxloglik<-maxloglik-(Nlist[vlist[i]]/2)*log(as.function(mp(paste0(a," l l-",2*b," l+",c)))(b/a))
    }
  }
  
  return(maxloglik)

}


#Recursive function that computes the maximum likelihood of the model in which r is the root
#Base case input:   partial_lik=0, 
#                   L=likelihood_matrix as defined in aftercollider_test
#                   r= index of the root
#                   not_neigh= c()
#                   G= subgraph as defined in aftercollider_test
#Output:            maximum likelihood of the model
root_likelihood<-function(partial_lik,L,r,not_neigh,G){
  r_neigh<-neighbors(G,r)
  not_neigh_index<-which(r_neigh==not_neigh)
  if(length(not_neigh)>0){
    if(length(r_neigh)>1){
      r_neigh<-r_neigh[-not_neigh_index]
    }
    else{
      return(partial_lik)
    }
  }
  for(i in c(1:length(r_neigh))){
    partial_lik<-partial_lik+L[r,r_neigh[i]]
    partial_lik<-root_likelihood(partial_lik,L,r_neigh[i],r,G)
  }
  return(partial_lik)
}



