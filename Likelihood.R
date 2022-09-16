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
  
  
  NN<-sum(Nlist[notvlist[]])
  MLAMBDA<-0
  
  polys<-list()
  derpolys<-list()
  for(i in c(1:length(notvlist))){
    a<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,u])
    b<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,v])
    c<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,v])
    
    MLAMBDA<-MLAMBDA+(Nlist[notvlist[i]]/NN)*(b/a)
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
  Par<-SOL$par
  
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
  return(list(ML=maxloglik,NVL=notvlist,VL=vlist,PAR=Par,Mlambda=MLAMBDA))
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
  
  
  NN<-sum(Nlist[notvlist[]])
  MLAMBDA<-0
  
  polys<-list()
  derpolys<-list()
  for(i in c(1:length(notvlist))){
    a<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,u])
    b<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,v])
    c<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,v])
    
    MLAMBDA<-MLAMBDA+(Nlist[notvlist[i]]/NN)*(b/a)
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
  Par<-SOL$par
  
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
  
  
  NN<-sum(Nlist[notwlist[]])
  MLAMBDA<-0
  
  polys<-list()
  derpolys<-list()
  for(i in c(1:length(notwlist))){
    a<-(Xlist[[notwlist[i]]][,v]%*%Xlist[[notwlist[i]]][,v])
    b<-(Xlist[[notwlist[i]]][,v]%*%Xlist[[notwlist[i]]][,w])
    c<-(Xlist[[notwlist[i]]][,w]%*%Xlist[[notwlist[i]]][,w])
    
    MLAMBDA<-MLAMBDA+(Nlist[notwlist[i]]/NN)*(b/a)
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
dir_ident<-function(Ilist,E){
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

