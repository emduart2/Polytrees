install.packages("mpoly")
library("mpoly")

#Computes the marginal likelihood of the model u->v
#INPUT: Xlist=list of samples
#       Ilist=list of intervened nodes
#       Nlist=list of sample sizes
#       u,v indeces of the nodes


n_poly<-function(Xlist,Ilist,Nlist,u,v){
  
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
    c<-(Xlist[[notvlist[i]]][,v]%*%Xlist[[notvlist[i]]][,v])
    b<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,v])
    a<-(Xlist[[notvlist[i]]][,u]%*%Xlist[[notvlist[i]]][,u])
    
    polys[[i]]<-mp(paste0(a," l l-",2*b," l+",c))
    derpolys[[i]]<-deriv(polys[[i]],"l")
  }
  
  print(polys)
  print(derpolys)
  llike<-function(l){
    ll<-0
    derll<-0
    for(i in c(1:length(notvlist))){
      ll<-ll-(Nlist[notvlist[i]]/2)*log(as.function(polys[[i]])(l))
    }
    return(ll)
  }
  
  
  derllike<-function(l){
    derll<-0
    for(i in c(1:length(notvlist))){
      derll<-derll-(Nlist[notvlist[i]]/2)*(as.function(derpolys[[i]])(l)/as.function(polys[[i]])(l))
    }
    return(derll)
  }

  plot(c(-500:500)/500,llike(c(-500:500)/500),type="l")
  maxloglik<-optim(0,method="Brent",llike,derllike,lower=-1,upper=1)$value
  
  if(length(vlist)>0){
    for(i in c(1:length(vlist))){
      c<-(Xlist[[vlist[i]]][,v]%*%Xlist[[vlist[i]]][,v])
      b<-(Xlist[[vlist[i]]][,u]%*%Xlist[[vlist[i]]][,v])
      a<-(Xlist[[vlist[i]]][,u]%*%Xlist[[vlist[i]]][,u])
      
      maxloglik<-maxloglik-(Nlist[vlist[i]]/2)*log(as.function(mp(paste0(a," l l-",2*b," l+",c)))(b/a)) 
    }
  }
  for(i in c(1:l)){
    maxloglik<-maxloglik-(Nlist[i]/2)*log((1/Nlist[i])*(Xlist[[i]][,u]%*%Xlist[[i]][,u]))
  }
  return(maxloglik)
}

