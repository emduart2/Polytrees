#Skeleton
ggplot(skeleton_04,aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+xlab(expression(paste("n"[O])))+
#Orientation
ggplot(or[or$method!="random",],aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+
  geom_hline(yintercept=mm, linetype="dashed")+xlab(expression(paste("n"[O])))+
#DAGs
ggplot(sda,aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+xlab(expression(paste("n"[O])))


