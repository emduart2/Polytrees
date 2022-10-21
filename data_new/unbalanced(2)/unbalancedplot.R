dev.new()
#Skeleton
p_s<-ggplot(skeleton_04,aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+xlab(expression(paste("n"[o])))+
  guides(fill=guide_legend(title=""))+facet_grid(cols = vars(c(purpose)))+
    theme(text = element_text(size = 20),legend.position = "bottom",
          legend.key.size = unit(1, 'cm'),legend.text = element_text(size=10))
#Orientation
p_o<-ggplot(or[or$method!="random",],aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+
  geom_hline(yintercept=mm, linetype="dashed")+xlab(expression(paste("n"[0])))+
      guides(fill=guide_legend(title=""))+ylab("")+facet_grid(cols = vars(c(purpose)))+
        theme(text = element_text(size = 20),legend.position = "bottom",
              legend.key.size = unit(0.7, 'cm'),legend.text = element_text(size=8))
dev.new()
p_s+p_o+plot_annotation(title = "n=500, p=1000, d=21, k=10") & 
  theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))


#DAGs
p_dag<-ggplot(sda,aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+xlab(expression(paste("n"[O])))+
  geom_hline(yintercept=m_dag, linetype="dashed")+ggtitle("n=500, p=1000, d=21, k=10, e=10")+
  guides(fill=guide_legend(title=""))

p_dag
