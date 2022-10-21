# load("data_new/unbalanced_withRandBase/unbalanceddata.RData")
load("data_new/unbalanced(2)/unbalanceddata.RData")

#Skeleton
skelplot = skeleton_04
skelplot$sdatasets = as.character(as.numeric(skelplot$sdatasets) / 1000)
p_s<-ggplot(skelplot,
            aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+
  xlab(expression(paste("w"[o])))+
  theme(legend.position = "bottom",legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(face="bold"))+
  guides(fill=guide_legend(title="")) +
  scale_fill_manual(values=colors_skel)
  # labs(title="p=500, n=1000, d=21, k=10")

#Orientation
breaks = c("P.1, refined, BIC",  "P.1, refined, IRC", "P.1, simple, BIC", "P.1, simple, IRC", 
           "P.2, refined, IRC","P.2, refined, BIC", "P.2, simple, IRC",  "P.2, simple, BIC")
orplot = or
orplot$sdatasets = as.character(as.numeric(orplot$sdatasets) / 1000)
p_o<-ggplot(orplot[orplot$method!="random" & orplot$kindOfIntervention == "perfect",],
            aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+
  geom_hline(yintercept=mm, linetype="dotted")+xlab(expression(paste("w"[o])))+
  theme(legend.position = "bottom",legend.key.size = unit(0.4, 'cm'))+
  theme(axis.title.x = element_text(vjust=-1))+
  guides(fill=guide_legend(title=""))+ylab("") +
  scale_fill_manual(values=colors_ort, breaks=breaks)
  # labs(title="p=500, n=1000, d=21, k=10")

fig_unb_skel_ort = p_s + p_o + plot_layout(nrow=1) + 
  plot_annotation(title="p=500, n=1000, d=21, k=10", 
                  theme=theme(plot.title=element_text(face="bold",hjust=0.5))); fig # 1000 by 500



#DAGs
sdaplot = sda
sdaplot$sdatasets = as.character(as.numeric(sdaplot$sdatasets) / 1000)
p_dag<-ggplot(sdaplot,aes(sdatasets,SHD,fill=factor(method)))+geom_boxplot()+
  xlab(expression(paste("w"[o])))+ labs(title = "p=500, n=1000, d=21, k=10",fill="Method") +
  ylab("SHD")+geom_hline(yintercept=m_dag, linetype="dotted")+
  theme(legend.position = "bottom", plot.title = element_text(face="bold",hjust=0.5))+
  scale_fill_manual(values = colors_ort, breaks=c("GIES","P.2, simple, IRC"))
p_dag


# library(gridExtra)
# dev.new()
# grid.arrange(p_s,p_o,p_dag, top = "n=500, p=1000, d=21, k=10",
#              layout_matrix = matrix(c(1,2,3), ncol=3, byrow=TRUE))


