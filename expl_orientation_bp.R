# as fig 1 left
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(35,100,250),
  interventionSize = c(1),
  ndatasets = c(10)+1,
  k = c(1:10),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  conservative = FALSE,
  estimateSkeleton = c("FALSE")
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l1 = l
allResOrientation <- append(allResOrientation, list(df_plot))
ggplot(df_plot, aes(totalSamples, SHD, fill=factor(method))) +
  geom_boxplot() + labs(title=plot_str)


# grid of interventions
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(35),
  interventionSize = c(1,3,5,7,9),
  ndatasets = c(1,3,5,7,10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  conservative = FALSE,
  estimateSkeleton = c("FALSE")
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l_grid_35 = l
allResOrientation <- append(allResOrientation, list(df_plot))


# fig 1 middle
ggplot(df_plot[df_plot$interventionSize == 1,], aes(ndatasets,SHD,fill=factor(method))) + 
  geom_boxplot() + labs(title=paste(plot_str,", interventionSize: 1"))

# fig 1 right
ggplot(df_plot[df_plot$ndatasets == 6,], aes(interventionSize,SHD,fill=factor(method))) + 
  geom_boxplot() + labs(title=paste(plot_str,", ndatasets: 6"))


# grid with more samples
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(70),
  interventionSize = c(1,3,5,7,9),
  ndatasets = c(1,3,5,7,10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  conservative = FALSE,
  estimateSkeleton = c("FALSE")
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l_grid_70 = l
allResOrientation <- append(allResOrientation, list(df_plot))
ggplot(df_plot, aes(NULL, SHD, fill=factor(method))) + geom_boxplot() +
  facet_grid(rows=vars(ndatasets), cols=vars(interventionSize), labeller=labelBoth) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  labs(title=plot_str)
       

# comparison to Yang figure 5
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(100,1000,2500),
  interventionSize = c(1),
  ndatasets = c(20)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  conservative = FALSE,
  estimateSkeleton = c("FALSE")
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l1 = l
ggplot(df_plot, aes(totalSamples,SHD,fill=factor(method))) +
  geom_boxplot() + labs(title=paste("Y fig5: ",plot_str))


# comparison to Hauser figure 10
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(100,1000,10000),
  interventionSize = c(4),
  ndatasets = c(4,12)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  conservative = FALSE,
  estimateSkeleton = c("FALSE")
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l2 = l
ggplot(df_plot, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(rows=vars(ndatasets), labeller = labelBoth) +
  labs(title=paste("HB fig10: ",plot_str))



# check behavior for more datasets with more samples 
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(400),
  interventionSize = c(1,2,4),
  ndatasets = c(1,3,5,7,10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  estimateSkeleton = c("FALSE")
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l_grid_400 = l

df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(70),
  interventionSize = c(1,2,4),
  ndatasets = c(1,3,5,7,10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  conservative = FALSE,
  estimateSkeleton = c("FALSE")
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l_grid_70 = l
Sys.time()


# compare all the grids
plot_interventionGrid = function(l){
  df_plot=l$df; plot_str=l$str
  ggplot(df_plot, aes(NULL, SHD, fill=factor(method))) + geom_boxplot() +
    facet_grid(rows=vars(ndatasets), cols=vars(interventionSize), labeller=labelBoth) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    labs(title=plot_str)
}
plot_interventionGrid(l_grid_400)
plot_interventionGrid(l_grid_1000)



# merge grids
df = rbind(l_grid_70_all$df,l_grid_400_all$df,l_grid_1000_all$df)
# for fig 1
df_plot = df[df$totalSamples==400 & df$interventionSize==2 &
               df$ndatasets %in% c(2,6,11),]
p11 = ggplot(df_plot, aes(ndatasets,SHD,fill=factor(method))) +geom_boxplot()
# for fig 2
df_plot = df[df$ndatasets==11 & df$interventionSize==2,]
p22 = ggplot(df_plot, aes(totalSamples,SHD,fill=factor(method)))+geom_boxplot()

fot1 = (p11 +labs(title=NULL,fill="Method")) +
  (p22 +labs(title=NULL, fill="Method") + theme(legend.position = "None")) +
  plot_layout(nrow = 2)
fot1




# fig1
# want to comment that p1 goes down but p3 up
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(500),
  interventionSize = c(2),
  ndatasets = c(1,10,20)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  conservative = FALSE
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l100 = l
fig100 <- ggplot(df_plot, aes(ndatasets,SHD,fill=factor(method)))+geom_boxplot()+
  labs(title=plot_str)

# fig2
# want to comment that although trend of fig1, still p1 better if large samples
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(100,500,1000),
  interventionSize = c(2),
  ndatasets = c(20)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  conservative = FALSE
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l101 = l
fig101 <- ggplot(df_plot, aes(totalSamples,SHD,fill=factor(method)))+ geom_boxplot()+
  labs(title=plot_str)

# merge figures
f103 = (fig100 + labs(title=NULL)+theme(legend.position = "NONE")) + 
  (fig101 + labs(title=NULL, fill="Procedure")) + plot_layout(nrow=1)
f103






# -------------------------- from here
# interventionSize = 1
# fig1
# want to comment that p1 goes down but p3 up
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(500),
  interventionSize = c(1),
  ndatasets = c(1,10,20)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l104 = l
# fig100 <- ggplot(df_plot, aes(ndatasets,SHD,fill=factor(method)))+geom_boxplot()+
#   labs(title=plot_str)

# fig2
# want to comment that although trend of fig1, still p1 better if large samples
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(100,500,1000),
  interventionSize = c(1),
  ndatasets = c(20)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l105 = l
# fig101 <- ggplot(df_plot, aes(totalSamples,SHD,fill=factor(method)))+ geom_boxplot()+
#   labs(title=plot_str)



# high dimensional (breaks my local machine!)
# fig1
# want to comment that p1 goes down but p3 up
df_params <- expand.grid(
  tsize = c(100),
  totalSamples = c(1000),
  interventionSize = c(1),
  ndatasets = c(25,50,100)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
)
# l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l106 = l
# fig100 <- ggplot(df_plot, aes(ndatasets,SHD,fill=factor(method)))+geom_boxplot()+
#   labs(title=plot_str)

# fig2
# want to comment that although trend of fig1, still p1 better if large samples
df_params <- expand.grid(
  tsize = c(100),
  totalSamples = c(500,1000,3000),
  interventionSize = c(1),
  ndatasets = c(25)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
)
# l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
l107 = l
# fig101 <- ggplot(df_plot, aes(totalSamples,SHD,fill=factor(method)))+ geom_boxplot()+
#   labs(title=plot_str)



# with final version
df_params <- expand.grid(
  tsize = c(100),
  totalSamples = c(500,1000,3000),
  interventionSize = c(1),
  ndatasets = c(25)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05
)
l <- orientationExploration(df_params, scoreFct = pcalg::shd, sFctName = "SHD")



# with final version
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(70),
  interventionSize = c(1),
  ndatasets = c(10)+1,
  k = c(1:3),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05
)
l <- orientationExploration(df_params, sFctName = "SHD")
