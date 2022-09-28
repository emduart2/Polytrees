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
       


# print time
# Sys.time()


# easy 
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(500),
  interventionSize = c(0),
  ndatasets = c(0)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  conservative = FALSE,
  estimateSkeleton = c("FALSE")
)
l <- orientationExploration(df_params,allResOrientation); df_plot=l$df; plot_str=l$str
ggplot(df_plot, aes(ndatasets,SHD,fill=factor(method))) + 
  geom_boxplot() + labs(title=paste(plot_str,", interventionSize: 1"))
