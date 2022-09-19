library(reshape2)
library(ggplot2)
library(patchwork)
# This code assumes, that a data.frame 'df' generated from Skeleton_exploration
# is in the local workspace.

# boxplot grid: x: interventionsize, y: ndatasets, y in each individual boxplot: SHD 
df_plot <- df[df$tsize == 300 & df$totalSamples == 200, ]
ggplot(df_plot, aes(NULL, SHD, fill=factor(method))) + geom_boxplot() +
  facet_grid(rows=vars(ndatasets), cols=vars(interventionSize),labeller=labeller(interventionSize=label_both, ndatasets=label_both))+
  labs(fill="Methods", title="nodes: 300, samples: 200") + ylim(0,0.5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  
# boxplot grid: x: sdatasets, y: interventionsize, y in each individual boxplot: SHD 
df_plot <- df
df_plot$sdatasets <- as.character(lapply(df$sdatasets, FUN=function(x){paste(x,collapse=", ")}))
ggplot(df_plot, aes(NULL, SHD, fill=factor(method))) + geom_boxplot() +
  facet_grid(rows=vars(interventionSize), cols=vars(sdatasets),labeller=labeller(interventionSize=label_both, sdatasets=label_both))+
  labs(fill="Methods", title="nodes: 100, samples: 100") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# create baseline boxplot (Eventually incorrect!)
warning("The following code still needs to be double-checked!")
df_baseline <- expand.grid(
  tsize = c(100,300),
  totalSamples = ceiling(200/c(2,3,6,9,12,15)),
  k = 1:20
)
vals <- foreach(
  p = df_baseline$tsize,
  s = df_baseline$totalSamples,
  k = df_baseline$k,
  .combine = 'c',
  .verbose = TRUE
) %do% {
  # generate graph & samples
  IS <- isetting(p, 1, 0, c(), s)
  G <- graph_from_adjacency_matrix(IS$gTrued)
  ID <- interventionalData(G,IS$L,IS$targetsI)
  
  # calculate SHD
  CL<-chowLiu(ID$Rs[[1]])
  E_e<-get.edgelist(CL)
  estimated_skeleton<-get.adjacency(CL)
  SHD<-sum(abs(estimated_skeleton-IS$gTrues))/(2*(p-1))

  return(SHD)
}
df_baseline$SHD <- vals
df_plot <- df_baseline
df_plot$totalSamples <- factor(df_baseline$totalSamples, levels=sort(unique(df_baseline$totalSamples),decreasing=TRUE))
ggplot(df_plot, aes(NULL, SHD)) + geom_boxplot() +
  facet_grid(cols=vars(tsize), rows=vars(totalSamples),labeller=labeller(tsize=label_both, totalSamples=label_both))+
  labs(title="Baseline with no interventions") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# boxplot as figure 6 in Yang
df_plot <- df
df_plot$ndatasets <- factor(df_plot$ndatasets)
ggplot(df_plot, aes(ndatasets, SHD, fill=factor(method))) + geom_boxplot() +
  labs(title="one-variable perfect interventions, in total 100 nodes") + ylim(0,0.475)
