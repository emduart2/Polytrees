# ----------------------
# THIS FILE IS DEPRECATED
# ----------------------

# setup
source("polyFunctions.R")
source("Likelihood.R")
source("expl_skeleton.R")
source("expl_orientation.R")
source("expl_both.R")
library(patchwork)
labelBoth = labeller(.default=label_both)


#------------figure 1 -----
p = 100
nSamples = 100
# left
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(70,100,200),
  interventionSize = c(1),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params)
l001 = l

# middle
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nSamples),
  interventionSize = c(1),
  ndatasets = c(2,11,21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params); 
l002 = l

# right
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nSamples),
  interventionSize = c(1,2,4),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params);
l003 = l

# plot figures
f001 = ggplot(l001$df, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot()
f002 = ggplot(l002$df, aes(ndatasets,SHD,fill=factor(method))) + geom_boxplot()
f003 = ggplot(l003$df, aes(interventionSize,SHD,fill=factor(method))) + geom_boxplot()

# figure 1: merge figure together
(f004 = (f001 +labs(title=NULL) +theme(legend.position = "none")) +
  (f002 +labs(title=NULL) +theme(legend.position = "none")) +
  (f003 +labs(title=NULL, fill="aggr fct")) +
  plot_layout(nrow = 1) + plot_annotation(title=l003$str))





# plot reference figure for figure 1
df_params <- expand.grid(
  tsize = c(100),
  totalSamples = c(100),
  interventionSize = c(1),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params);
ggplot(l$df, aes(totalSamples,SHD,fill=factor(method)))+geom_boxplot() +
  labs(title=l$str)









# simulate data
df_params <- expand.grid(
  tsize = c(100),
  totalSamples = c(300),
  interventionSize = c(1,2,4),
  ndatasets = c(20,40,60,80,100),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  conservative = FALSE
)
res <- skeletonExploration(df_params,allResults)
df = res$df; allResults = res$allResults
res_plot <- prepare_df_plot(df)
df_plot = res_plot$df; plot_str = res_plot$str

# consistency
ggplot(df_plot,aes(totalSamples,SHD,fill=factor(method)))+geom_boxplot()+
  labs(title=plot_str)

# slide 1
ggplot(df_plot,aes(ndatasets,SHD,fill=factor(method))) + geom_boxplot() +
  labs(title=plot_str)

# slide 2
ggplot(df_plot, aes(interventionSize,SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=vars(kindOfIntervention),rows=vars(ndatasets),
             labeller=labeller(kindOfIntervention=label_both,interventionSize=label_both,ndatasets=label_both))+
  labs(title=plot_str)

# slide 3 left
ggplot(df_plot, aes(NULL, SHD, fill=factor(method))) + geom_boxplot() +
  facet_grid(rows=vars(ndatasets), cols=vars(interventionSize), labeller=labellerAll)+
  labs(title=plot_str)

# slide 3 right
ggplot(df_plot, aes(NULL, SHD, fill=factor(method))) + geom_boxplot() +
  facet_grid(rows=vars(ndatasets), cols=vars(interventionSize), labeller=labellerAll)+
  labs(title=plot_str)

# slide 4
ggplot(df_plot, aes(ndatasets,SHD, fill=factor(method))) + geom_boxplot()+
  facet_grid(cols=vars(interventionSize), labeller=labellerAll)+labs(title=plot_str)










# df <- rbind(allResults[[11]]$df_vals, df)
df_plot = df[df$kindOfIntervention == "perfect",]
df_plot$ndatasets <- factor(df_plot$ndatasets, levels=sort(unique(df_plot$ndatasets),decreasing = FALSE))
df_plot$interventionSize <- factor(df_plot$interventionSize, levels=sort(unique(df_plot$interventionSize),decreasing = FALSE))
# df_plot = df_plot[df_plot$ndatasets %in% factor(c(20,40,60,80,100)),]
ggplot(df_plot, aes(ndatasets,SHD,fill=factor(method)))+ geom_boxplot() +
  facet_grid(cols=vars(interventionSize),labeller = labelBoth) +
  labs(title="nodes: 100, totalSamples: 1000, kindOfIntervention: perfect, conservative=TRUE")

#
df_plot = df[df$kindOfIntervention == "perfect",]
df_plot$ndatasets <- factor(df_plot$ndatasets, levels=sort(unique(df_plot$ndatasets),decreasing = FALSE))
df_plot$interventionSize <- factor(df_plot$interventionSize, levels=sort(unique(df_plot$interventionSize),decreasing = FALSE))
ggplot(df_plot, aes(NULL,SHD,fill=factor(method)))+ geom_boxplot() +
  facet_grid(cols=vars(interventionSize),rows=vars(ndatasets),labeller = labelBoth) +
  labs(title="nodes: 40, totalSamples: 300, kindOfIntervention: perfect, conservative=FALSE")


# mini benchmark for mean vs median
w = dnorm(-10:9)
w = w/sum(w)
benchmark(
  weighted.mean(sample(20,20,replace=TRUE),w=w),
  weighted.median(sample(20,20,replace=TRUE),w=w),
  replications=1000
)

# as in Hauser
df_plot = df
df_plot$ndatasets <- factor(df_plot$ndatasets, levels=sort(unique(df_plot$ndatasets),decreasing = FALSE))
df_plot$interventionSize <- factor(df_plot$interventionSize, levels=sort(unique(df_plot$interventionSize),decreasing = FALSE))
ggplot(df_plot, aes(ndatasets,SHD,fill=factor(method)))+geom_boxplot()+
  facet_grid(cols=vars(kindOfIntervention),labeller=labellerAll)+
  labs(title="nodes: 40, interventionsize: 1, samples: 300, intervention: perfect")

















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

# boxplot: samplesize
df_plot <- df
df_plot$totalSamples <- factor(df_plot$totalSamples, levels=sort(unique(df_plot$totalSamples),decreasing = FALSE))
ggplot(df_plot, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=vars(kindOfIntervention),labeller=labeller(kindOfIntervention=label_both))+
  labs(title="nodes: 10, all one-node interventions")

# boxplot: samplesize with grid
df_plot <- df
df_plot$totalSamples <- factor(df_plot$totalSamples, levels=sort(unique(df_plot$totalSamples),decreasing = FALSE))
ggplot(df_plot, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=vars(kindOfIntervention),rows=vars(interventionSize),
             labeller=labeller(kindOfIntervention=label_both,interventionSize=label_both))+
  labs(title="nodes: 100, number of interventions: 10")



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


#
df_plot <- df
df_plot$totalSamples <- factor(df_plot$totalSamples, levels=sort(unique(df_plot$totalSamples),decreasing = FALSE))
df_plot$interventionSize <- factor(df_plot$interventionSize, levels=sort(unique(df_plot$interventionSize),decreasing = FALSE))
ggplot(df_plot, aes(interventionSize,SHD,fill=factor(method))) + geom_boxplot() +
  facet_grid(cols=vars(kindOfIntervention),
             labeller=labeller(kindOfIntervention=label_both,interventionSize=label_both))+
  labs(title="nodes: 100, samples: 300, number of interventions: 10")
