library(patchwork)
source("polyFunctions.R")
source("Likelihood.R")
source("expl_orientation.R")
source("expl_skeleton.R")
source("expl_both.R")


# ot: but with interventionSize = 1
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(500),
  interventionSize = c(1),
  ndatasets = c(1,10,20)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  estimateSkeleton = c("FALSE")
)
l104 <- orientationExploration(df_params); 
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(100,500,1000),
  interventionSize = c(1),
  ndatasets = c(20)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  estimateSkeleton = c("FALSE")
)
l105 <- orientationExploration(df_params); 



# both: with increasing datasets
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(200),
  interventionSize = c(1),
  ndatasets = c(1,5,10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l303 <- bothExploration(df_params,methods)


# both: grid
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(400),
  interventionSize = c(1,2,4),
  ndatasets = c(1,5,10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l304 <- bothExploration(df_params,methods)

  
  
# both: increasing number nodes
df_params <- expand.grid(
  tsize = c(10,40),
  totalSamples = c(40),
  interventionSize = c(1),
  ndatasets = c(10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"))
l305 <- bothExploration(df_params,methods)
df_params <- expand.grid(
  tsize = c(80),
  totalSamples = c(40),
  interventionSize = c(1),
  ndatasets = c(10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"))
l306 <- bothExploration(df_params,methods)


# figures
done = Sys.time()
detach("package:mpoly", unload = TRUE)
(fig101 = ggplot(l104$df, aes(ndatasets,SHD,fill=factor(method)))+ geom_boxplot()+
  ggplot(l105$df, aes(totalSamples,SHD,fill=factor(method)))+ geom_boxplot()+
   plot_layout(nrow=1,guides="collect") + plot_annotation(title=l104$str))

(fig303 = 
  ggplot(l303$df, aes(ndatasets, shd_all,fill=factor(method))) +geom_boxplot() + 
  ggplot(l303$df, aes(ndatasets, time_all, fill = factor(method)))+geom_boxplot() +
  plot_layout(nrow=1, guides="collect") + plot_annotation(title=l303$str) )

(fig304 = ggplot(l304$df, aes(NULL,shd_skel,fill=factor(method))) + geom_boxplot()+
  facet_grid(rows=ggplot2::vars(ndatasets),cols=ggplot2::vars(interventionSize),labeller = labelBoth)+
  labs(title=l304$str) +  theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank(),
                                axis.ticks.x=element_blank()))
(fig305 = 
  ggplot(l305$df, aes(tsize, shd_all,fill=factor(method))) +
  geom_boxplot() + theme(legend.position=NULL)+
  ggplot(l305$df, aes(tsize, time_all, fill = factor(method)))+
  geom_boxplot() + plot_layout(nrow=1)+plot_annotation(title=l305$str))




# from here ---------
# figure 1 but up to 20% interventions
# figure 1 left
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(70,200,500),
  interventionSize = c(1),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params); df_plot=l$df; plot_str=l$str
l001 = l

# figure 1 middle
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(70),
  interventionSize = c(1),
  ndatasets = c(2,11,21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params); df_plot=l$df; plot_str=l$str
l002 = l

# figure 1 right
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(70),
  interventionSize = c(1,2,4,10),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params); df_plot=l$df; plot_str=l$str
l003 = l

# plot
f001 = ggplot(l001$df, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot()
f002 = ggplot(l002$df, aes(ndatasets,SHD,fill=factor(method))) + geom_boxplot()
f003 = ggplot(l003$df[l003$df$interventionSize %in% c(1,2,4),], aes(interventionSize,SHD,fill=factor(method))) + geom_boxplot()

# figure 1: merge figure together
(f004 = (f001 +labs(title=NULL) +theme(legend.position = "none")) +
    (f002 +labs(title=NULL) +theme(legend.position = "none")) +
    (f003 +labs(title=NULL, fill="aggr fct")) +
    plot_layout(nrow = 1))
save(l001,l002,l003,f001,f002,f003,f004, file="fig1_data")



# plot of comparison to GIES but with non-perfect interventions
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(200),
  interventionSize = c(1),
  ndatasets = c(10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect","imperfect","inhibitory"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l315 <- bothExploration(df_params,methods);
(f315 = ggplot(l315$df, aes(kindOfIntervention, shd_all, fill=factor(method)))+
                geom_boxplot())



# check whether simple procedures are faster
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(50,200,500),
  interventionSize = c(1),
  ndatasets = c(10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("mean","1"),c("mean","1simp"),c("mean","2"),c("mean","2simp"),c("mean","3"),c("mean","3simp"))
l316 <- bothExploration(df_params,methods);
(f316 = ggplot(l316$df, aes(totalSamples, time_all, fill=factor(method)))+ geom_boxplot()+
   labs(title=l316$str))

# add random interventions to f315
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(200),
  interventionSize = c(1),
  ndatasets = c(10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("random"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l317 = bothExploration(df_params,methods)
(f317 = ggplot(rbind(l315$df,l317$df), aes(kindOfIntervention, shd_all, fill=factor(method)))+
    geom_boxplot() + labs(title=l317$str))


# -------- still need to compute
# fig 1 with SID





#---- large Scale final ----



