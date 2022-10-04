# The idea is that all parameters that determine the setting of the graph are put
#   into df_params. Which functions we use and how we evaluate them gets passed
#   to the exploration function directly.
# The result is a list with a field data.frame that can be nicely plotted with 
#   ggplot.
library(ggarrange)

# explore skeleton
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(70,100,200),
  interventionSize = c(1),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l1 <- explore_skeleton(df_params)
ggplot(l1$df, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot() +
  labs(title=l1$str)


# explore only orientation (i.e. in methods_all set c("gtruth",...) to work with 
#   ground Truth skeleton)
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(10,20),
  interventionSize = c(1),
  ndatasets = c(3),
  k = c(1:3),
  sdatasets = list(c()),
  kindOfIntervention = c("perfekt"),
  ensureDiff = TRUE,
  alpha = c(0.05),
  use_dags = FALSE,
  dag_nbh = 0
)
l2 = explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID), 
  sFctNames_all = c("SHD","SID"),
  methods_all = list(c("mean","1"),c("mean","3"),c("gtruth","3")),
  pw_methods_all = c("BIC","TEST"))
ggplot(l2$df, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot()+
  ggplot(l2$df, aes(totalSamples,SID,fill=factor(method))) + geom_boxplot()+
  ggplot(l2$df, aes(totalSamples,time,fill=factor(method))) + geom_boxplot()+
  plot_layout(nrow=1) + plot_annotation(title=l2$str)


# explore both skeleton and exploration and compare with GIES
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(50),
  interventionSize = c(1),
  ndatasets = c(3),
  k = c(1:3),
  sdatasets = list(c()),
  kindOfIntervention = c("perfekt"),
  ensureDiff = TRUE,
  alpha = c(0.05),
  use_dags = FALSE,
  dag_nbh = 0
)
l3 = explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID), 
  sFctNames_all = c("SHD","SID"),
  methods_all = list(c("GIES","GIES"),c("mean","1"),c("mean","3")),
  pw_methods_all = c("BIC"))
ggplot(l3$df, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot()+
  ggplot(l3$df, aes(totalSamples,SID,fill=factor(method))) + geom_boxplot()+
  ggplot(l3$df, aes(totalSamples,time,fill=factor(method))) + geom_boxplot()+
  plot_layout(nrow=1) + plot_annotation(title=l3$str)
