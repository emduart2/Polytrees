# this file contains code to recreate the figures of the paper

#---- figure 1-----
p = 20
nSamples = 200
# left
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(70,100,200),
  interventionSize = c(1),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l01 <- skeletonExploration(df_params)

# middle
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nSamples),
  interventionSize = c(1),
  ndatasets = c(2,11,21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l02 <- skeletonExploration(df_params); 


# right
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nSamples),
  interventionSize = c(1,2,4),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l03 <- skeletonExploration(df_params);


# plot figures
f01 = ggplot(l01$df, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot()
f02 = ggplot(l02$df, aes(ndatasets,SHD,fill=factor(method))) + geom_boxplot()
f03 = ggplot(l03$df, aes(interventionSize,SHD,fill=factor(method))) + geom_boxplot()

# figure 1: merge figure together
(f04 = (f01 +labs(title=NULL) +theme(legend.position = "none")) +
    (f02 +labs(title=NULL) +theme(legend.position = "none")) +
    (f03 +labs(title=NULL, fill="aggr fct")) +
    plot_layout(nrow = 1) + plot_annotation(title=l003$str))