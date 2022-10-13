# this file contains code to recreate the figures of the paper

#---- figure 1-----
p = 20
nSamples = 70
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
l01 <- explore_skeleton(df_params)

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
l02 <- explore_skeleton(df_params); 


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
l03 <- explore_skeleton(df_params);


# plot figures
f01 = ggplot(l01$df, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot()
f02 = ggplot(l02$df, aes(ndatasets,SHD,fill=factor(method))) + geom_boxplot()
f03 = ggplot(l03$df, aes(interventionSize,SHD,fill=factor(method))) + geom_boxplot()

# figure 1: merge figure together
(fig1 = (f01 +labs(title=NULL) +theme(legend.position = "none")) +
    (f02 +labs(title=NULL) +theme(legend.position = "none")) +
    (f03 +labs(title=NULL, fill="aggr fct")) +
    plot_layout(nrow = 1) + plot_annotation(title=l03$str))





#---- figure 2 -----
kmax = 5
# left
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(70,200),
  interventionSize = c(1),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = FALSE,
  dag_nbh = 0
)
l1 <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","1"),c("mean","1simp"),c("mean","2"),c("mean","2simp"),c("mean","3"),c("mean","3simp")),
  pw_methods_all = c("BIC"))

# right
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(100),
  interventionSize = c(1),
  ndatasets = c(1,21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = FALSE,
  dag_nbh = 0
)
l2 <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","1"),c("mean","1simp"),c("mean","2"),c("mean","2simp"),c("mean","3"),c("mean","3simp")),
  pw_methods_all = c("BIC"))

ggplot(df1, aes(totalSamples, SHD, fill=factor(method))) + geom_boxplot() +
  ggplot(df2, aes(ndatasets, SHD, fill=factor(method))) + geom_boxplot() + 
  plot_layout(nrow = 1, guides="collect") + plot_annotation(title=l2$str)

