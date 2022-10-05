# ----------------------
# THIS FILE IS DEPRECATED
# ----------------------



library(patchwork)


# get a feeling for the algorithm
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(200),
  interventionSize = c(1),
  ndatasets = c(10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("mean","1"),c("mean","3"))
l301 <- bothExploration(df_params,methods); df_plot=l$df; plot_str=l$str

f300 = ggplot(l200$df, aes(totalSamples, shd_all,fill=factor(method))) +
  geom_boxplot() + labs(title=l$str)
f301 = ggplot(l200$df, aes(totalSamples, time_all, fill = factor(method)))+
  geom_boxplot()
f302 = ggplot(l200$df, aes(totalSamples, time_ort, fill = factor(method)))+
  geom_boxplot()
f303 = (f300+theme(legend.position="NONE")) + (f301+theme(legend.position="NONE")) +
  f302 + plot_layout(nrow=1)


# compare runtimes
df_params <- expand.grid(
  tsize = c(60),
  totalSamples = c(350),
  interventionSize = c(1),
  ndatasets = c(10)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l302 <- bothExploration(df_params,methods); df_plot=l$df; plot_str=l$str
f302 <- ggplot(l302$df, aes(totalSamples, shd_all,fill=factor(method))) +
  geom_boxplot() + labs(title=l302$str) + 
  ggplot(l302$df, aes(totalSamples, time_all, fill = factor(method)))+
  geom_boxplot() + plot_layout(nrow=1,guides='collect')
f302



# test against DAGs
df_params <- expand.grid(
  tsize = c(11),
  totalSamples = c(5000),
  interventionSize = c(1),
  ndatasets = c(5)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l310 <- bothExploration(df_params,methods); df_plot=l$df; plot_str=l$str
(f310 <- ggplot(l310$df, aes(totalSamples, shd_all,fill=factor(method))) +
    geom_boxplot() + 
    ggplot(l310$df, aes(totalSamples, time_all, fill = factor(method)))+
    geom_boxplot() + plot_layout(nrow=1,guides='collect')+
    plot_annotation(title=paste("nbh: 2, ",l310$str)))


df_params <- expand.grid(
  tsize = c(11),
  totalSamples = c(5000),
  interventionSize = c(1),
  ndatasets = c(5)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l311 <- bothExploration(df_params,methods); df_plot=l$df; plot_str=l$str
(f311 <- ggplot(l311$df, aes(totalSamples, shd_all,fill=factor(method))) +
    geom_boxplot() + 
    ggplot(l311$df, aes(totalSamples, time_all, fill = factor(method)))+
    geom_boxplot() + plot_layout(nrow=1,guides='collect')+
    plot_annotation(title=paste("nbh: 5, ",l311$str)))



df_params <- expand.grid(
  tsize = c(11),
  totalSamples = c(5000),
  interventionSize = c(1),
  ndatasets = c(5)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l312 <- bothExploration(df_params,methods); df_plot=l$df; plot_str=l$str
(f312 <- ggplot(l312$df, aes(totalSamples, shd_all,fill=factor(method))) +
    geom_boxplot() + 
    ggplot(l312$df, aes(totalSamples, time_all, fill = factor(method)))+
    geom_boxplot() + plot_layout(nrow=1,guides='collect')+
    plot_annotation(title=paste("nbh: 3.5, ",l312$str)))



df_params <- expand.grid(
  tsize = c(11),
  totalSamples = c(5000),
  interventionSize = c(1),
  ndatasets = c(5)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l313 <- bothExploration(df_params,methods); df_plot=l$df; plot_str=l$str
(l313 <- ggplot(l313$df, aes(totalSamples, shd_all,fill=factor(method))) +
    geom_boxplot() + 
    ggplot(l313$df, aes(totalSamples, time_all, fill = factor(method)))+
    geom_boxplot() + plot_layout(nrow=1,guides='collect')+
    plot_annotation(title=paste(l313$str)))



# try high-dimensional setting
df_params <- expand.grid(
  tsize = c(40),
  totalSamples = c(20,40),
  interventionSize = c(1),
  ndatasets = c(5)+1,
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l314 <- bothExploration(df_params,methods); 
(f314 <- ggplot(l314$df, aes(totalSamples, shd_all,fill=factor(method))) +
    geom_boxplot() + 
    ggplot(l314$df, aes(totalSamples, time_all, fill = factor(method)))+
    geom_boxplot() + plot_layout(nrow=1,guides='collect')+
    plot_annotation(title=paste("nbh: 3.5, ",l314$str)))



# check agains dags
# check high-dim dags with different nbh
df_params <- expand.grid(
  tsize = c(50),
  totalSamples = c(20,50,100),
  interventionSize = c(2),
  ndatasets = c(5),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = c(0.05),
  use_dags = TRUE,
  dag_nbh = 4
)
l10 = explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives), 
  sFctNames_all = c("SHD","SID","TP","FP"),
  methods_all = list(c("GIES","GIES"),c("mean","1"),c("mean","3")),
  pw_methods_all = c("BIC"))
t1 = Sys.time()

df_params <- expand.grid(
  tsize = c(50),
  totalSamples = c(20,50,100),
  interventionSize = c(2),
  ndatasets = c(5),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = c(0.05),
  use_dags = TRUE,
  dag_nbh = 10
)
l11 = explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID,true_positives, false_positives), 
  sFctNames_all = c("SHD","SID","TP","FP"),
  methods_all = list(c("GIES","GIES"),c("mean","1"),c("mean","3")),
  pw_methods_all = c("BIC"))
t2 = Sys.time()


l=l10
ggplot(l$df, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot() +
  ggplot(l$df, aes(totalSamples,SID,fill=factor(method))) + geom_boxplot() +
  ggplot(l$df, aes(totalSamples,TP,fill=factor(method))) + geom_boxplot() +
  ggplot(l$df, aes(totalSamples,FP,fill=factor(method))) + geom_boxplot() +
  ggplot(l$df, aes(totalSamples,time,fill=factor(method))) + geom_boxplot()+
  plot_layout(ncol=3) + plot_annotation(title=l$str)


l = l11
ggplot(l$df, aes(totalSamples,SHD,fill=factor(method))) + geom_boxplot() +
  ggplot(l$df, aes(totalSamples,SID,fill=factor(method))) + geom_boxplot() +
  ggplot(l$df, aes(totalSamples,TP,fill=factor(method))) + geom_boxplot() +
  ggplot(l$df, aes(totalSamples,FP,fill=factor(method))) + geom_boxplot() +
  ggplot(l$df, aes(totalSamples,time,fill=factor(method))) + geom_boxplot()+
  plot_layout(ncol=3) + plot_annotation(title=l$str)





#-------- from here on with new explore function ----------
# figure 2 with comparison
# left
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(500),
  interventionSize = c(2),
  ndatasets = c(2,11,21),
  k = c(1:1),
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
  methods_all = list(c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
  pw_methods_all = c("BIC","TEST"))


# right
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(100,500,1000),
  interventionSize = c(2),
  ndatasets = c(21),
  k = c(1:20),
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
  pw_methods_all = c("BIC","TEST"))




# plot figures
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(500),
  interventionSize = c(2),
  ndatasets = c(2,21),
  k = c(1:20),
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
  methods_all = list(c("gtruth","PearlAll"),c("gtruth","PearlObs"),c("gtruth","1")),
  pw_methods_all = c("BIC"))



ggplot(l1$df, aes(ndatasets,SHD,fill=factor(method))) + geom_boxplot()

