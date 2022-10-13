save_dir = "local_computations/"

# overnight calculations
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(500),
  interventionSize = c(2),
  ndatasets = c(2,11,21),
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
  methods_all = list(c("mean","1simp"),c("mean","3simp"),c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l1, file = paste(save_dir,"01-01.Rds",sep=""))


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
  methods_all = list(c("mean","1simp"),c("mean","3simp"),c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l2, file = paste(save_dir,"01-02.Rds",sep=""))

Sys.time()



# overday calculations
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(500),
  interventionSize = c(2),
  ndatasets = c(2,11,21),
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
  methods_all = list(c("mean","2"),c("mean","2simp")),
  pw_methods_all = c("TEST"))
saveRDS(l1, file = paste(save_dir,"01-03.Rds",sep=""))

df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(500),
  interventionSize = c(2),
  ndatasets = c(2,11,21),
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
  methods_all = list(c("mean","1"),c("mean","1simp"),c("mean","2"),c("mean","2simp"),c("mean","3"),c("mean","3simp")),
  pw_methods_all = c("BIC"))
saveRDS(l1, file = paste(save_dir,"01-04.Rds",sep=""))




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
l1 <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","2"),c("mean","2simp")),
  pw_methods_all = c("TEST"))
saveRDS(l1, file = paste(save_dir,"01-05.Rds",sep=""))

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
l1 <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","1"),c("mean","1simp"),c("mean","2"),c("mean","2simp"),c("mean","3"),c("mean","3simp")),
  pw_methods_all = c("BIC"))
saveRDS(l1, file = paste(save_dir,"01-06.Rds",sep=""))

Sys.time()
