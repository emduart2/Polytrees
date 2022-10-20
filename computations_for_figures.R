# This script contains the function calls to generate the data for the 
# figures of the paper. 
# WARNING: The execution of this script takes long. It is recommended to not
#   execute this script without utilizing multiple parallelized cores.
save_dir = "data_local_comp/"


# ---- figure 1a (high-dim skeleton recovery) ----
# left
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(200,500,1000),
  interventionSize = c(1),
  ndatasets = c(20),
  k = c(1:10),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"1a_01.Rds",sep=""))

# middle``
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(500),
  interventionSize = c(1),
  ndatasets = c(2,11,21),
  k = c(1:10),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"1a_02.Rds",sep=""))

# right
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(500),
  interventionSize = c(1,2,4),
  ndatasets = c(20),
  k = c(1:10),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"1a_03.Rds",sep=""))




# ---- figure 1b (low- and high-dim polytrees and DAGs with GIES)
# low-dim polytrees
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(200,500,1000),
  interventionSize = c(1),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","TP","FP","TN","FN"),
  methods_all = list(c("mean","1simp"),c("mean","3simp"),c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,"1b_ld_pt.Rds",sep=""))

# low-dim dags
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(200,500,1000),
  interventionSize = c(1),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 3
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","TP","FP","TN","FN"),
  methods_all = list(c("mean","1simp"),c("mean","3simp"),c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,"1b_ld_dags.Rds",sep=""))


# high-dim polytrees
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(200,500,1000),
  interventionSize = c(2,10), # we later removed interventionSize=2
  ndatasets = c(11),
  k = c(1:10),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","1simp"),c("mean","3simp")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,"1b_hd_pt.Rds",sep=''))
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(200,500,1000),
  interventionSize = c(2,10), # we later removed interventionSize=2
  ndatasets = c(11),
  k = c(1:10),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,"1b_hd_pt_GIES.Rds",sep=''))


# high-dim dags
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(500,1000,2000,3000), # we later removed 3000 samples
  interventionSize = c(10),
  ndatasets = c(21),
  k = c(1:10),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 5
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","1simp"),c("mean","3simp")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,"1b_hd_dags.Rds",sep=''))
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(500,1000,2000,3000), # we later removed 3000 samples
  interventionSize = c(10),
  ndatasets = c(21),
  k = c(1:10),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 5
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,"1b_hd_dags_GIES.Rds",sep=''))