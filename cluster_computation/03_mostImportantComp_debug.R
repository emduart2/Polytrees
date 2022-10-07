# paramters
kmax = 2
save_dir = 'cluster_computation/'  # with final '/' in the end
experiment_id = "03"

# -- parallellisation below here--



#---- fig 2 but high-dim ----
# - pick one of interventionsize (whatever look sbetter)
# - analyse the trends of varying totalSamples & ndatasets
# varying ndatases
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11),
  interventionSize = c(1,2),
  ndatasets = c(2,3),
  k = c(1:kmax),
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
saveRDS(l, file = paste(save_dir,experiment_id,"-01.Rds",sep=''))

# varying ndatases with GIES
# - check runtime of GIES 
# - if it terminates, we want to say that we are faster
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11),
  interventionSize = c(1,2),
  ndatasets = c(2,3),
  k = c(1:kmax),
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
saveRDS(l, file = paste(save_dir,experiment_id,"-01_GIES.Rds",sep=''))


# varying samplesize
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11,12),
  interventionSize = c(1,2),
  ndatasets = c(3),
  k = c(1:kmax),
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
saveRDS(l, file = paste(save_dir,experiment_id,"-03.Rds",sep=''))



#---- fig 2 but with all intervention types ----
# - we want to see if there is any difference with different intervention types
# left
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11),
  interventionSize = c(1),
  ndatasets = c(2,3),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect","imperfect","inhibitory","random"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","1"),c("mean","1simp"),c("mean","2"),c("mean","2simp"),c("mean","3"),c("mean","3simp")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-04.Rds",sep=''))

# right
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11,12),
  interventionSize = c(1),
  ndatasets = c(3),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect","imperfect","inhibitory","random"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","1"),c("mean","1simp"),c("mean","2"),c("mean","2simp"),c("mean","3"),c("mean","3simp")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-05.Rds",sep=''))




# #---- potential fig 3: high-dim DAG setting for varying nsamples ----
# - pick one interventionsize (whatever looks better)
# - compare runtime of algorithms, ours should be faster than GIES
# nbh = 1
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11),
  interventionSize = c(1,2),
  ndatasets = c(3),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 1
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","1simp"),c("mean","3simp")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-06.Rds",sep=''))

# nbh = 1, GIES
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11),
  interventionSize = c(1,2),
  ndatasets = c(3),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 1
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-06_GIES.Rds",sep=''))


# nbh = 5
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(11),
  interventionSize = c(1,2),
  ndatasets = c(3),
  k = c(1:kmax),
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
saveRDS(l, file = paste(save_dir,experiment_id,"-07.Rds",sep=''))

# nbh = 5, GIES
df_params <- expand.grid(
  tsize = c(10),
  totalSamples = c(11),
  interventionSize = c(1,2),
  ndatasets = c(3),
  k = c(1:kmax),
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
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-07_GIES.Rds",sep=''))



# nbh = 10
df_params <- expand.grid(
  tsize = c(15),
  totalSamples = c(11),
  interventionSize = c(1,2),
  ndatasets = c(3),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 10
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","1simp"),c("mean","3simp")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-08.Rds",sep=''))

# nbh = 10, GIES
df_params <- expand.grid(
  tsize = c(15),
  totalSamples = c(11),
  interventionSize = c(1,2),
  ndatasets = c(3),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 10
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-08_GIES.Rds",sep=''))
