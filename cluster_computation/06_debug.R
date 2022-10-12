# paramters
kmax = 2
save_dir = ''  # with final '/' in the end
experiment_id = "06"

# -- parallellisation below here-- 


#---- figure 2 low-dim with correct runtime and diff intervSize -----
# - compare runtime of simple vs refined
# - compare how TEST behaves for different interventionsizes
# left
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
  methods_all = list(c("gtruth","1"),c("gtruth","1simp"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
  pw_methods_all = c("BIC","TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-01.Rds",sep=''))


# right
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
  methods_all = list(c("gtruth","1"),c("gtruth","1simp"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
  pw_methods_all = c("BIC","TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-02.Rds",sep=''))



#----- figure 1 high-dim -------
# left
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11,12),
  interventionSize = c(1),
  ndatasets = c(3),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,experiment_id,"-03.Rds",sep=''))

# middle
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11),
  interventionSize = c(1),
  ndatasets = c(2),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,experiment_id,"-04.Rds",sep=''))

# right
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11),
  interventionSize = c(1,2),
  ndatasets = c(3),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,experiment_id,"-05.Rds",sep=''))


# ------- figure 2 high-dim ----------
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11),
  interventionSize = c(1),
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
  methods_all = list(c("gtruth","1"),c("gtruth","1simp"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
  pw_methods_all = c("BIC","TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-06.Rds",sep=''))


# right
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11),
  interventionSize = c(1),
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
  methods_all = list(c("gtruth","1"),c("gtruth","1simp"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
  pw_methods_all = c("BIC","TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-07.Rds",sep=''))