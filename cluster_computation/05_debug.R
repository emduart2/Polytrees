# paramters
kmax = 2
save_dir = ''  # with final '/' in the end
experiment_id = "05"

# -- parallellisation below here-- 



#---- fig 2 but high-dim ----
# compute values for GIES which we forgot for one setting in computation 03
# - analyse the trends of varying totalSamples & ndatasets
# varying samplesize
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11,12),
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
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-03_GIES.Rds",sep=''))



#---- optional: fig2 but with BIC ----
# - compare against TEST from computation03 in terms of runtime, SHD, SID
# varying ndatases
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
  methods_all = list(c("mean","1simp"),c("mean","3simp")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-01.Rds",sep=''))

# varying samplesize
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
  methods_all = list(c("mean","1simp"),c("mean","3simp")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-03.Rds",sep=''))




# #---- potential fig 3: high-dim DAG setting for varying nsamples ----
# - see for which samplesize our algos achieve optimal SHD
# - see if SID gets better beyong optimal SHD for our algos
# nbh = 1
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11,12),
  interventionSize = c(1),
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
  methods_all = list(c("mean","1simp"),c("mean","3simp"),c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-06.Rds",sep=''))



# nbh = 5
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11,12),
  interventionSize = c(1),
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
  methods_all = list(c("mean","1simp"),c("mean","3simp"),c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-07.Rds",sep=''))



# nbh = 10
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(11,12),
  interventionSize = c(1),
  ndatasets = c(3),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 2
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","1simp"),c("mean","3simp"),c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-08.Rds",sep=''))
