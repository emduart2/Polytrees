# notes: 
# - skeletonExploration now includes baseline plots (hence replot fig1)
# - high-dim: 1000 nodes, 21 datasets (2%), interventionsize 1
kmax = 2
save_dir = 'testres/'

#---- fig 1 low-dim with all intervention types ----
p = 3
nSamples = 20
# left
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(10,11),
  interventionSize = c(1),
  ndatasets = c(2),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect","imperfect","inhibitory","random"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"l001.Rds",sep=''))

# middle
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nSamples),
  interventionSize = c(1),
  ndatasets = c(2,3),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect","imperfect","inhibitory","random"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params); 
saveRDS(l, file = paste(save_dir,"l002.Rds",sep=''))

# right
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nSamples),
  interventionSize = c(1,2),
  ndatasets = c(2),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect","imperfect","inhibitory","random"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params);
saveRDS(l, file = paste(save_dir,"l003.Rds",sep=''))



#---- fig 1 but high-dim ----
# only take one of totalSamples for figure middle&right
p = 3
# left
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(10,11),
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
saveRDS(l, file = paste(save_dir,"l004.Rds",sep=''))

# middle
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(10,11),
  interventionSize = c(1),
  ndatasets = c(2,3),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params); 
saveRDS(l, file = paste(save_dir,"l005.Rds",sep=''))

# right
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(10,11),
  interventionSize = c(1,2),
  ndatasets = c(2),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params);
saveRDS(l, file = paste(save_dir,"l006.Rds",sep=''))




#---- fig 2 but high-dim ----
# for left, only take one of totalSamples
# left
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(10,11),
  interventionSize = c(2),
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
  methods_all = list(c("mean","1"),c("mean","1simp"),c("mean","2"),c("mean","2simp"),c("mean","3"),c("mean","3simp")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,"l007.Rds",sep=''))

# right
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(10,11),
  interventionSize = c(2),
  ndatasets = c(2),
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
  methods_all = list(c("mean","1"),c("mean","1simp"),c("mean","2"),c("mean","2simp"),c("mean","3"),c("mean","3simp")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,"l008.Rds",sep=''))



#---- fig 2 but with all intervention types ----
# left
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(10),
  interventionSize = c(2),
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
saveRDS(l, file = paste(save_dir,"l009.Rds",sep=''))

# right
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(10,11),
  interventionSize = c(1),
  ndatasets = c(2),
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
saveRDS(l, file = paste(save_dir,"l010.Rds",sep=''))


#---- potential fig 3: compare GIES to our algo in high-dim for varying nsamples ----
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(10,11),
  interventionSize = c(1),
  ndatasets = c(2),
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
  methods_all = list(c("GIES","GIES"),c("mean","1"),c("mean","3")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,"l011.Rds",sep=''))


#---- potential fig 3: compare GIES to our algo in high-dim for varying ndatasets ----
# pick one samplesize
df_params <- expand.grid(
  tsize = c(3),
  totalSamples = c(10,11),
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
  methods_all = list(c("GIES","GIES"),c("mean","1"),c("mean","3")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,"l012.Rds",sep=''))



# #---- potential fig 3: high-dim DAG setting for varying nsamples ----
# with nbh = 1
df_params <- expand.grid(
  tsize = c(4),
  totalSamples = c(11,12),
  interventionSize = c(1),
  ndatasets = c(2),
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
  methods_all = list(c("GIES","GIES"),c("mean","1"),c("mean","3")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,"l013.Rds",sep=''))

# with nbh = 3
df_params <- expand.grid(
  tsize = c(4),
  totalSamples = c(11),
  interventionSize = c(1),
  ndatasets = c(2),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 3
)
l014 <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("GIES","GIES"),c("mean","1"),c("mean","3")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,"l014.Rds",sep=''))

