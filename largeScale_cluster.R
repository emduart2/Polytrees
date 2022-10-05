# install packages and load packages
source("setup.R")

# notes: 
# - skeletonExploration now includes baseline plots (hence replot fig1)
# - high-dim: 1000 nodes, 21 datasets (2%), interventionsize 1
kmax = 50
save_dir = ''

#---- fig 1 low-dim with all intervention types ----
p = 20
nSamples = 200
# left
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(70,100,200),
  interventionSize = c(1),
  ndatasets = c(21),
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
  ndatasets = c(2,11,21),
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
  interventionSize = c(1,2,4),
  ndatasets = c(21),
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
p = 1000
# left
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(100,200,500,1000),
  interventionSize = c(1),
  ndatasets = c(21),
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
  totalSamples = c(100,200,500,1000),
  interventionSize = c(1),
  ndatasets = c(11,21,31),
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
  totalSamples = c(100,200,500,1000),
  interventionSize = c(1,5,10),
  ndatasets = c(21),
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
  tsize = c(1000),
  totalSamples = c(100,200,500),
  interventionSize = c(2),
  ndatasets = c(11,21,31),
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
  tsize = c(1000),
  totalSamples = c(100,200,500),
  interventionSize = c(2),
  ndatasets = c(21),
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
  tsize = c(20),
  totalSamples = c(500),
  interventionSize = c(2),
  ndatasets = c(2,11,21),
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
  tsize = c(20),
  totalSamples = c(100,500,1000),
  interventionSize = c(1),
  ndatasets = c(21),
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
  tsize = c(1000),
  totalSamples = c(100,200,500,1000),
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
  tsize = c(1000),
  totalSamples = c(100,200,500),
  interventionSize = c(1),
  ndatasets = c(11,21,31),
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
# # with nbh = ?
# df_params <- expand.grid(
#   tsize = c(1000),
#   totalSamples = c(100,200,500),
#   interventionSize = c(1),
#   ndatasets = c(21),
#   k = c(1:kmax),
#   sdatasets = list(c()),
#   kindOfIntervention = c("perfect"),
#   ensureDiff = TRUE,
#   alpha = 0.05,
#   use_dags = TRUE,
#   dag_nbh = 10
# )
# l013 <- explore(
#   df_params,
#   scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives), 
#   sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
#   methods_all = list(c("GIES","GIES"),c("mean","1"),c("mean","3")),
#   pw_methods_all = c("BIC"))
# 
# # with nbh = ?
# df_params <- expand.grid(
#   tsize = c(1000),
#   totalSamples = c(100,200,500),
#   interventionSize = c(1),
#   ndatasets = c(11,51,101),
#   k = c(1:kmax),
#   sdatasets = list(c()),
#   kindOfIntervention = c("perfect"),
#   ensureDiff = TRUE,
#   alpha = 0.05,
#   use_dags = TRUE,
#   dag_nbh = 10
# )
# l014 <- explore(
#   df_params,
#   scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives), 
#   sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
#   methods_all = list(c("GIES","GIES"),c("mean","1"),c("mean","3")),
#   pw_methods_all = c("BIC"))


