save_dir = ""
p = 500
nsamples_solo = 500
nsamples_range = c(200,500,1000)
nds_solo = 21
nds_range = c(2,11,21)
kmax = 10
intervSize = 1
start = Sys.time()


#---- high-dim figure 2 -----
# varying ndatasets
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nsamples_solo),
  interventionSize = c(intervSize),
  ndatasets = c(nds_range),
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
  scoreFct_all = list(pcalg::shd, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","TP","FP","TN","FN"),
  methods_all = list(c("gtruth","1simp"),c("gtruth","1"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
  pw_methods_all = c("TEST","BIC"))
saveRDS(l1, file = paste(save_dir,"02-01.Rds",sep=""))

# varying ndatasets GIES
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nsamples_solo),
  interventionSize = c(intervSize),
  ndatasets = c(nds_range),
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
  scoreFct_all = list(pcalg::shd, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","TP","FP","TN","FN"),
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l1, file = paste(save_dir,"02-01_GIES.Rds",sep=""))


# varying nsamples
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nsamples_range),
  interventionSize = c(intervSize),
  ndatasets = c(nds_solo),
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
  scoreFct_all = list(pcalg::shd, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","TP","FP","TN","FN"),
  methods_all = list(c("gtruth","1simp"),c("gtruth","1"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
  pw_methods_all = c("TEST","BIC"))
saveRDS(l1, file = paste(save_dir,"02-02.Rds",sep=""))


# varying nsamples GIES
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nsamples_range),
  interventionSize = c(intervSize),
  ndatasets = c(nds_solo),
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
  scoreFct_all = list(pcalg::shd, true_positives, false_positives, true_negatives, false_negatives), 
  sFctNames_all = c("SHD","TP","FP","TN","FN"),
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l1, file = paste(save_dir,"02-02_GIES.Rds",sep=""))




#---- high-dim figure 1 ----
# left
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nsamples_range),
  interventionSize = c(1),
  ndatasets = c(nds_solo),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"02-03.Rds",sep=""))

# middle
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nsamples_solo),
  interventionSize = c(1),
  ndatasets = c(nds_range),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"02-04.Rds",sep=""))

# right
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nsamples_solo),
  interventionSize = c(1,2,4),
  ndatasets = c(nds_solo),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"02-05.Rds",sep=""))






Sys.time() - start