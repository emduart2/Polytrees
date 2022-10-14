

# ---- all: low dimensional DAGs ----
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
saveRDS(l, file = "data/lowdim_dags.Rds")


# ---- all: low dimensional polytrees ----
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
saveRDS(l, file = "data/lowdim_polytrees.Rds")
