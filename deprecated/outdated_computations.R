

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






# ---- orientation high-dimensional ??? -----
# ll1 = readRDS("cluster_computation/03_results_with500/03-05.Rds")
ll1 = readRDS("cluster_computation/04_results/04-05.Rds")$df
df_ssBIC = ll1[ll1$kindOfIntervention == "perfect",]

ll2 = readRDS("cluster_computation/04_results/04-06.Rds")$df

tmp = readRDS("local_computations/01-02.Rds")$df # lowdimensional
tmp$time = tmp$time_s
tmp$time_s = NULL



f1 = ggplot(rbind(dfBIC, dfTESTGIES), aes(totalSamples, SHD, fill=method)) +
  geom_boxplot() + labs(title = "low-dim Polytrees")

ndsBIC = readRDS("cluster_computation/03_results_with500/03-04.Rds")
nssBIC = readRDS("cluster_computation/03_results_with500/03-05.Rds")




nds vary: 1,2,3,3simp * TEST
