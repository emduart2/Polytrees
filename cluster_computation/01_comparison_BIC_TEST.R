# figure 2 with comparison BIC vs TEST
# left
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
  methods_all = list(c("gtruth","1"),c("gtruth","1simp"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
  pw_methods_all = c("BIC","TEST"))
saveRDS(l1, file = "l101")


# right
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
  methods_all = list(c("gtruth","1"),c("gtruth","1simp"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
  pw_methods_all = c("BIC","TEST"))
saveRDS(l2, file="l102")