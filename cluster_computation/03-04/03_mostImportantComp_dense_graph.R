# paramters
kmax = 10
save_dir ="/dss/dsshome1/lxc06/ge45bop2/mkdir/Intervention/Parallel/Results/" # with final '/' in the end
experiment_id = "04"

# -- parallellisation below here--
p<-500

ww<-as.numeric(Sys.getenv("SLURM_PROCID"))+1

if(ww==1){
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(500),
  interventionSize = c(1,10),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 20
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","3simp")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-10.Rds",sep=''))
}else if(ww==2){

# right
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(500),
  interventionSize = c(1,10),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 20
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-10_GIES.Rds",sep=''))
}else if(ww==3){

# nbh = 30, GIES
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(500),
  interventionSize = c(1,10),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 30
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","3simp")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-11.Rds",sep=''))
}else if(ww==3){
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(500),
  interventionSize = c(1,10),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05,
  use_dags = TRUE,
  dag_nbh = 30
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-11_GIES.Rds",sep=''))
}