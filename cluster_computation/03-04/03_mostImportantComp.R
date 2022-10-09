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
saveRDS(l, file = paste(save_dir,experiment_id,"-04.Rds",sep=''))
}else if(ww==2){

# right
df_params <- expand.grid(
  tsize = c(p),
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
saveRDS(l, file = paste(save_dir,experiment_id,"-05.Rds",sep=''))
}else if(ww==3){

# nbh = 10, GIES
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
  dag_nbh = 10
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-08_GIES.Rds",sep=''))
}