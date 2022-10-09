# paramters
kmax = 10
save_dir ="/dss/dsshome1/lxc06/ge45bop2/mkdir/Intervention/Parallel/Results/" # with final '/' in the end
experiment_id = "04"

# -- parallellisation below here-- 
ww<-as.numeric(Sys.getenv("SLURM_PROCID"))+1

p<-500

#---- fig 2 but high-dim ----
# - pick one of interventionsize (whatever look sbetter)
# - analyse the trends of varying totalSamples & ndatasets
# varying ndatases
if(ww==1){
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(500),
  interventionSize = c(2,10),
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
  methods_all = list(c("mean","1simp"),c("mean","3simp")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-01.Rds",sep=''))
}else if(ww==2){
# varying ndatases with GIES
# - check runtime of GIES 
# - if it terminates, we want to say that we are faster
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(500),
  interventionSize = c(2,10),
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
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-01_GIES.Rds",sep=''))
}else if(ww==3){

# varying samplesize
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(200,500,1000),
  interventionSize = c(2,10),
  ndatasets = c(11),
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
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-03.Rds",sep=''))
}else if(ww==4){

# #---- potential fig 3: high-dim DAG setting for varying nsamples ----
# - pick one interventionsize (whatever looks better)
# - compare runtime of algorithms, ours should be faster than GIES
# nbh = 1
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
  dag_nbh = 1
)
l <- explore(
  df_params,
  scoreFct_all = list(pcalg::shd, SID, true_positives, false_positives, true_negatives, false_negatives),
  sFctNames_all = c("SHD","SID","TP","FP","TN","FN"),
  methods_all = list(c("mean","1simp"),c("mean","3simp")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-06.Rds",sep=''))

}
