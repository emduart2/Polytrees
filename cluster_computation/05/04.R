# paramters
kmax = 10
save_dir ="/dss/dsshome1/lxc06/ge45bop2/mkdir/Intervention/Parallel/Results/" # with final '/' in the end
experiment_id = "05"

# -- parallellisation below here-- 


ww<-as.numeric(Sys.getenv("SLURM_PROCID"))+1

if(ww==1){

#---- fig 2 but high-dim ----
# compute values for GIES which we forgot for one setting in computation 03
# - analyse the trends of varying totalSamples & ndatasets
# varying samplesize
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(200,500,1000),
  interventionSize = c(10),
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
  methods_all = list(c("GIES","GIES")),
  pw_methods_all = c("TEST"))
saveRDS(l, file = paste(save_dir,experiment_id,"-03_GIES.Rds",sep=''))
}else if(ww==2){

#---- optional: fig2 but with BIC ----
# - compare against TEST from computation03 in terms of runtime, SHD, SID
# varying ndatases
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(500),
  interventionSize = c(10),
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
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-01.Rds",sep=''))

}else if(ww==3){

# varying samplesize
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(200,500,1000),
  interventionSize = c(10),
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
  pw_methods_all = c("BIC"))
saveRDS(l, file = paste(save_dir,experiment_id,"-03.Rds",sep=''))

}else if(ww==4){

# #---- potential fig 3: high-dim DAG setting for varying nsamples ----
# - see for which samplesize our algos achieve optimal SHD
# - see if SID gets better beyong optimal SHD for our algos
# nbh = 1
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(500,1000,2000,3000),
  interventionSize = c(10),
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

