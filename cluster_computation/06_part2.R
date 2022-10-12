# paramters
kmax = 10
save_dir = ''  # with final '/' in the end
experiment_id = "06"

# -- parallellisation below here-- 
ww<-as.numeric(Sys.getenv("SLURM_PROCID"))+1

if(ww == 1){
  # right
  df_params <- expand.grid(
    tsize = c(500),
    totalSamples = c(200,500,1000),
    interventionSize = c(10),
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
    methods_all = list(c("gtruth","1"),c("gtruth","1simp"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
    pw_methods_all = c("BIC","TEST"))
  saveRDS(l, file = paste(save_dir,experiment_id,"-07.Rds",sep=''))
}
