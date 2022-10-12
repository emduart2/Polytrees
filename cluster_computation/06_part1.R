# paramters
kmax = 10
save_dir = ''  # with final '/' in the end
experiment_id = "06"

# -- parallellisation below here-- 
ww<-as.numeric(Sys.getenv("SLURM_PROCID"))+1


if(ww == 1){
  #---- figure 2 low-dim with correct runtime and diff intervSize -----
  # - compare runtime of simple vs refined
  # - compare how TEST behaves for different interventionsizes
  # left
  df_params <- expand.grid(
    tsize = c(20),
    totalSamples = c(200),
    interventionSize = c(1,2,4),
    ndatasets = c(2,11,21),
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
  saveRDS(l, file = paste(save_dir,experiment_id,"-01.Rds",sep=''))
} 

if(ww == 2){
  # right
  df_params <- expand.grid(
    tsize = c(20),
    totalSamples = c(70,200,500),
    interventionSize = c(1,2,4),
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
  saveRDS(l, file = paste(save_dir,experiment_id,"-02.Rds",sep=''))
}





# these are not that computationally difficult
if(ww == 3){
  #----- figure 1 high-dim -------
  # left
  df_params <- expand.grid(
    tsize = c(500),
    totalSamples = c(200,500,1000),
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
  saveRDS(l, file = paste(save_dir,experiment_id,"-03.Rds",sep=''))
  
  # middle
  df_params <- expand.grid(
    tsize = c(500),
    totalSamples = c(500),
    interventionSize = c(1),
    ndatasets = c(2,11,21),
    k = c(1:kmax),
    sdatasets = list(c()),
    kindOfIntervention = c("perfect"),
    ensureDiff = TRUE,
    use_dags = FALSE,
    dag_nbh = 0
  )
  l <- explore_skeleton(df_params)
  saveRDS(l, file = paste(save_dir,experiment_id,"-04.Rds",sep=''))
  
  # right
  df_params <- expand.grid(
    tsize = c(500),
    totalSamples = c(500),
    interventionSize = c(1,2,4),
    ndatasets = c(21),
    k = c(1:kmax),
    sdatasets = list(c()),
    kindOfIntervention = c("perfect"),
    ensureDiff = TRUE,
    use_dags = FALSE,
    dag_nbh = 0
  )
  l <- explore_skeleton(df_params)
  saveRDS(l, file = paste(save_dir,experiment_id,"-05.Rds",sep=''))
}


if(ww == 4){
  # ------- figure 2 high-dim ----------
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
    methods_all = list(c("gtruth","1"),c("gtruth","1simp"),c("gtruth","2"),c("gtruth","2simp"),c("gtruth","3"),c("gtruth","3simp")),
    pw_methods_all = c("BIC","TEST"))
  saveRDS(l, file = paste(save_dir,experiment_id,"-06.Rds",sep=''))
}

