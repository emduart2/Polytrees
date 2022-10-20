# This script contains the function calls to generate the data for the 
# figures of the supplementary  material of this paper. 
# WARNING: The execution of this script takes long. It is recommended to not
#   execute this script without utilizing multiple parallelized cores.
save_dir = "data/"

# ---- runtime skel rec ----
# left
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(200,500,1000),
  interventionSize = c(1),
  ndatasets = c(20),
  k = c(1:10),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"1a_01_runtime.Rds",sep=""))

# middle``
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(500),
  interventionSize = c(1),
  ndatasets = c(2,11,21),
  k = c(1:10),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"1a_02_runtime.Rds",sep=""))

# right
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(500),
  interventionSize = c(1,2,4),
  ndatasets = c(20),
  k = c(1:10),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"1a_03_runtime.Rds",sep=""))


# ---- low-dim pooled ----
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(500),
  interventionSize = c(1,2,4,6,8,10),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect","flipped"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"supp_pooled.Rds",sep=""))


# ---- additional comp for different kinds of interventions ----
df_params <- expand.grid(
  tsize = c(500),
  totalSamples = c(200),
  interventionSize = c(10),
  ndatasets = c(21),
  k = c(1:20),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect","flipped"),
  ensureDiff = TRUE,
  use_dags = FALSE,
  dag_nbh = 0
)
l <- explore_skeleton(df_params)
saveRDS(l, file = paste(save_dir,"supp_pooled.Rds",sep=""))