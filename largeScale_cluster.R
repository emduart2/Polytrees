# TODO add install packages and load packages
source("setup.R")

# notes: 
# - skeletonExploration now includes baseline plots (hence replot fig1)
# - high-dim: 1000 nodes, 51 datasets (5%), interventionsize 1
kmax = 50

#---- fig 1 low-dim with all intervention types ----
p = 20
nSamples = 200
# left
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(70,100,200),
  interventionSize = c(1),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect","imperfect","inhibitory","random"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params)
l001 = l

# middle
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nSamples),
  interventionSize = c(1),
  ndatasets = c(2,11,21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params); 
l002 = l

# right
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(nSamples),
  interventionSize = c(1,2,4),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params);
l003 = l



#---- fig 1 but high-dim ----
# only take one of totalSamples for figure middle&right
p = 1000
# left
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(100,200,500,800,1000),
  interventionSize = c(1),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params)
l004 = l

# middle
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(100,200,500),
  interventionSize = c(1),
  ndatasets = c(2,11,21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params); 
l005 = l

# right
df_params <- expand.grid(
  tsize = c(p),
  totalSamples = c(100,200,500),
  interventionSize = c(1,2,4),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE
)
l <- skeletonExploration(df_params);
l006 = l




#---- fig 2 but high-dim ----
# for left, only take one of totalSamples
# left
df_params <- expand.grid(
  tsize = c(1000),
  totalSamples = c(100,200,500),
  interventionSize = c(2),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05
)
l <- orientationExploration(df_params);
l007 = l

# right
df_params <- expand.grid(
  tsize = c(1000),
  totalSamples = c(100,200,500),
  interventionSize = c(2),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05
)
l <- orientationExploration(df_params);
l008 = l



#---- fig 2 but with all intervention types ----
# left
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(500),
  interventionSize = c(2),
  ndatasets = c(2,11,21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect","imperfect","inhibitory","random"),
  ensureDiff = TRUE,
  alpha = 0.05
)
l <- orientationExploration(df_params);
l009 = l

# right
df_params <- expand.grid(
  tsize = c(20),
  totalSamples = c(100,500,1000),
  interventionSize = c(1),
  ndatasets = c(21),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect","imperfect","inhibitory","random"),
  ensureDiff = TRUE,
  alpha = 0.05
)
l <- orientationExploration(df_params);
l010 = l


#---- potential fig 3: compare GIES to our algo in high-dim for varying nsamples ----
df_params <- expand.grid(
  tsize = c(1000),
  totalSamples = c(100,200,500,1000),
  interventionSize = c(1),
  ndatasets = c(51),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l312 <- bothExploration(df_params,methods);


#---- potential fig 3: compare GIES to our algo in high-dim for varying ndatasets ----
# pick one samplesize
df_params <- expand.grid(
  tsize = c(1000),
  totalSamples = c(100,200,500),
  interventionSize = c(1),
  ndatasets = c(11,51,101),
  k = c(1:kmax),
  sdatasets = list(c()),
  kindOfIntervention = c("perfect"),
  ensureDiff = TRUE,
  alpha = 0.05
)
methods = list(c("GIES"),c("mean","1"),c("mean","3"))
l312 <- bothExploration(df_params,methods)


#---- potential fig 3: high-dim DAG setting for varying nsamples ----
# TODO if we want to do that, we will need to add a DAG flag ...