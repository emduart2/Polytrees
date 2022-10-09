# load custom functions
source("/dss/dsshome1/lxc06/ge45bop2/mkdir/Intervention/Parallel/polyFunctions.R")
source("/dss/dsshome1/lxc06/ge45bop2/mkdir/Intervention/Parallel/Likelihood.R")
source("/dss/dsshome1/lxc06/ge45bop2/mkdir/Intervention/Parallel/explore_performance.R")
source("/dss/dsshome1/lxc06/ge45bop2/mkdir/Intervention/Parallel/eval_fcts.R")
library(pcalg)
library(patchwork)
library(igraph)
library(sets)
library(mpoly)

# setup parallel cores
print("opened parallel threads")
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK"
)
doParallel::registerDoParallel(cl = my.cluster)


# execute main file
print("started main execution")
s = Sys.time()
source("/dss/dsshome1/lxc06/ge45bop2/mkdir/Intervention/Parallel/04/04_part1.R")
(time_elapsed = Sys.time() - s)
print("finished main execution")

# close parallel cores
parallel::stopCluster(cl = my.cluster)
print("closed parallel threads")