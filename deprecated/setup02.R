# load custom functions
source("polyFunctions.R")
source("Likelihood.R")
source("explore_performance.R")
source("eval_fcts.R")
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
source("02_cluster.R")
(time_elapsed = Sys.time() - s)
print("finished main execution")

# close parallel cores
parallel::stopCluster(cl = my.cluster)
print("closed parallel threads")