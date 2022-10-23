# File to setup parallel computation.
# WARNING: if you want to use parallelisation, only execute 01. Once you are done,
#   execute 02


# 01 setup parallel cores
print("opened parallel threads")
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK"
)
doParallel::registerDoParallel(cl = my.cluster)


# 02 close parallel cores
parallel::stopCluster(cl = my.cluster)
print("closed parallel threads")