library(foreach)
library(doParallel)

# setup parallel cores
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK"
)
doParallel::registerDoParallel(cl = my.cluster)

# parameter
k = 10

# with parallel
s = Sys.time()
foreach(i = 1:k) %dopar% sum(sort(runif(1e7)));
withParallel = Sys.time() - s

# without parallel
s = Sys.time()
for(i in 1:k) {sum(sort(runif(1e7)))};
withoutParallel = Sys.time() - s

# print results
cat("With parallelisation: ",withParallel)
cat("Without parallelisation:", withoutParallel)

# close parallel cores
parallel::stopCluster(cl = my.cluster)



library(foreach)
library(doMPI)

cl <- startMPIcluster()  # use verbose = TRUE for detailed worker message output
registerDoMPI(cl)

system.time(
  foreach(i = 1:167) %dopar% sum(sort(runif(1e7)))
)

closeCluster(cl)
mpi.quit()



p = 10
nds = 3
iss = 1
sd= c()
totalSmpl = 50
ensureDiff=TRUE


