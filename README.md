### R Code for paper "Learning Linear Gaussian Polytree Models with Interventions"

This repository contains the code implementing the methods of the paper "Learning Linear Gaussian Polytree Models with Interventions". It serves two purposes: (1) providing the code for generating the Figures&Tables of the paper and (2) providing functionality to further explore the properties of the algorithms and to compare it to other methods. To get an overview of (2), start at "tutorial.R"!

For (1):
- figures_main.R: Figure 1 and Table 1
- sachsdata.R: Table 2
- figures_suppl.R: Figures of the supplementary material

For (2):
- tutorial.R




### Parallelisation
Parallelisation of the exploration functions can currently only be performed on unix systems! To enable parallelisation, execute
```
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK"
)
doParallel::registerDoParallel(cl = my.cluster)
```

To disable parallelisation, execute
```
parallel::stopCluster(cl = my.cluster)
```
This code can also be found in "setup_parallel.R".
