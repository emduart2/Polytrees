### R Code for paper "Learning Linear Gaussian Polytree Models with Interventions"

This repository contains the code to create the figures of the paper "Learning Linear Gaussian Polytree Models with Interventions" and to further explore the presented algorithms.

For an overview of the functionality for further exploring the algorithms, start at "tutorial.R". Further content of the files:
- figures_from_paper.R: figure 1 and table 1
- sachsdata.R: table 2
- computations_for_figures.R: script to generate the data for figures_from_paper.R. Results can be found in folder data


### Parallelisation
Parallelisation of the exploration functions can only be performed on unix systems! To enable parallelisation, execute
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
