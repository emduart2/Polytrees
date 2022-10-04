# This file contains functions for evaluating the performance of predictions

# computes the SHD of two skeletons given as (symmetric) adjacency
# matrices. If number of nodes larger three and the skeletons are both
# trees, the result lies in [0,1]
SHD_skeleton <- function(S1,S2){
  return(sum(abs(S1-S2))/2)
}

# computes SID of two DAGs given as (symmetric) adjacency matrices. 
SID <- function(S1,S2){
  return(SID::structIntervDist(S1,S2)$sidUpperBound)
}