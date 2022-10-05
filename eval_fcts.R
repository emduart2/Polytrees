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

# true positives
true_positives <- function(gtrue,gest){
  adj_true = as(gtrue,"matrix")
  adj_est = as(gest,"matrix")
  return(sum((adj_est==1) & (adj_true==1)))
}

# false positives
false_positives <- function(gtrue,gest){
  adj_true = as(gtrue,"matrix")
  adj_est = as(gest,"matrix")
  return(sum((adj_est==1) & (adj_true==0)))
}


# true negatives
true_negatives <- function(gtrue,gest){
  adj_true = as(gtrue,"matrix")
  adj_est = as(gest,"matrix")
  return(sum((adj_est==0) & (adj_true==0)))
}

# false negatives
false_negatives <- function(gtrue,gest){
  adj_true = as(gtrue,"matrix")
  adj_est = as(gest,"matrix")
  return(sum((adj_est==0) & (adj_true==1)))
}



# optional, might also be calculated in the dataframe from the 4 above
# precision
precision <- function(gtrue,est){
  adj_true = as(gtrue,"matrix")
  adj_est = as(gest,"matrix")
  tp = sum((adj_est==1) & (adj_true==1))
  fp = sum((adj_est==1) & (adj_true==0))
  return(tp / (tp + fp))
}

# recall
recall <- function(gtrue,est){
  adj_true = as(gtrue,"matrix")
  adj_est = as(gest,"matrix")
  tp = sum((adj_est==1) & (adj_true==1))
  p = sum(adj_true==1)
  return(tp / p)
}