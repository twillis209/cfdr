source('functions.R')
library(data.table)
library(parallel)
library(tictoc)
library(profvis)

# Using aster_pid data with p-value threshold 1e-7 
load('aster_pid_for_vl.RData')

set.seed(27)

p<-profvis({
  v<-vl(p,q, indices=ind[[1]], mode=2, fold=folds[[1]])
})
