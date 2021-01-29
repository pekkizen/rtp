
library(cubature)
library(Rcpp)

# rm(list = ls()) # clear environment
sourceCpp("fun.cpp")
source("tfisher.R")
source("testfun.R")
source("mutoss.R")
source("rank.R")