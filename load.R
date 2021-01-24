
library(cubature)
library(Rcpp)

rm(list = ls()) # clear environment
sourceCpp("func.cpp")
source("tfisher.R")
source("testfun.R")
source("mutoss.R")
source("rank.R")