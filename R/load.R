
library(cubature)
library(microbenchmark)
library(Rcpp)

rm(list = ls()) # clear environment
sourceCpp("./src/fun.cpp")
source("./R/tfisher.R")
source("./R/bench.R")
source("./R/plot.R")
source("./R/tests.R")
source("./R/pvalues.R")
source("./R/mutoss.R")
source("./R/rank.R")