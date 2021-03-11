
# rtp  

Methods for combining p-values.

## P-value functions

**p.fisher(p)**  
Standard Fisher's method using all p-values.  

**p.art(K, p)**  
Augmented truncated product function from Vsevolozhskaya et al.  

**p.tfisher.soft(tau, p)**  
Soft thresholding TFisher method from Zhang et al.
This is implemented by a fast C++ TFisher function.

**p.tfisher.tpm(tau, p)**  
Truncated product method (tpm) from Zaykin et al.
solved by TFisher function.  

**p.rtp(K, p, stepscale = 1)**  
Rank truncated product p-value by integrand gamma PDF x beta CDF integrated by Riemann sum.
This is a fast, reliable and accurate, enough or more,  function selected from all functions in the C++ source code.

### Reference implementations of p.rtp

**p.rtp.mutoss(K, p)**  
Rank truncated product p-value from
https://rdrr.io/cran/mutoss/src/R/Rank_Truncated.R. This implements
exactly equations from Dudbridge and Koeleman.

**p.rtp.qbeta(K, p, abstol = 1e-4, reltol = 1e-2)**  
Rank truncated product p-value by inverse beta CDF method
from Vsevolozhskaya et al.

**p.rtp.dbeta.cuba(K, p)**  
Rank truncated product p-value by integrand beta PDF x (1 - gamma CDF) integrated by library cubature function pcubature.
The integrand implements Dudbridge and Koeleman equations.
This gives high accuracy "exact" reference values very reliably, but quite slow.

## Documentation

R-style documentation is coming soon or later. Statistical and computational methods are somewhat commented in the program code.

## Installation

Package rtp is built and tested with R version 4.0.3, Rcpp version 1.0.6 and
g++ compiler version 8.3.0 in Rtools 4.0 on Windows system.

```R
install.packages("devtools") # if not installed
library(devtools)
install_github("pekkizen/rtp")
```

## Examples

Generate 199 random p-values, add one small value 1e-7 and calculate
combined p-value by different methods.

```R
K <- 10       # 10 smallest p-values
L <- 200
p <- c(1e-7, runif(L - 1))
tau <- K / L  # p-values less than 0.05
p.art(K, p)
p.rtp(K, p)
p.tfisher.soft(tau, p)
p.tfisher.tpm(tau, p)
p.fisher(p)

# Comparing methods has also dedicated function
pvalues.methods(K = 10, L = 200, small = 1e-7, seed = 0)

# The integrand functions can be examined by
plot.quantile.integrands(K = 10, L = 200, small = 1e-1, seed = 0)
plot.BxG.integrand(K = 10, L = 200, small = 1e-3, seed = 0)
plot.GxB.integrand(K = 10, L = 200, small = 1e-3, seed = 0)
```

## References

Zaykin DV, Zhivotovsky LA, Westfall PH and Weir BS. (2002).
Truncated product method for combining p-values
Genetic Epidemiology 22: 170-185. ftp://statgen.ncsu.edu/pub/zaykin/tpm  

Dudbridge F. and Koeleman B. P. (2003).
Rank truncated product of p-values, with application to genomewide association scans.
Genet. Epidemiol. 25, 360–366. doi: 10.1002/gepi.10264  

Olga A. Vsevolozhskaya, Fengjiao Hu and Dmitri V. Zaykin. (2019). Detecting Weak
Signals by Combining Small P-Values in Genetic Association Studies.
Front. Genet. 10:1051. doi: 10.3389/fgene.2019.01051.
https://github.com/dmitri-zaykin/Total_Decor.  

Zhang, Hong; Tong, Tiejun; Landers, John; Wu, Zheyang. TFisher: A powerful truncation
and weighting procedure for combining p-values. Ann. Appl. Stat. 14 (2020), no. 1, 178--201.
doi:10.1214/19-AOAS1302.
https://CRAN.R-project.org/package=TFisher

https://en.wikipedia.org/wiki/Order_statistic  
Order statistics sampled from a uniform distribution.  

https://en.wikipedia.org/wiki/Gamma_distribution  
Characterization using shape α and rate β.

https://en.wikipedia.org/wiki/Beta_distribution  
Probability density function.

https://en.wikipedia.org/wiki/Beta_function  Incomplete_beta_function
