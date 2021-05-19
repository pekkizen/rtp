
# rtp  

Methods for combining p-values.

## P-value functions

| Parameters |                                     |
|-----------|----------------------------------------|
|p          |Vector of p-values.                      |
|K          |Number of smallest p-values in RTP 
|tau        |Threshold p-value for TFisher methods.  |
|stepscale  |Multiplier for default integration      |  
|           |                                     |

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
Rank truncated product p-value by integrand gamma PDF x beta CDF or beta PDF x (1-gamma CDF) and Riemann sum integral over [0, inf) or [0, 1].
The second integrand is used for  K / number of all p-values proportion greater than 1/3.
Parameter stepscale affects the integration step size.
E.g. stepscale = 0.5 halves the default step size and increases accuracy (which mostly is unnecessary).

**p.rtp.simulated(K, p, rounds = 100000, seed = 0)**  
Rank truncated product p-value by Monte Carlo simulation. This is a slow function and not good for actual use, but the program code in the rtp.cpp file demonstrates what rtp p-values are really measuring.
seed is used to initialize random number generator. If 0, randomly seeded.

### Reference implementations of Rank Truncated Product method

**p.rtp.mutoss(K, p)**  
Rank truncated product p-value from
https://rdrr.io/cran/mutoss/src/R/Rank_Truncated.R. This implements exactly equations from Dudbridge and Koeleman. Integrated over [0, 1] by 1000 equidistant function evaluations.

**p.rtp.qbeta(K, p, abstol = 1e-4, reltol = 1e-2)**  
Rank truncated product p-value by inverse beta CDF method from Vsevolozhskaya et al.
Integrated over [0, 1] by R integrate.

**p.rtp.dbeta.cuba(K, p)**  
Rank truncated product p-value by integrand beta PDF x (1 - gamma CDF) integrated by package cubature function pcubature over [0, 1] with very small tolerances.
The integrand implements Dudbridge and Koeleman equations.
This gives high accuracy "exact" reference values very reliably, but quite slow.

## Documentation

R-style documentation may come soon or later. Statistical and computational methods are somewhat commented in the program code.

## Installation

Package rtp is built and tested with R version 4.0.4, Rcpp version 1.0.6 and
g++ compiler version 8.3.0 in Rtools 4.0 on Windows system.

```R
install.packages("devtools") # if not installed
library(devtools)
install_github("pekkizen/rtp")
library(rtp)
```

## Examples

```R
# Generate 199 random p-values, add one small value 1e-7 and  
# calculate combined p-value by different methods.
K <- 10       # 10 smallest p-values
L <- 200
p <- c(1e-7, runif(L - 1))
tau <- K / L  # p-values less than 0.05
p.art(K, p)
p.rtp(K, p)
p.tfisher.soft(tau, p)
p.tfisher.tpm(tau, p)
p.fisher(p)

# Comparing methods has also a dedicated function
pvalues.methods(K = 10, L = 200, small = 1e-7)

# Comparing accuracy of the rtp functions
pvalues.rtp(K=10, L=200, small=1e-6, seed=0)

# The four integrand functions can be examined by
plot.quantile.integrands(K = 10, L = 200, small = 1e-1)
plot.BxG.integrand(K = 10, L = 200, small = 1e-3)
plot.GxB.integrand(K = 10, L = 200, small = 1e-3)

# Integrand and p-value functions benchmarks
bench.integrands(K = 10, L = 200, small = 1e-5)
bench.pvalues(K = 10, L = 200, small = 1e-5)
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
The probability density function.
The cumulative distribution function for positive
integers α and β = 1.

https://en.wikipedia.org/wiki/Beta_distribution  
Probability density function.

https://en.wikipedia.org/wiki/Beta_function  
https://dlmf.nist.gov/8.17#E5. Formula 8.17.5  
Regularized incomplete beta function = beta CDF function.
