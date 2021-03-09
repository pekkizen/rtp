
# rtp  

Methods for combining p-values.

## P-value functions

**p.fisher(p)**  
Standard Fisher's method using all p-values.  

**p.art(K, p)**  
Augmented RTP function from Vsevolozhskaya et al (2019).  

**p.tfisher.soft(tau, p)**  
Soft thresholding TFisher method from Zhang et al (2020).  

**p.tfisher.tpm(tau, p)**  
Truncated product method (tpm) from Zaykin et al (2002)
solved by TFisher algorithm.  

**p.rtp(K, p, stepscale = 1)**  
Rank Truncated Product method implemented by gamma PDF x beta CDF
integrated by Riemann sums.
This is a fast and robust function selected from
the other functions presented in R/C++ code.

**p.rtp.mutoss(K, p)**  
Rank Truncated Product p-value from
rdrr.io/cran/mutoss/src/R/Rank_Truncated.R. This uses exactly equations
from Dudbridge and Koeleman (2003).

**p.rtp.qbeta(K, p, abstol = 1e-4, reltol = 1e-2)**  
Rank Truncated Product p-value by
inverse beta CDF/Quantile function method from Vsevolozhskaya et al (2019).

## Documentation

R-style documentation is coming soon or later. Statistical and computational methods
are somewhat commented in the program code.

## Installation

Package rtp is build and tested by Rcpp version 1.0.6, R version 4.0.3 and g++
compiler in Rtools 4.0 on Windows system.

```R
install.packages("devtools") # if not installed
library(devtools)
install_github("pekkizen/rtp")
```

## Example

Generate 100 random p-values with one small value 1e-6 and calculate
combined p-value by different methods.

```R
library(rtp)
K <- 10
L <- 100
p <- c(1e-6, runif(L - 1))
tau <- K / L
p.art(K, p)
p.rtp(K, p)
p.tfisher.soft(tau, p)
p.tfisher.tpm(tau, p)
p.fisher(p)
```

## References

[1] Zaykin DV, Zhivotovsky LA, Westfall PH, and Weir BS. (2002).
Truncated product method for combining p-values
Genetic Epidemiology 22: 170-185.
ftp://statgen.ncsu.edu/pub/zaykin/tpm  

[2] Dudbridge, F., and Koeleman, B. P. (2003).
Rank truncated product of p-values, with application to genomewide association scans.
Genet. Epidemiol. 25, 360–366. doi: 10.1002/gepi.10264  

[3] Olga A. Vsevolozhskaya, Fengjiao Hu and Dmitri V. Zaykin. (2019). Detecting Weak Signals by Combining Small P-Values in Genetic Association Studies.
Front. Genet. 10:1051. doi: 10.3389/fgene.2019.01051.  
www.frontiersin.org/articles/10.3389/fgene.2019.01051/full. www.frontiersin.org/articles/10.3389/fgene.2019.01051/full#supplementary-material.  
R code: github.com/dmitri-zaykin/Total_Decor.  

[4] Zhang, Hong; Tong, Tiejun; Landers, John; Wu, Zheyang. TFisher: A powerful truncation and weighting procedure for combining p-values. Ann. Appl. Stat. 14 (2020), no. 1, 178--201. doi:10.1214/19-AOAS1302.  
projecteuclid.org/euclid.aoas/1587002670.
projecteuclid.org/euclid.aoas/1587002670#supplemental.  
R code: CRAN.R-project.org/package=TFisher

https://en.wikipedia.org/wiki/Gamma_distribution  
Characterization using shape α and rate β.

https://en.wikipedia.org/wiki/Beta_distribution  
Probability density function.

https://en.wikipedia.org/wiki/Order_statistic  
Order statistics sampled from a uniform distribution.  
