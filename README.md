# rtp

## P-value functions

**p.fisher(p)**  
Standard Fisher's method using all p-values.  

**p.art(K, p)**   
Augmented RTP function from Vsevolozhskaya et al (2019).  

**p.tfisher.soft(tau, p)**   
Soft thresholding TFisher method from Zhang et al (2020)  

**p.tfisher.tpm(tau, p)**   
Truncated product method (tpm) from Zaykin et al (2002).    
Solved by TFisfer algorithm  


**p.rtp(K, p)**   
Rank Truncated Product p-value.  
This is a fast and robust function selected from  
the other functions presented here. The other  
functions are not worse in any significant way.  

**ranktruncated(K, p)**   
Rank Truncated Product p-value  
rdrr.io/cran/mutoss/src/R/Rank_Truncated.R  

**p.rtp.qbeta.integrate(K, p)**    
Rank Truncated Product p-value  
Inverse beta CDF/Quantile function method from Vsevolozhskaya et al (2019).

## References  
[1] Zaykin DV, Zhivotovsky LA, Westfall PH, and Weir BS. (2002).   
Truncated product method for combining p-values.   
Genetic Epidemiology 22: 170-185  
ftp://statgen.ncsu.edu/pub/zaykin/tpm  

[2] Dudbridge, F., and Koeleman, B. P. (2003).   
Rank truncated product of p-values, with application to genomewide association scans.   
Genet. Epidemiol. 25, 360–366. doi: 10.1002/gepi.10264  

[3] Olga A. Vsevolozhskaya, Fengjiao Hu and Dmitri V. Zaykin. (2019).  
Detecting Weak Signals by Combining Small P-Values in Genetic Association Studies.  
Front. Genet. 10:1051. doi: 10.3389/fgene.2019.01051.  
Plus Supplementary Material.  
R: github.com/dmitri-zaykin/Total_Decor. 

[4] Zhang, Hong; Tong, Tiejun; Landers, John; Wu, Zheyang.  
TFisher: A powerful truncation and weighting procedure for combining p-values.  
Ann. Appl. Stat. 14 (2020), no. 1, 178--201. doi:10.1214/19-AOAS1302.  
Plus Supplement to TFisher.  
R: CRAN.R-project.org/package=TFisher.  

https://en.wikipedia.org/wiki/Gamma_distribution. Characterization using shape α and rate β  

https://en.wikipedia.org/wiki/Beta_distribution. Probability density function 

https://en.wikipedia.org/wiki/Order_statistic. Order statistics sampled from a uniform distribution  
