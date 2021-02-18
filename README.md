# rtp

Standard Fisher's method using all p-values.  
p.fisher <- function(p)
 
p.art is Augmented RTP function from Vsevolozhskaya et al (2019).  
p.art <- function(K, p)

Soft thresholding TFisher method from Zhang et al (2020)
p.tfisher.soft <- function(tau, p)

runcated product method (tpm) from Zaykin et al (2002) by TFisfer.
Solved by TFisfer algorithm
p.tfisher.tpm <- function(tau, p)

p.rtp returns Rank Truncated Product p-value.
This is a fast and robust function selected from
the other functions presented here. The other
functions are not worse in any significant way.
p.rtp <- function(K, p, tol = 1e-10, stepscale = 1)
