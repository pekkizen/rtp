
# tfisher R-functions are adapted from Zhang et al (2020).
# Here are R functions and equivalent C/C++ accelerated functions.

# Returns TFisher statistic.
stat.tfisher <- function(p, tau1, tau2) {
    sum(-2 * log(p[p <= tau1] / tau2))
}

# Package TFisher p.tfisher for independant p-values. Zhang et al (2020).
p.tfisher <- function(q, L, tau1, tau2) {
    lt <- log(tau1 / tau2)
    p <- sum(pchisq(q + 2 * (1:L) * lt, 2 * (1:L)) * dbinom(1:L, L, tau1))
    if (q > 0) p <- p + dbinom(0, L, tau1)
    1 - p
}

# Standard Fisher's method using all p-values.
p.fisher <- function(p) {
    lw <- sum(log(p))
    L <- length(p)
    gammaSurv(-lw, L)
}

p.fisher.R <- function(p) {
    lw <- sum(log(p))
    L <- length(p)
    1 - pchisq(-2 * lw, 2 * L)
}

# Soft thresholding TFisher's method in Zhang et al (2020).
p.tfisher.soft <- function(tau, p) {
    lw <- stat.tfisher(p, tau, tau)
    L <- length(p)
    pTFisherCpp(lw, L, tau, tau, tol = 1e-14)
}

p.tfisher.soft.R <- function(tau, p) {
    lw <- stat.tfisher(p, tau, tau)
    L <- length(p)
    p.tfisher(lw, L, tau, tau)
}

# Truncated product method (tpm) in Zaykin et al (2002).
p.tfisher.tpm <- function(tau, p) {
    lw <- stat.tfisher(p, tau, 1)
    L <- length(p)
    pTFisherCpp(lw, L, tau, 1, tol = 1e-14)
}

p.tfisher.tpm.R <- function(tau, p) {
    lw <- stat.tfisher(p, tau, 1)
    L <- length(p)
    p.tfisher(lw, L, tau, 1)
}