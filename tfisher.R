
# tfisher R-functions are adapted from Zhang et al (2020) /
# CRAN.R-project.org/package=TFisher
# Here are R functions and equivalent C++ accelerated functions.

# stat.tfisher returns TFisher statistic.
stat.tfisher <- function(p, tau1, tau2) {
    sum(-2 * log(p[p <= tau1] / tau2))
}

# Package TFisher p.tfisher for independant p-values.
# This returns the right tail.
# Equivalent (~14 digits) to function fun.cpp/pTFisher.
p.tfisherR <- function(q, L, tau1, tau2) {
    lt <- log(tau1 / tau2)
    p <- sum(pchisq(q + 2 * (1:L) * lt, 2 * (1:L)) * dbinom(1:L, L, tau1))
    if (q > 0) p <- p + dbinom(0, L, tau1)
    1 - p
}

# Soft thresholding TFisher method.
p.tfisher.soft <- function(tau, p, tol = 1e-14) {
    lw <- stat.tfisher(p, tau, tau)
    L <- length(p)
    tfisher(lw, L, tau, tau, tol)
}

# Soft thresholding TFisher method by pure R.
p.tfisher.softR <- function(tau, p) {
    lw <- stat.tfisher(p, tau, tau)
    L <- length(p)
    p.tfisherR(lw, L, tau, tau)
}

# Truncated product method (tpm) in Zaykin et al (2002) by TFisfer.
p.tfisher.tpm <- function(tau, p) {
    lw <- stat.tfisher(p, tau, 1)
    L <- length(p)
    tfisher(lw, L, tau, 1)
}

# Truncated product method by pure R.
p.tfisher.tpmR <- function(tau, p) {
    lw <- stat.tfisher(p, tau, 1)
    L <- length(p)
    p.tfisherR(lw, L, tau, 1)
}