
# This file has the four possible integrands as R code.
# Gamma functions use default rate = 1.

# fBetaQuantile.R is rtp integrand over Beta probabilities in [0, 1].
# This is the integrand in Vsevolozhskaya et al.
fBetaQ.R <- function(p, lw, K, L) {
    b <- qbeta(p, K + 1, L - K)
    g <- log(b) * K - lw
    1 - pgamma(g, K)
    # pgamma(g, K, lower.tail = FALSE) # better accuracy for small p's
}

# fGammaQ.R is rtp integrand over Gamma probabilities in [0, 1].
# fGammaQ.R(1 - p, ) is close to fBetaQ.R(p, ).
fGammaQ.R <- function(p, lw, K, L) {
    g <- qgamma(p, K)
    b <- exp((g + lw) / K)
    pbeta(b, K + 1, L - K)
}

# fBetaD is rtp integrand Beta PDF x (1 - Gamma CDF) over b in [0, 1].
# This implement the equations in Dudbridge and Koeleman.
fBetaD.R <- function(b, lw, K, L) {
    g <- K * log(b) - lw
    dbeta(b, K + 1, L - K) * pgamma(g, K, lower.tail = FALSE)
}

# fGammaD is rtp integrand Gamma PDF x Beta CDF over g in [0, inf).
fGammaD.R <- function(g, lw, K, L) {
    b <- exp((g + lw) / K)
    dgamma(g, K) * pbeta(b, K + 1, L - K)
}