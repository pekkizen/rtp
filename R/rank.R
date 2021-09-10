
# Standard Fisher's method using all p-values.
#  Solved by Gamma distribution. RTP for K == L.
p.fisher <- function(p) {
    n <- length(p)
    e <- exp(1)
    lw <- log(prod(p * e)) - n
    if (is.infinite(lw)) {
        lw <- sum(log(p))
    }
    pgamma(-lw, n, lower.tail = FALSE)
}

# single returns the probability of getting by random
# one or more p-values <= minimum observed p-value.
# RTP for K == 1.
p.single <- function(p) {
    L <- length(p)
    return(pbinom(0, L, min(p), lower.tail = FALSE))
}

# p.art is slightly modified Art function in Vsevolozhskaya et al.
p.art <- function(K, q) {
    L <- length(q)
    p <- sort(q, partial = c(1:K))[1:K]
    lw <- sum(log(p[1:K - 1]))
    pk <- p[K]
    d <- (K - 1) * (digamma(L + 1) - digamma(K))
    q <- qgamma(1 - pbeta(pk, K, L - K + 1), shape = d)
    ak <- (K - 1) * log(pk) - lw + q
    1 - pgamma(ak, shape = K + d - 1)
}

# Reference integration by library cubature function pcubature.
# The integral is divided into two parts at the approx.
# highest point of the integrand fBetaD.
p.rtp.dbeta.cuba <- function(K, p, tol = 1e-4) {
    r <- init(K, p, 1) # Rcpp
    if (r <= 1) {
        return(r)
    }
    top <- K / length(p) * 0.5 # very approx.

    I <- cubature::pcubature(fBetaD, 0, top, tol = tol)$integral
    I <- I + cubature::pcubature(fBetaD, top, 1, tol = tol)$integral
    min(1, I)
}

# RPT p-value by inverse Beta CDF and R integrate function.
# The inverse beta CDF function method from Vsevolozhskaya et al.
p.rtp.qbeta <- function(K, p, tol = 1e-4) {
    L <- length(p)
    lw <- stat.rtp(K, p) # Rcpp
    p.rtp.qbeta.lw(lw, K, L, tol = tol)
}

p.rtp.qbeta.lw <- function(lw, K, L, tol = 1e-5) {
    f <- function(x) fBetaQ.R(x, lw, K, L)
    if (tol < 1e-12) tol <- 1e-12
    integrate(f, 0, 1, rel.tol = tol)$value
}

# RPT p-value by simulation.
p.rtp.simulated <- function(K, p, rounds = 5e5, seed = 0) {
    if (seed > 0) set.seed(seed) # goes to Rcpp
    rtpSimulated(K, p, rounds) # Rcpp
}

meanRTP <- function(K, L) {
    K - K * (digamma(K + 1) - digamma(L + 1))
}
varRTP <- function(K, L) {
    K + (K * K) * (trigamma(K + 1) - trigamma(L + 1))
}