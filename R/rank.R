
# Standard Fisher's method using all p-values.
#  Solved by Gamma distribution. RTP for K == L.
p.fisher <- function(p) {
    lw <- sum(log(p))
    L <- length(p)
    pgamma(-lw, L, lower.tail = FALSE)
    # pchisq(-lw * 2, L * 2, lower.tail = FALSE)
}

# sidak returns the probability of getting one or
# more p-values = minimum observed p-value.
# RTP for K == 1.
p.sidak <- function(p) {
    L <- length(p)
    return(pbinom(0, L, min(p), lower.tail = FALSE))
}

# p.art is slightly modified Art function in Vsevolozhskaya et al.
p.art <- function(K, p) {
    L <- length(p)
    p <- sort(p, partial = c(1:K))[1:K]
    lw <- sum(log(p[1:K - 1]))
    pk <- p[K]
    d <- (K - 1) * (digamma(L + 1) - digamma(K))
    q <- qgamma(1 - pbeta(pk, K, L - K + 1), shape = d)
    ak <- (K - 1) * log(pk) - lw + q
    1 - pgamma(ak, shape = K + d - 1)
}

# p.rtp returns Rank Truncated Product p-value.
p.rtp <- function(K, p, tol = 1e-10, stepscale = 1) {
    q <- p
    rtpRiema(K, q, tol, stepscale)
}

# Reference integration by library cubature function pcubature.
# The integration is divided into two parts at the approx.
# highest point of the integrand fBetaD.
p.rtp.dbeta.cuba <- function(K, p) {
    L <- length(p)
    if (K == 1) {
        return(p.sidak(p))
    }
    if (K == L) {
        return(p.fisher(p))
    }
    e <- init(K, p, 1) # Rcpp
    if (e <= 1) {
        return(e)
    }
    top <- fBetaDtop() # Rcpp
    tol <- 1e-14
    I <- cubature::pcubature(fBetaD, 0, top, tol = tol)$integral # fBetaD Rcpp
    I + cubature::pcubature(fBetaD, top, 1, tol = tol)$integral
}

# RPT p-value by Beta density and adaptive Simpson's 1/3 rule.
p.rtp.dbeta.asimp <- function(K, p, abstol = 1e-8, reltol = 1e-4) {
    rtpDbetaAsimp(K, p, abstol, reltol)
}

# RPT p-value by Beta density and Riemann sum integral.
p.rtp.dbeta.riema <- function(K, p, tol = 1e-10, stepscale = 1) {
    rtpDbetaRiema(K, p, tol, stepscale)
}

# RPT p-value by Gamma density and fixed step Simpson's 1/3 rule.
p.rtp.dgamma.simp <- function(K, p, tol = 1e-10, stepscale = 1) {
    rtpDgammaSimp(K, p, tol, stepscale)
}

# RPT p-value by Gamma density and Riemann sum integral.
p.rtp.dgamma.riema <- function(K, p, tol = 1e-10, stepscale = 1) {
    rtpDgammaRiema(K, p, tol, stepscale)
}

# stat.rpt returns RPT method test statistic.
stat.rtp <- function(K, p) {
    sum(log(sort(p, partial = c(1:K))[1:K]))
}

# RPT p-value by inverse Beta CDF and R integrate function.
# The inverse beta CDF function method from Vsevolozhskaya et al.
p.rtp.qbeta <- function(K, p, abstol = 1e-6, reltol = 1e-3) {
    L <- length(p)
    lw <- stat.rtp(K, p)
    f <- function(u) fBetaQ.R(u, lw, K, L)

    integrate(f, 0, 1, abs.tol = abstol, rel.tol = reltol)$value
}

beta.skewness <- function(a, b) {
    2 * (b - a) * sqrt(a + b + 1) / ((a + b + 2) * sqrt(a * b))
}

# RPT p-value by simulation.
p.rtp.simulated <- function(K, p, rounds = 100000, seed = 0) {
    if (seed > 0) set.seed(seed)
    rtpSimulated(K, p, rounds)
}