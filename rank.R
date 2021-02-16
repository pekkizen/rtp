
# Gamma functions use default rate = 1.

# Standard Fisher's method using all p-values.
p.fisher <- function(p) {
    lw <- sum(log(p))
    L <- length(p)
    pgamma(-lw, L, lower.tail = FALSE)
}

# p.art is Augmented RTP function from Vsevolozhskaya et al (2019).
p.art <- function(K, p) {
    L <- length(p)
    lw <- sum(log(p[1:K - 1]))
    pk <- p[K]
    d <- (K - 1) * (digamma(L + 1) - digamma(K))
    q <- qgamma(1 - pbeta(pk, K, L - K + 1), shape = d)
    ak <- (K - 1) * log(pk) - lw + q
    1 - pgamma(ak, shape = K + d - 1)
}

# p.rtp is a robust and fast function for practically
# accurate Rank Truncated Product p-value.
p.rtp <- function(K, p, tol = 1e-10, stepscale = 1) {
    L <- length(p)
    if (K > L) {
        return(-1)
    }
    if (K == L) {
        return(p.fisher(p))
    }
    if (K == 1) {
        return(1 - (1 - p[1])^L)
    }
    riemannGamma(K, p, tol, stepscale)
}

# fBetaQuantile.R is integrand over Beta probabilities in [0, 1].
# This is the integrand in Vsevolozhskaya et al. (2019).
fBetaQ.R <- function(p, lw, K, L) {
    b <- qbeta(p, K + 1, L - K)
    g <- log(b) * K - lw
    1 - pgamma(g, K)
    # pgamma(g, K, lower.tail = FALSE)
}

# fGammaQ.R is integrand over Gamma probabilities in [0, 1].
# fGammaQ(1 - u) is very close to fBetaQ(u).
fGammaQ.R <- function(p, lw, K, L) {
    g <- qgamma(p, K)
    b <- exp((g + lw) / K)
    pbeta(b, K + 1, L - K)
}

# fBetaD.R is integrand over Beta density in [0, 1]. This is
# equivalent to the integrand equation in Dudbridge and Koeleman (2003).
fBetaD.R <- function(b, lw, K, L) {
    g <- K * log(b) - lw
    dbeta(b, K + 1, L - K) * pgamma(g, K, lower.tail = FALSE)
}

# fGammaD.R is integrand over Gamma density in [0, inf).
fGammaD.R <- function(g, lw, K, L) {
    b <- exp((g + lw) / K)
    dgamma(g, K) * pbeta(b, K + 1, L - K)
}

# fBetaD ----------------------------------------- fBetaD

# Reference integration by library cubature function pcubature.
p.rtp.dbeta.cuba <- function(K, p, tol = 1e-15) {
    if (init(K, p) < 0) {
        return(-1)
    }
    top <- fBetaDtop()

    pcubature(fBetaD, 0, top, tol = tol)$integral +
        pcubature(fBetaD, top, 1, tol = tol)$integral
}

# RPT p-value by Beta density and adaptive Simpson's 1/3.
p.rtp.dbeta.simp.a <- function(K, p, abstol = 1e-7, reltol = 1e-3) {
    simpsonAdaBeta(K, p, abstol, reltol)
}

# RPT p-value by Beta density and Riemann sum integration.
p.rtp.dbeta.riema <- function(K, p, tol = 1e-10, stepscale = 1) {
    riemannBeta(K, p, tol, stepscale)
}

# RPT p-value by Beta density and R integrate function.
p.rtp.dbeta.integrate <- function(K, p, abstol = 1e-3, reltol = 1e-2) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    f <- function(u) fBetaD.R(u, lw, K, L)

    # This may fail for small K/L.
    init(K, p)
    top <- fBetaDtop()
    integrate(f, 0, top, abs.tol = abstol, rel.tol = reltol)$value +
        integrate(f, top, 1, abs.tol = abstol, rel.tol = reltol)$value
}

#  fGammaD --------------------------------------- fGammaD

# RPT p-value by Gamma density and fixed step Simpson's 1/3.
p.rtp.dgamma.simp <- function(K, p, tol = 1e-10, stepscale = 1) {
    simpsonGamma(K, p, tol, stepscale)
}

# RPT p-value by Gamma density and Riemann sum integration.
p.rtp.dgamma.riema <- function(K, p, tol = 1e-10, stepscale = 1) {
    riemannGamma(K, p, tol, stepscale)
}

# RPT p-value by Gamma density and R integrate function.
p.rtp.dgamma.integrate <- function(K, p, abstol = 1e-2, reltol = 1e-2) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    f <- function(u) fGammaD.R(u, lw, K, L)
    end <- K + 25 * sqrt(K) # mean + 25 x SD
    # integrate(f, 0, Inf, ..)  is slower and less accurate.

    integrate(f, 0, end, abs.tol = abstol, rel.tol = reltol)$value
}

# fBetaQ ----------------------------------------- fBetaQ

# RPT p-value by inverse beta CDF and adaptive Simpson's 1/3.
p.rtp.qbeta.simp.a <- function(K, p, abstol = 1e-7, reltol = 1e-3) {
    simpsonAdaBetaQ(K, p, abstol, reltol, depth = 25)
}

# RPT p-value by inverse Beta CDF and R integrate function.
# Inverse beta CDF/Quantile function method from Vsevolozhskaya et al (2019).
p.rtp.qbeta.integrate <- function(K, p, abstol = 1e-4, reltol = 1e-2) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    f <- function(u) fBetaQ.R(u, lw, K, L)

    integrate(f, 0, 1, abs.tol = abstol, rel.tol = reltol)$value
}

# fGammaQ ---------------------------------------- fGammaQ

# RPT p-value by inverse Gamma CDF and adaptive Simpson's 1/3.
p.rtp.qgamma.simp.a <- function(K, p, abstol = 1e-7, reltol = 1e-3) {
    simpsonAdaGammaQ(K, p, abstol, reltol, depth = 25)
}

# RPT p-value by inverse Gamma CDF and R integrate function.
p.rtp.qgamma.integrate <- function(K, p, abstol = 1e-4, reltol = 1e-2) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    f <- function(u) fGammaQ.R(u, lw, K, L)

    integrate(f, 0, 1, abs.tol = abstol, rel.tol = reltol)$value
}