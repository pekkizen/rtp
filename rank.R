

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

# fBetaQuantile is integrand over Beta probabilities in [0, 1].
# This is the integrand in Vsevolozhskaya et al. (2019).
fBetaQuantile <- function(u, lw, K, L) {
    b <- qbeta(u, K + 1, L - K)
    g <- log(b) * K - lw
    1 - pgamma(g, K)
}

# fBetaDensity is integrand over Beta density in [0, 1]. This is
# equivalent to the integrand equation in Dudbridge and Koeleman (2003).
fBetaDensity <- function(u, lw, K, L) {
    g <- K * log(u) - lw
    dbeta(u, K + 1, L - K) * (1 - pgamma(g, K))
}

# fGammaDensity is integrand over Gamma density in [0, inf).
fGammaDensity <- function(g, lw, K, L) {
    u <- exp((g + lw) / K)
    dgamma(g, K) * pbeta(u, K + 1, L - K)
}

# Reference integration by library cubature function pcubature.
p.rpt.dbeta.cuba <- function(K, p, tol = 1e-2) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    initCpp(lw, K, L)
    top <- fBetaTop(lw, K, L)

    I <- pcubature(fBetaCpp, 0, top, tol = tol / 2)$integral
    I + pcubature(fBetaCpp, top, 1, tol = tol / 2)$integral
}

# RPT p-value by Beta density and adaptive Cpp Simpson's 1/3.
p.rtp.dbeta.simp.a <- function(K, p, abstol = 1e-7, reltol = 1e-3) {
    L <- length(p)
    lw <- sum(log(p[1:K]))

    simpsonAdaBetaCpp(lw, K, L, abstol, reltol, depth = 25)
}

# RPT p-value by Gamma density and fixed step Cpp Simpson's 1/3.
p.rtp.dgamma.simp <- function(K, p, tol = 1e-12, stepscale = 1) {
    L <- length(p)
    lw <- sum(log(p[1:K]))

    simpsonGammaCpp(lw, K, L, tol, stepscale)
}

# RPT p-value by Beta density and Riemann sum Cpp integration.
p.rtp.dbeta.riema <- function(K, p, tol = 1e-12, stepscale = 1) {
    L <- length(p)
    lw <- sum(log(p[1:K]))

    riemannBetaCpp(lw, K, L, tol, stepscale)
}

# RPT p-value by Gamma density and Riemann sum Cpp integration.
p.rtp.dgamma.riema <- function(K, p, tol = 1e-12, stepscale = 1) {
    L <- length(p)
    lw <- sum(log(p[1:K]))

    riemannGammaCpp(lw, K, L, tol, stepscale)
}

# RPT p-value by beta density and R integrate function.
p.rtp.dbeta.integrate <- function(K, p, abstol = 1e-4) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    f <- function(u) fBetaDensity(u, lw, K, L)

    # This fails too often, e.g K=5, L=1000:
    # the integral is probably divergent
    # integrate(f, 0, 1, abs.tol = abstol, rel.tol = 1e-2)$value

    # This also may fail for small K/L.
    top <- fBetaTop(lw, K, L)
    I <- integrate(f, 0, top, abs.tol = abstol, rel.tol = 1e-2)$value
    I + integrate(f, top, 1, abs.tol = abstol, rel.tol = 1e-2)$value
}

# RPT p-value by inverse beta CDF and R integrate function.
# Inverse betaCDF/Quantile function method from Vsevolozhskaya et al (2019).
p.rtp.qbeta.integrate <- function(K, p, abstol = 1e-5) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    f <- function(u) fBetaQuantile(u, lw, K, L)

    integrate(f, 0, 1, abs.tol = abstol, rel.tol = 1e-2)$value
}

# RPT p-value by inverse beta CDF and adaptive Cpp Simpson.
p.rtp.qbeta.simpa <- function(K, p, abstol = 1e-6, reltol = 1e-4) {
    L <- length(p)
    lw <- sum(log(p[1:K]))

    simpsonAdaBetaQuantileCpp(lw, K, L, abstol, reltol, depth = 25)
}