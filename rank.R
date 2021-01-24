

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

# fQuantileCpp is integrand by Beta quantile function qbeta.
# This is the integrand in Vsevolozhskaya et al. (2019).
fQuantile <- function(u, lw, K, L) {
    b <- qbeta(u, K + 1, L - K)
    g <- log(b) * K - lw
    1 - pgamma(g, K)
}

# fDensity is integrand by Beta density x Gamma survival. This is
# equivalent to the integrand equation in Dudbridge and Koeleman (2003).
fDensity <- function(u, lw, K, L) {
    g <- K * log(u) - lw
    dbeta(u, K + 1, L - K) * (1 - pgamma(g, K))
}

# Reference integration by library cubature function pcubature.
p.rpt.dbeta.cuba <- function(K, p, tol = 1e-2) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    initCpp(lw, K, L)
    top <- fDenTop(lw, K, L)

    I <- pcubature(fDensityCpp, 0, top, tol = tol / 2)$integral
    I + pcubature(fDensityCpp, top, 1, tol = tol / 2)$integral
}

# RPT p-value by Beta density and adaptive Cpp Simpson
p.rtp.dbeta.simpa <- function(K, p, abstol = 1e-6, reltol = 1e-3) {
    L <- length(p)
    lw <- sum(log(p[1:K]))

    adaSimpDensityCpp(lw, K, L, abstol, reltol, depth = 25)
}

# RPT p-value by Beta density and Riemann sum Cpp integration.
p.rtp.dbeta.riema <- function(K, p, tol = 1e-12, steps = 6) {
    L <- length(p)
    lw <- sum(log(p[1:K]))

    riemannCpp(lw, K, L, tol, steps)
}

# RPT p-value by beta density and R integrate function.
p.rtp.dbeta.integrate <- function(K, p, abstol = 1e-4) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    f <- function(u) fDensity(u, lw, K, L)

    # This fails too often, e.g K=5, L=1000:
    # the integral is probably divergent
    # integrate(f, 0, 1, abs.tol = abstol, rel.tol = 1e-2)$value

    # This also may fail for small K/L.
    top <- fDenTop(lw, K, L)
    I <- integrate(f, 0, top, abs.tol = abstol, rel.tol = 1e-2)$value
    I + integrate(f, top, 1, abs.tol = abstol, rel.tol = 1e-2)$value
}

# RPT p-value by inverse betaCDF and R integrate function.
# Inverse betaCDF/Quantile function method from Vsevolozhskaya et al (2019).
p.rtp.qbeta.integrate <- function(K, p, abstol = 1e-5) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    f <- function(u) fQuantile(u, lw, K, L)

    integrate(f, 0, 1, abs.tol = abstol, rel.tol = 1e-2)$value
}

# RPT p-value by beta inverse CDF and adaptive Cpp Simpson.
p.rtp.qbeta.simpa <- function(K, p, abstol = 1e-6, reltol = 1e-4) {
    L <- length(p)
    lw <- sum(log(p[1:K]))

    adaSimpQuantileCpp(lw, K, L, abstol, reltol, depth = 25)
}