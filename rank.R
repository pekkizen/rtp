

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

# fBetaQuantile.R is integrand over Beta probabilities in [0, 1].
# This is the integrand in Vsevolozhskaya et al. (2019).
fBetaQuantile.R <- function(p, lw, K, L) {
    b <- qbeta(p, K + 1, L - K)
    g <- log(b) * K - lw
    1 - pgamma(g, K)
}

# fGammaQuantile.R is integrand over Gamma probabilities in [0, 1].
# fGammaQuantile(1 - u) is very close to fBetaQuantile(u).
fGammaQuantile.R <- function(p, lw, K, L) {
    g <- qgamma(p, K)
    b <- exp((g + lw) / K)
    pbeta(b, K + 1, L - K)
}

# fBeta.R is integrand over Beta density in [0, 1]. This is
# equivalent to the integrand equation in Dudbridge and Koeleman (2003).
fBeta.R <- function(b, lw, K, L) {
    g <- K * log(b) - lw
    dbeta(b, K + 1, L - K) * (1 - pgamma(g, K))
}

# fGamma.R is integrand over Gamma density in [0, inf).
fGamma.R <- function(g, lw, K, L) {
    b <- exp((g + lw) / K)
    dgamma(g, K) * pbeta(b, K + 1, L - K)
}

# Reference integration by library cubature function pcubature.
p.rpt.dbeta.cuba <- function(K, p, tol = 1e-2) {
    init(K, p)
    top <- fBetaTop()

    I <- pcubature(fBeta, 0, top, tol = tol / 2)$integral
    I + pcubature(fBeta, top, 1, tol = tol / 2)$integral
}

# RPT p-value by Beta density and adaptive Simpson's 1/3.
p.rtp.dbeta.simp.a <- function(K, p, abstol = 1e-7, reltol = 1e-3) {
    simpsonAdaBeta(K, p, abstol, reltol, depth = 25)
}

# RPT p-value by Gamma density and fixed step Simpson's 1/3.
p.rtp.dgamma.simp <- function(K, p, tol = 1e-8, stepscale = 1) {
    simpsonGamma(K, p, tol, stepscale)
}

# RPT p-value by Beta density and Riemann sum Cpp integration.
p.rtp.dbeta.riema <- function(K, p, tol = 1e-8, stepscale = 1) {
    riemannBeta(K, p, tol, stepscale)
}

# RPT p-value by Gamma density and Riemann sum Cpp integration.
p.rtp.dgamma.riema <- function(K, p, tol = 1e-8, stepscale = 1) {
    riemannGamma(K, p, tol, stepscale)
}

# RPT p-value by beta density and R integrate function.
p.rtp.dbeta.integrate <- function(K, p, abstol = 1e-5) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    f <- function(u) fBeta.R(u, lw, K, L)

    # This fails too often.
    # integrate(f, 0, 1, abs.tol = abstol, rel.tol = 1e-2)$value

    # This also may fail for small K/L.
    init(K, p)
    top <- fBetaTop()
    I <- integrate(f, 0, top, abs.tol = abstol, rel.tol = 1e-2)$value
    I + integrate(f, top, 1, abs.tol = abstol, rel.tol = 1e-2)$value
}

# RPT p-value by inverse beta CDF and R integrate function.
# Inverse betaCDF/Quantile function method from Vsevolozhskaya et al (2019).
p.rtp.qbeta.integrate <- function(K, p, abstol = 1e-5, reltol = 1e-2) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    f <- function(u) fBetaQuantile.R(u, lw, K, L)

    integrate(f, 0, 1, abs.tol = abstol, rel.tol = reltol)$value
}

# RPT p-value by inverse Gamma CDF and R integrate function.
p.rtp.qgamma.integrate <- function(K, p, abstol = 1e-5, reltol = 1e-2) {
    L <- length(p)
    lw <- sum(log(p[1:K]))
    f <- function(u) fGammaQuantile.R(u, lw, K, L)

    integrate(f, 0, 1, abs.tol = abstol, rel.tol = reltol)$value
}

# RPT p-value by inverse beta CDF and adaptive Cpp Simpson.
p.rtp.qbeta.simp.a <- function(K, p, abstol = 1e-6, reltol = 1e-3) {
    simpsonAdaBetaQuantile(K, p, abstol, reltol, depth = 25)
}

# RPT p-value by inverse Gamma CDF and adaptive Cpp Simpson.
p.rtp.qgamma.simp.a <- function(K, p, abstol = 1e-6, reltol = 1e-3) {
    simpsonAdaGammaQuantile(K, p, abstol, reltol, depth = 25)
}