

# bench.pvalues(K=10, L=100, times=1000)
bench.pvalues <- function(K, L, small = 1, tau = 0, times = 200) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    if (small == 1) {
        small <- 1 / L * 1e-3
    }
    p <- c(small, runif(L - 1))
    if (tau == 0) tau <- K / L
    # f <- function(u) drtp(u, K, L, tol = 1e-12, stepscale = 1)

    res <- microbenchmark::microbenchmark(
        p.rtp.mutoss = p.rtp.mutoss(K, p),
        # p.rtp.simulated = p.rtp.simulated(K, p, rounds = 20000),
        p.rtp.qbeta = p.rtp.qbeta(K, p),
        p.art = p.art(K, p),
        p.tfisher.softR = p.tfisher.softR(tau, p),
        p.tfisher.soft = p.tfisher.soft(tau, p),
        p.rtp.garpR = p.rtp.garpR(K, p),
        p.rtp.garp = p.rtp.garp(K, p),
        p.rtp.dgamma = p.rtp.dgamma(K, p),
        p.rtp.dbeta = p.rtp.dbeta(K, p),
        p.rtp.dbeta.asimp = p.rtp.dbeta.asimp(K, p, tol = 1e-3),
        p.rtp.dbeta.cuba = p.rtp.dbeta.cuba(K, p),
        # p.fisher = p.fisher(p),
        # p.fisherR = p.fisherR(p),

        fNull = baseNull(1),
        times = times
    )
    print(res, "us", signif = 3)
}

# bench.integrands(K=10, L=200)
# R adds ~1000 ns baseline cost for Rcpp functions, fNull.
bench.integrands <- function(K, L, small = 1e-1) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    p <- c(small, runif(L - 1))
    lw <- stat.rtp(K, p)
    init(K, p, 0)
    init(K, p, 1)

    res <- microbenchmark::microbenchmark(
        setup = {
            u <- runif(1)
            g <- K * log(u) - lw
        },
        QbetaR = fBetaQ.R(u, lw, K, L),
        QGammaR = fGammaQ.R(u, lw, K, L),
        DBetaR = fBetaD.R(u, lw, K, L),
        DGammaR = fGammaD.R(g, lw, K, L),
        DBeta = fBetaD(u),
        DGamma = fGammaD(g),
        fNull = baseNull(u),
        times = 5000
    )
    print(res, unit = "us", signif = 3)
}

# bench.select(K=10, L=1000, times=2000)
bench.select <- function(K, L, times = 2000) {
    res <- microbenchmark::microbenchmark(
        setup = {
            p <- c(runif(L))
            # p <- sort(p, decreasing = TRUE)
            # p <- sort(p)
        },
        kSelect = uniSel(K, p),
        simpleSelect = simpleSel(K, p),
        nth_element = nth_element(K, p),
        sort_partial = sort(p, partial = c(1:K)),
        sort_full = sort(p),
        fNull = baseNull(1),
        times = times
    )
    print(res, unit = "us", signif = 3)
}