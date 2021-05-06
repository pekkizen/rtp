

# bench.pvalues(K=10, L=100)
bench.pvalues <- function(K, L, small = 1, tau = 0) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    if (small == 1) {
        small <- K / L * 0.01
    }
    p <- c(small, runif(L - 1))
    if (tau == 0) tau <- K / L

    res <- microbenchmark::microbenchmark(
        p.rtp.mutoss = p.rtp.mutoss(K, p),
        p.rtp.qbeta = p.rtp.qbeta(K, p),
        p.art = p.art(K, p),
        p.tfisher.softR = p.tfisher.softR(tau, p),
        p.tfisher.soft = p.tfisher.soft(tau, p),
        p.rtp = rtpDgammaRiema(K, p),

        # p.rtp.dbeta.asimp = rtpDbetaAsimp(K, p),
        # p.rtp.dbeta.riema = rtpDbetaRiema(K, p),
        # p.rtp.dgamma.simp = rtpDgammaSimp(K, p),
        #  p.rtp.dbeta.cuba = p.rtp.dbeta.cuba(K, p),
        times = 500
    )
    print(res, "us", signif = 3)
}

# bench.integrands(K=5, L=100)
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
        Qbeta = fBetaQ(u),
        QbetaR = fBetaQ.R(u, lw, K, L),
        QGamma = fGammaQ(u),
        QGammaR = fGammaQ.R(u, lw, K, L),
        DGamma = fGammaD(g),
        DGammaR = fGammaD.R(u, lw, K, L),
        DBeta = fBetaD(u),
        DBetaR = fBetaD.R(u, lw, K, L),
        fNull = baseNull(u),
        times = 10000
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
        uniSelect = uniSel(K, p),
        simpleSelect = simpleSel(K, p),
        nth_element = nth_element(K, p),
        sort_partial = sort(p, partial = c(1:K)),
        sortFull = sort(p),
        fNull = baseNull(1),
        times = times
    )
    print(res, unit = "us", signif = 3)
}