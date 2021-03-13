

# bench.integrals(K=10, L=100, plot = FALSE, seed=0)
bench.integrals <- function(K, L, seed = 0, small = 1e-1, plot = FALSE) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    p <- c(small, runif(L - 1))
    tau <- K / L

    res <- microbenchmark::microbenchmark(
        p.rtp.mutoss = p.rtp.mutoss(K, p),
        p.rtp.qbeta = p.rtp.qbeta(K, p),
        p.art = p.art(K, p),
        p.rtp = rtpDgammaRiema(K, p, stepscale = 1),
        p.tfisher.soft = p.tfisher.soft(tau, p),
        p.rtp.dbeta.asimp = rtpDbetaAsimp(K, p),
        p.rtp.dbeta.riema = rtpDbetaRiema(K, p),
        p.rtp.dgamma.simp = rtpDgammaSimp(K, p, stepscale = 1),
        #  p.rtp.dbeta.cuba = p.rtp.dbeta.cuba(K, p),
        times = 500
    )
    print(res, "us", signif = 3)
    if (plot) {
        plot.new()
        boxplot(res,
            unit = "us",
            font.main = 1, cex.main = 1,
            xlab = "XXyyyy = Integrand function and integrating function",
            font.xlab = 1,
            main = paste(
                "Rank truncated p-value integrating time\n",
                "microseconds, logarithmic time",
                " p = ", format(pval, digits = 2)
            )
        )
    }
}

# bench.integrands(K=5, L=100, plot = FALSE)
# R adds ~1000 ns baseline cost for Rcpp functions, fNull.
bench.integrands <- function(K, L, small = 1e-1, unit = "us") {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    p <- c(small, runif(L - 1))
    lw <- stat.rtp(K, p)
    init(K, p)

    res <- microbenchmark::microbenchmark(
        setup = {
            u <- runif(1)
            g <- K * log(u) - lw
        },
        Qbeta = fBetaQ(u),
        QbetaR = fBetaQ.R(u, lw, K, L),
        QGamma = fGammaQ(u),
        DGamma = fGammaD(g),
        DBeta = fBetaD(u),
        DBetaR = fBetaD.R(u, lw, K, L),
        fNull = baseNull(u),
        times = 5000
    )
    print(res, unit, signif = 3)
}

# bench.select(K=10, L=1000, times=2000)
bench.select <- function(K, L, unit = "us", times = 2000) {
    res <- microbenchmark::microbenchmark(
        setup = {
            p <- c(runif(L) * 1.0)
            # p <- sort(p, decreasing = TRUE)
            # p <- sort(p)
        },
        uniSelect = uniSel(K, p),
        simpleSelect = simpleSel(K, p),
        nth_element = nth_element(K, p),
        sort_partial = sort(p, partial = c(1:K)),
        # sortFull = sort(p),
        times = times
    )
    print(res, unit, signif = 4)
}