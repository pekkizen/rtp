

# bench.integrals(K=10, L=100, plot = FALSE, seed=0)
bench.integrals <- function(K, L, seed = 0, small = 1e-1, plot = FALSE) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    p <- c(small, runif(L - 1))
    pval <- p.rtp.dbeta.cuba(K, p)

    res <- microbenchmark::microbenchmark(
        setup = {
            p <- c(small, runif(L - 1))
        },
        unit = "us",
        DBmutos = p.rtp.mutoss(K, p),
        QBinte = p.rtp.qbeta(K, p),
        # DBsimA = pRtpDbetaAsimp(K, p),
        DBriem = pRtpDbetaRiema(K, p),
        DGsimp = pRtpDgammaSimp(K, p, stepscale = 1),
        DGriem = pRrtpDgammaRiema(K, p, stepscale = 1),

        # DBcuba = p.rtp.dbeta.cuba(K, p),
        # tau <- K / L
        # TFish = p.tfisher.softR(tau, p),
        # TFish = p.tfisher.soft(tau, p),
        # Art = p.art(K, p),
        times = 1000
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
# R adds ~900 ns baseline cost for Rcpp functions.
bench.integrands <- function(K, L, small = 1e-1, unit = "us", plot = FALSE) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    p <- c(small, runif(L - 1))
    init(K, p)
    lw <- sum(log(p[1:K]))

    title <- "Microseconds, logarithmic time"
    if (unit == "ns") title <- "nanoseconds, logarithmic time"

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
        times = 5000
    )
    print(res, unit, signif = 3)
    if (plot) {
        plot.new()
        boxplot(res,
            font.main = 1, cex.main = 1,
            xlab = "",
            font.xlab = 1,
            unit = unit,
            main = c("Integrands", title)
        )
    }
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