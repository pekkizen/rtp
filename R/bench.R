

# bench(K=10, L=100, small=1e-5, seed=0, tau=0.05)
benchIntegrals <- function(K, L, small, tau = 0.05) {
    err <- check(K, L, small)
    if (err != "") {
        return(err)
    }
    # p <- p.gen(K, L, small, seed)
    p <- c(small, runif(L - 1))
    q <- sort(p)
    lw <- stat.tfisher(p, tau, tau)
    pval <- p.rtp.dbeta.cuba(K, p, 1e-10)

    res <- microbenchmark(
        setup = {
            p <- c(small, runif(L - 1))
            # p <- sort(c(small, runif(L - 1)))
        },
        unit = "us",
        # DBmutos = function() ranktruncated(K, p),
        QBinte = p.rtp.qbeta.integrate(K, p),
        DBsimA = simpsonAdaBeta(K, p),
        DBriem = riemannBeta(K, p),
        DGsimp = simpsonGamma(K, p, stepscale = 1),
        DGriem = riemannGamma(K, p, stepscale = 1),

        # DBcuba = p.rtp.dbeta.cuba(K, p),
        # TFish = p.tfisher.softR(tau, p),
        # TFish = p.tfisher.soft(tau, p),
        # TFish = tfisher(lw, L, tau, tau),
        # Art = p.art(K, p),
        times = 500
    )
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
    print(res, "us", signif = 3)
}

# benchIntegrands(K=5, L=100, small=1e-3)
# R adds ~900 ns baseline cost (fNull) for these C-functions.
benchIntegrands <- function(K, L, small, unit = "us") {
    err <- check(K, L, small)
    if (err != "") {
        return(err)
    }
    p <- p.gen(K, L, small, 0)
    init(K, p)
    lw <- sum(log(p[1:K]))

    title <- "Microseconds, logarithmic time"
    if (unit == "ns") title <- "nanoseconds, logarithmic time"

    res <- microbenchmark(
        setup = {
            u <- runif(1)
            g <- K * log(u) - lw
        },
        Qbeta = fBetaQ(u),
        QbetaR = fBetaQ.R(u, lw, K, L),
        DGamma = fGammaD(g),
        DBeta = fBetaD(u),
        DBetaR = fBetaD.R(u, lw, K, L),
        fNull = fNull(u),
        times = 5000
    )
    boxplot(res,
        font.main = 1, cex.main = 1,
        xlab = "",
        font.xlab = 1,
        unit = unit,
        main = c("Integrands", title)
    )
    print(res, unit, signif = 3)
}

# benchSelect(K=25, L=1000, times=2000)
benchSelect <- function(K, L, small = 1e-1, unit = "us", times = 2000) {
    res <- microbenchmark(
        setup = {
            p <- c(runif(L) * 1.0)
            # p <- sort(p, decreasing = TRUE)
            # p <- sort(p)
        },
        uniSel = uniSel(K, p),
        simpleSel = simpleSel(K, p),
        nth_element = nth_element(K, p),
        sortPart = sort(p, partial = c(1:K)),
        sortFull = sort(p),
        times = times
    )
    print(res, unit, signif = 4)
}