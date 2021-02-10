
# p <- p.gen(L = 100, small = 1e-4, seed = 0)
p.gen <- function(L, small, seed = 0) {
    if (seed > 0) set.seed(seed)
    sort(c(small, runif(L - 1)))
}

# pvalues(K=10, L=100, small=1e-6, seed=0)
pvalues <- function(K, L, small, seed) {
    if (small < 1e-300) {
        return("Invalid: small < 1e-300")
    }
    if (K >= L) {
        return("Invalid: K >= L")
    }
    if (K > 5000) {
        return("Invalid: K > 5000")
    }
    p <- p.gen(L, small, seed)
    p1 <- p.rtp.qbeta.integrate(K, p)
    p2 <- p.rtp.qbeta.simp.a(K, p)
    p12 <- p.rtp.qgamma.simp.a(K, p)
    p7 <- p.rtp.dbeta.integrate(K, p)
    p8 <- p.rtp.dbeta.simp.a(K, p)
    p4 <- p.rtp.dbeta.riema(K, p, stepscale = 1)
    p3 <- p.rpt.dbeta.cuba(K, p, tol = 1e-6)
    p10 <- p.art(K, p)
    p9 <- p.tfisher.soft((K + 1) / (L + 1), p)
    p6 <- ranktruncated(K, p)
    p11 <- p.rtp.dgamma.riema(K, p, stepscale = 1)
    p13 <- p.rtp.dgamma.simp(K, p, stepscale = 1)
    pe <- p.rpt.dbeta.cuba(K, p, tol = 1e-14) # "exact" reference

    w <- function(s) writeLines(s)
    fmt <- function(p) {
        format(p,
            width = wi,
            digits = d, nsmall = nsmall, drop0trailing = droptrail,
            justify = c("left")
        )
    }
    ferr <- function(p) {
        format(-log10(abs(p - pe) / pe), digits = d)
    }
    d <- 6
    wi <- 10
    nsmall <- 8
    droptrail <- FALSE

    f1 <- fmt(p1)
    f2 <- fmt(p2)
    f3 <- fmt(p3)
    f4 <- fmt(p4)
    f6 <- fmt(p6)
    f7 <- fmt(p7)
    f8 <- fmt(p8)
    f9 <- fmt(p9)
    f10 <- fmt(p10)
    f11 <- fmt(p11)
    f12 <- fmt(p12)
    f13 <- fmt(p13)
    d <- 1
    e1 <- format(abs(p1 - pe), digits = d)
    e2 <- format(abs(p2 - pe), digits = d)
    e3 <- format(abs(p3 - pe), digits = d)
    e4 <- format(abs(p4 - pe), digits = d)
    e6 <- format(abs(p6 - pe), digits = d)
    e7 <- format(abs(p7 - pe), digits = d)
    e8 <- format(abs(p8 - pe), digits = d)
    e9 <- format(abs(p9 - pe), digits = d)
    e10 <- format(abs(p10 - pe), digits = d)
    e11 <- format(abs(p11 - pe), digits = d)
    e12 <- format(abs(p12 - pe), digits = d)
    e13 <- format(abs(p13 - pe), digits = d)
    d <- 2
    d1 <- ferr(p1)
    d2 <- ferr(p2)
    d3 <- ferr(p3)
    d4 <- ferr(p4)
    d6 <- ferr(p6)
    d7 <- ferr(p7)
    d8 <- ferr(p8)
    d9 <- ferr(p9)
    d10 <- ferr(p10)
    d11 <- ferr(p11)
    d12 <- ferr(p12)
    d13 <- ferr(p13)
    tab <- "   "
    w("")
    w("                                   Abs     # of correct")
    w("----Function---------P-value-------error---nonzero digits")
    w(paste("p.rtp.qbeta.integ  ", f1, tab, e1, "  ", d1))
    w(paste("p.rtp.qbeta.simp.a ", f2, tab, e2, "  ", d2))
    w(paste("p.rtp.qgamm.simp.a ", f12, tab, e12, "  ", d12))
    w(paste("p.rtp.dbeta.integ  ", f7, tab, e7, "  ", d7))
    w(paste("p.rpt.dbeta.cuba   ", f3, tab, e3, "  ", d3))
    # w(paste("mutoss/ranktrunca  ", f6, tab, e6, "  ", d6))
    w("")
    w(paste("p.rtp.dbeta.simp.a ", f8, tab, e8, "  ", d8))
    w(paste("p.rtp.dbeta.riema  ", f4, tab, e4, "  ", d4))
    w(paste("p.rtp.dgamma.simp  ", f13, tab, e13, "  ", d13))
    w(paste("p.rtp.dgamma.riema ", f11, tab, e11, "  ", d11))
    # w(paste("p.tfisher.soft     ", f9, tab, e9, tab, "  ", d9))
    # w(paste("p.art              ", f10, tab, e10, tab, "  ", d10))
}

# plotQuantile(K=10, L=100, small=1e-4, seed=0)
plotQuantile <- function(K, L, small, seed, xmax = 1) {
    p <- p.gen(L, small, seed)
    lw <- sum(log(p[1:K]))
    pval <- p.rtp.qbeta.integrate(K, p)

    f1 <- function(u) fBetaQuantile.R(u, lw, K, L)
    f2 <- function(u) fGammaQuantile.R(u, lw, K, L)

    plot(f1,
        type = "l", ylab = "", font.main = 1,
        col = "blue", cex.main = 1,
        xlim = c(0, xmax), ylim = c(0, 1), yaxp = c(0, 1, 2),

        axis(1,
            at = pretty(c(0, xmax)),
        ),
        xlab = paste(
            "Area under the integrands is ",
            "p = ", format(pval, digits = 2)
        ),
        main = c(
            "Integrand by inverse CDF (quantile function)",
            "by Beta inverse (blue) and Gamma inverse (red)"
        )
    )
    par(new = TRUE)
    plot(f2,
        type = "l", ylab = "", col = "red",
        xlim = c(0, xmax), ylim = c(0, 1),
        xlab = "",
        axes = FALSE,
    )

    abline(v = 0, lty = 1, col = "#757474", lwd = 0.25)
    abline(v = 1, lty = 1, col = "#757474", lwd = 0.25)
    abline(v = 0.5, lty = 1, col = "#757474", lwd = 0.25)
    abline(h = 0, lty = 1, col = "#757474", lwd = 0.25)
    abline(h = 0.5, lty = 1, col = "#757474", lwd = 0.25)
}

# plotBxGintegrand(K=10, L=100, small=1e-4, seed=0)
plotBxGintegrand <- function(K, L, small, seed, xmin = -1, xmax = 0) {
    p <- p.gen(L, small, seed)
    pval <- p.rpt.dbeta.cuba(K, p)
    init(K, p)
    top <- fBetaTop()
    lw <- sum(log(p[1:K]))
    left <- exp(lw / K) * 1.5

    right <- K / (L - 1)
    hight <- dbeta(K / (L - 1), K + 1, L - K)

    sd <- betaSD(K + 1, L - K)
    if (xmax <= 0) xmax <- min(1.0, top + 8 * sd)
    if (xmin < 0) xmin <- max(0, top - 6 * sd)

    f1 <- function(u) 1 - pgamma(K * log(u) - lw, K)
    f2 <- function(u) fBeta.R(u, lw, K, L)
    f3 <- function(u) dbeta(u, K + 1, L - K)

    plot(f3,
        type = "l", lwd = 1, font.main = 1, cex.main = 1, col = "darkgreen",
        xlim = c(xmin, xmax), ylim = c(0, hight),
        ylab = "", yaxt = "n", xaxt = "n",
        xlab = c(
            paste(
                "Area under the blue integrand is p = ",
                format(pval, digits = 2)
            ),
            "Area under the green Beta density is 1"
        ),
        main = c(
            "BxG integrand (blue) = Beta density (green) x Gamma survival (red)",
            paste(
                " K = ", format(K, digits = 5),
                " L = ", format(L, digits = 5)
            ),
            paste(
                " Approx. integrand location = ",
                format(top, digits = 2)
            )
        )
    )
    axis(1,
        col.ticks = "#3a3939", col.axis = "#3a3939",
        at = pretty(c(xmin, xmax))
    )
    axis(2,
        col.ticks = "darkgreen", col.axis = "darkgreen",
    )
    abline(v = left, lty = 2, col = "red", lwd = 0.5)
    abline(v = K / (L - 1), lty = 2, col = "darkgreen", lwd = 0.5)
    abline(v = top, lty = 2, col = "blue", lwd = 0.5)
    abline(h = 0, lty = 1, col = "#757474", lwd = 0.25)
    abline(v = 0, lty = 1, col = "#757474", lwd = 0.25)

    par(new = TRUE)
    plot(f2,
        type = "l", lwd = 1, col = "blue",
        xlim = c(xmin, xmax), ylim = c(0, hight),
        ylab = "", xlab = "", axes = FALSE,
    )
    par(new = TRUE)
    plot(f1,
        type = "l", lwd = 1, col = "red",
        xlim = c(xmin, xmax), ylim = c(0, 1),
        ylab = "", xlab = "", axes = FALSE,
    )
    axis(4,
        col.ticks = "red", col.axis = "red",
        at = c(0, 0.5, 1)
    )
}

# plotIntegrandLocation(K=10, L=100, small=1e-4, seed=0)
plotIntegrandLocation <- function(K, L, small, seed, xmin = -1, xmax = 0) {
    p <- p.gen(L, small, seed)
    lw <- sum(log(p[1:K]))
    init(K, p)
    pval <- p.rpt.dbeta.cuba(K, p)
    iTop <- fBetaTop()
    bTop <- K / (L - 1)
    bSD <- betaSD(K + 1, L - K)
    gTop <- exp((K - 1 + lw) / K)
    gSD <- exp((sqrt(K) + lw) / K)

    if (xmax <= 0) xmax <- min(1.0, iTop + 6 * bSD)
    if (xmin < 0) xmin <- max(0.0, iTop - 6 * bSD)

    f1 <- function(u) 1 - pgamma(K * log(u) - lw, K)
    f2 <- function(u) fBeta.R(u, lw, K, L)
    f3 <- function(u) dbeta(u, K + 1, L - K)

    plot(f1,
        type = "l", font.main = 1, lwd = 1, cex.main = 1, col = "red",
        xlim = c(xmin, xmax), # yaxp = c(0, 1, 2),
        axis(1,
            at = c(xmin, xmax),
        ),
        xlab = "Blue and green curves are vertically not in actual size",
        ylab = "",
        main = c(
            paste(
                "Gamma survival (red), ", "Beta density (green), ",
                "Beta x Gamma integrand (blue)."
            ),
            paste(
                "K = ", format(K, digits = 5),
                " L = ", format(L, digits = 5),
                " p = ", format(pval, digits = 2)
            )
        )
    )
    abline(v = iTop, col = "blue", lty = 2, lwd = 0.5)
    abline(v = bTop, col = "#088308", lty = 2, lwd = 0.5)
    abline(v = gTop, col = "red", lty = 2, lwd = 0.5)
    abline(h = 0, lty = 1, col = "gray", lwd = 1)
    par(new = TRUE)
    plot(f3,
        type = "l", lwd = 1, xlim = c(xmin, xmax),
        yaxt = "n", xaxt = "n", xlab = "", ylab = "",
        col = "#088308",
    )
    par(new = TRUE)
    plot(f2,
        type = "l", lwd = 1, xlim = c(xmin, xmax),
        yaxt = "n", xaxt = "n", xlab = "", ylab = "",
        col = "blue",
    )
}

# plotGxBintegrand(K=10, L=100, small=1e-4, seed=0)
plotGxBintegrand <- function(K, L, small = 1e-1, seed = 0, xmin = -1, xmax = 0) {
    p <- p.gen(L, small, seed)
    pval <- p.rpt.dbeta.cuba(K, p, 1e-10)
    lw <- sum(log(p[1:K]))
    init(K, p)
    hight <- dgamma(K - 1, K)

    f3 <- function(g) dgamma(g, K)
    f2 <- function(g) fGamma.R(g, lw, K, L)
    f1 <- function(g) {
        b <- exp((g + lw) / K)
        pbeta(b, K + 1, L - K)
    }

    bTop <- K * log(K / (L - 1)) - lw
    if (xmax <= 0) xmax <- floor(bTop + 4 * sqrt(K))
    if (xmin < 0) xmin <- max(0, floor(K - 1 - 3 * sqrt(K)))
    gTop <- K - 1
    iTop <- (gTop + 2 * bTop) / 3

    plot(f2,
        type = "l", lwd = 1, font.main = 1, cex.main = 1, col = "blue",
        xlim = c(xmin, xmax), ylim = c(0, hight),
        ylab = "", yaxt = "n", xaxt = "n",
        xlab = c(
            paste(
                "Area under the blue integrand is p = ",
                format(pval, digits = 2)
            ),
            "Area under the green Gamma density is 1"
        ),
        main = c(
            "GxB integrand (blue) = Gamma density (green) x Beta CDF (red)",
            paste(
                " K = ", format(K, digits = 5),
                " L = ", format(L, digits = 5)
            )
        )
    )
    axis(1,
        col.ticks = "#3a3939", col.axis = "#3a3939",
        at = pretty(c(xmin, xmax))
    )
    axis(2,
        at = pretty(c(0, hight / 2, hight)),
        col.ticks = "blue", col.axis = "blue",
    )
    abline(v = gTop, lty = 2, col = "darkgreen", lwd = 0.5)
    abline(v = iTop, lty = 2, col = "blue", lwd = 0.5)
    abline(v = bTop, lty = 2, col = "red", lwd = 0.5)
    abline(h = 0, lty = 1, col = "#757474", lwd = 0.25)
    # abline(v = 0, lty = 1, col = "#757474", lwd = 0.25)

    par(new = TRUE)
    plot(f1,
        type = "l", lwd = 1, col = "red",
        xlim = c(xmin, xmax), ylim = c(0, 1),
        ylab = "", xlab = "", axes = FALSE,
    )
    axis(4,
        at = c(0, 0.5, 1),
        col.ticks = "red", col.axis = "red",
    )
    par(new = TRUE)
    plot(f3,
        type = "l", lwd = 1, col = "darkgreen",
        xlim = c(xmin, xmax), ylim = c(0, hight),
        xlab = "", ylab = "", axes = FALSE,
    )
}

# bench(K=10, L=100, small=1e-5, seed=0)
bench <- function(K, L, small, seed) {
    library(microbenchmark)
    p <- p.gen(L, small, seed)
    pval <- p.rpt.dbeta.cuba(K, p, 1e-10)

    QBinte <- function() p.rtp.qbeta.integrate(K, p)
    QGinte <- function() p.rtp.qgamma.integrate(K, p)
    DBinte <- function() p.rtp.dbeta.integrate(K, p)

    QBsimA <- function() p.rtp.qbeta.simp.a(K, p)
    QGsimA <- function() p.rtp.qgamma.simp.a(K, p)

    DBsimA <- function() simpsonAdaBeta(K, p)
    #  DBsimA <- function() p.rtp.dbeta.simp.a(K, p)

    DBriem <- function() riemannBeta(K, p)
    # DBriem <- function() p.rtp.dbeta.riema(K, p)

    DGriem <- function() riemannGamma(K, p, stepscale = 1)
    # DGriema <- function() p.rtp.dgamma.riema(K, p, stepscale = 1)

    # DGsimp <- function() p.rtp.dgamma.simp(K, p, stepscale = 1)
    DGsimp <- function() simpsonGamma(K, p, stepscale = 1)

    DBcuba <- function() p.rpt.dbeta.cuba(K, p)
    # TFish <- function() p.tfisher.soft.R(0.05, p)
    TFish <- function() p.tfisher.soft(0.05, p)
    DBmutos <- function() ranktruncated(K, p)
    Art <- function() p.art(K, p)

    res <- microbenchmark(
        unit = "us",
        # DBmutos(),
        # QBinte(),
        QBsimA(),
        QGsimA(),
        # DBinte(),
        DBsimA(),
        DBriem(),
        DGsimp(),
        DGriem(),
        # DBcuba(),
        # TFish(),
        # Art(),
        times = 500
    )
    # plotQuantile(K, L, small, seed)
    boxplot(res,
        unit = "us",
        font.main = 1, cex.main = 1,
        xlab = "XXyyyy = Integrand function and integrating function",
        font.xlab = 1,
        # names = c(
        #     "Qinteg", "Qsimpa", "Dinteg", "Dcuba",
        #     "Dsimpa", "Driema"
        # ),
        main = paste(
            "Rank truncated p-value integrating time\n",
            "microseconds, logarithmic time",
            " p = ", format(pval, digits = 2)
        )
    )
    print(res, "us", signif = 3)
}

# benchIntegrands(K=1, L=100, small=1e-5, seed=0, unit="us")
# R adds ~900 ns baseline cost (fNull) for these C-functions.
benchIntegrands <- function(K, L, small, seed, unit = "us") {
    library(microbenchmark)
    p <- p.gen(L, small, seed)
    lw <- sum(log(p[1:K]))
    # initCpp(lw, K, L)
    init(K, p)
    # a <- 10
    # b <- 100
    # u <- runif(1)
    # g <- K * log(u) - lw

    Qgamma <- function(u) fGammaQuantile(u)
    Qbeta <- function(u) fBetaQuantile(u)
    QgammaR <- function(u) fGammaQuantile.R(u, lw, K, L)
    QbetaR <- function(u) fBetaQuantile.R(u, lw, K, L)
    pgamm <- function(g) pgamma(g, K)
    fBetaR <- function(u) fBeta.R(u, lw, K, L)

    title <- "Microseconds, logarithmic time"
    if (unit == "ns") title <- "nanoseconds, logarithmic time"
    set.seed(1)

    res <- microbenchmark(
        setup = {
            u <- runif(1)
            g <- K * log(u) - lw
        },
        Qbeta(u),
        Qgamma(u),
        fGamma(g),
        # fBetaR(u),
        fBeta(u),
        fNull(u),

        # QgammaR(u),
        # QbetaR(u),
        # pgamm(g),
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