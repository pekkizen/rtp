
# p <- p.gen(L = 100, small = 1e-4, seed = 0)
p.gen <- function(L, small, seed = 0) {
    if (seed > 0) set.seed(seed)
    sort(c(small, runif(L - 1)))
}

# pvalues(K=5, L=100, small=1e-6, seed=0)
pvalues <- function(K, L, small, seed) {
    # tol <- defTol
    p <- p.gen(L, small, seed)
    p1 <- p.rtp.qbeta.integrate(K, p)
    p2 <- p.rtp.qbeta.simp.a(K, p)
    p7 <- p.rtp.dbeta.integrate(K, p)
    p8 <- p.rtp.dbeta.simp.a(K, p)
    p4 <- p.rtp.dbeta.riema(K, p, stepscale = 1)
    p3 <- p.rpt.dbeta.cuba(K, p)
    p10 <- p.art(K, p)
    p9 <- p.tfisher.soft(K / L, p)
    p6 <- ranktruncated(K, p)
    p11 <- p.rtp.dgamma.riema(K, p, stepscale = 1)
    p13 <- p.rtp.dgamma.simp(K, p, stepscale = 1)

    pe <- p.rpt.dbeta.cuba(K, p, tol = 1e-14) # "exact" reference

    d <- 7
    f1 <- format(p1, digits = d)
    f2 <- format(p2, digits = d)
    f3 <- format(p3, digits = d)
    f4 <- format(p4, digits = d)
    f5 <- format(pe, digits = d)
    f6 <- format(p6, digits = d)
    f7 <- format(p7, digits = d)
    f8 <- format(p8, digits = d)
    f9 <- format(p9, digits = d)
    f10 <- format(p10, digits = d)
    f11 <- format(p11, digits = d)
    f13 <- format(p13, digits = d)
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
    e13 <- format(abs(p13 - pe), digits = d)
    d <- 2
    d1 <- format(-log10(abs(p1 - pe) / pe), digits = d)
    d2 <- format(-log10(abs(p2 - pe) / pe), digits = d)
    d3 <- format(-log10(abs(p3 - pe) / pe), digits = d)
    d4 <- format(-log10(abs(p4 - pe) / pe), digits = d)
    d6 <- format(-log10(abs(p6 - pe) / pe), digits = d)
    d7 <- format(-log10(abs(p7 - pe) / pe), digits = d)
    d8 <- format(-log10(abs(p8 - pe) / pe), digits = d)
    d9 <- format(-log10(abs(p9 - pe) / pe), digits = d)
    d10 <- format(-log10(abs(p10 - pe) / pe), digits = d)
    d11 <- format(-log10(abs(p11 - pe) / pe), digits = d)
    d13 <- format(-log10(abs(p13 - pe) / pe), digits = d)
    tab <- "  \t"
    writeLines("                                                   # of correct")
    writeLines("----Function-----------P-value--------Abs error----non zero digits")
    writeLines(paste("p.rtp.qbeta.integrate ", f1, tab, e1, tab, d1))
    writeLines(paste("p.rtp.qbeta.simp.a    ", f2, tab, e2, tab, d2))
    writeLines(paste("p.rtp.dbeta.integrate ", f7, tab, e7, tab, d7))
    writeLines(paste("p.rpt.dbeta.cuba      ", f3, tab, e3, tab, d3))

    writeLines(paste("p.rtp.dbeta.riema     ", f4, tab, e4, tab, d4))

    writeLines(paste("p.rtp.dbeta.simp.a    ", f8, tab, e8, tab, d8))
    writeLines(paste("p.rtp.dgamma.simp     ", f13, tab, e13, tab, d13))
    writeLines(paste("p.rtp.dgamma.riema    ", f11, tab, e11, tab, d11))

    # writeLines(paste("mutoss/ranktruncated  ", f6, tab, e6, tab, d6))
    # writeLines(paste("p.tfisher.soft        ", f9, tab, e9, tab, d9))
    # writeLines(paste("p.art                 ", f10, tab, e10, tab, d10))
}

# plotQuantile(K=10, L=100, small=1e-2, seed=0)
plotQuantile <- function(K, L, xmax = 1, small, seed) {
    p <- p.gen(L, small, seed)
    lw <- sum(log(p[1:K]))
    pval <- p.rtp.qbeta.integrate(K, p)

    f1 <- function(u) fBetaQuantile(u, lw, K, L)
    f2 <- function(u) fGammaQuantile(u, lw, K, L)

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
        axes = FALSE,
    )

    abline(v = 0, lty = 1, col = "#757474", lwd = 0.25)
    abline(v = 1, lty = 1, col = "#757474", lwd = 0.25)
    abline(v = 0.5, lty = 1, col = "#757474", lwd = 0.25)
    abline(h = 0, lty = 1, col = "#757474", lwd = 0.25)
    abline(h = 0.5, lty = 1, col = "#757474", lwd = 0.25)
}

# plotBxGintegrand(K=10, L=100, small=1e-2, seed=0)
plotBxGintegrand <- function(K, L, xmin = -1, xmax = 0, small, seed) {
    p <- p.gen(L, small, seed)
    pval <- p.rpt.dbeta.cuba(K, p)
    lw <- sum(log(p[1:K]))
    top <- fBetaTop(lw, K, L)
    left <- exp(lw / K) * 1.5
    # left <- exp((K - 1 + lw) / K) # Gamma mode

    right <- K / (L - 1)
    hight <- dbeta(K / (L - 1), K + 1, L - K)

    sd <- betaSD(K + 1, L - K)
    if (xmax <= 0) xmax <- min(1.0, top + 8 * sd)
    if (xmin < 0) xmin <- max(0, top - 6 * sd)

    f1 <- function(u) 1 - pgamma(K * log(u) - lw, K)
    f2 <- function(u) fBetaDensity(u, lw, K, L)
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

# plotIntegrandLocation(K=10, L=100, small=1e-1, seed=0)
plotIntegrandLocation <- function(K, L, xmin = -1, xmax = 0, small, seed) {
    p <- p.gen(L, small, seed)
    lw <- sum(log(p[1:K]))
    pval <- p.rpt.dbeta.cuba(K, p)

    top <- fBetaTop(lw, K, L)
    left <- exp(lw / K) * 1.5
    right <- K / (L - 1)
    # left <- exp((K - 1 + lw) / K) # Gamma mode


    sd <- betaSD(K + 1, L - K)
    if (xmax <= 0) xmax <- min(1.0, top + 6 * sd)
    if (xmin < 0) xmin <- max(0.0, top - 6 * sd)

    f <- function(u) 1 - pgamma(K * log(u) - lw, K)
    f2 <- function(u) fBetaDensity(u, lw, K, L)
    f3 <- function(u) dbeta(u, K + 1, L - K)

    plot(f,
        type = "l", font.main = 1, lwd = 1, cex.main = 1, col = "red",
        xlim = c(xmin, xmax), yaxp = c(0, 1, 1),
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
            ),
            paste(
                "Left = ", format(left, digits = 2),
                " Integrand = ", format(top, digits = 2),
                " Right = ", format(right, digits = 2)
            )
        )
    )
    points(
        x = c(max(0, top - 3 * sd), top + 3 * sd, xmax),
        y = c(0, 0, 0.02), type = "p",
        pch = c(4, 4, 25), lwd = 2,
        col = c("blue", "blue", "blue"), bg = c("blue")
    )
    abline(v = left, col = "red", lty = 2, lwd = 0.5)
    abline(v = right, col = "#088308", lty = 2, lwd = 0.5)
    abline(v = top, col = "blue", lty = 2, lwd = 0.5)
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
plotGxBintegrand <- function(K, L, xmin = -1, xmax = 0, small = 1e-1, seed = 0) {
    p <- p.gen(L, small, seed)
    pval <- p.rpt.dbeta.cuba(K, p, 1e-10)
    lw <- sum(log(p[1:K]))
    hight <- dgamma(K - 1, K)

    f3 <- function(g) dgamma(g, K)
    f2 <- function(g) fGammaDensity(g, lw, K, L)
    f1 <- function(g) {
        u <- exp((g + lw) / K)
        pbeta(u, K + 1, L - K)
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
    abline(v = 0, lty = 1, col = "#757474", lwd = 0.25)

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


# bench(K=5, L=100, small=1e-5, seed=0)
bench <- function(K, L, small, seed) {
    library(microbenchmark)
    p <- p.gen(L, small, seed)
    Qinteg <- function() p.rtp.qbeta.integrate(K, p)
    Binteg <- function() p.rtp.dbeta.integrate(K, p)
    QsimpA <- function() p.rtp.qbeta.simp.a(K, p)
    BsimpA <- function() p.rtp.dbeta.simp.a(K, p)
    Briema <- function() p.rtp.dbeta.riema(K, p)
    Griema <- function() p.rtp.dgamma.riema(K, p, stepscale = 1)
    Gsimp <- function() p.rtp.dgamma.simp(K, p, stepscale = 1)
    Bcuba <- function() p.rpt.dbeta.cuba(K, p)
    # TFish <- function() p.tfisher.soft.R(0.05, p)
    TFish <- function() p.tfisher.soft(0.05, p)
    Bmutos <- function() ranktruncated(K, p)
    Art <- function() p.art(K, p)

    res <- microbenchmark(
        unit = "us",
        # Bmutos(),
        Qinteg(),
        QsimpA(),
        # Dinteg(),
        BsimpA(),
        Briema(),
        Gsimp(),
        Griema(),
        # Bcuba(),
        # TFish(),
        # Art(),
        times = 500
    )
    boxplot(res,
        font.main = 1, cex.main = 1,
        xlab = "Xyyyyy = Integrand function and integrating function",
        font.xlab = 1,
        # names = c(
        #     "Qinteg", "Qsimpa", "Dinteg", "Dcuba",
        #     "Dsimpa", "Driema"
        # ),
        main = paste(
            "Rank truncated p-value integrating time\n",
            "(microseconds, logarithmic time)"
        )
    )
    print(res, "us", signif = 3)
}