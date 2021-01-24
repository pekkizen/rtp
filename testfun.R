
p.gen <- function(n, small, seed = 0) {
    if (seed > 0) set.seed(seed)
    sort(c(small, runif(n - 1)))
}

# pvalues(K=10, L=100, small=1e-8, seed=0)
pvalues <- function(K, L, small, seed) {
    # tol <- defTol
    p <- p.gen(L, small, seed)
    p1 <- p.rtp.qbeta.integrate(K, p)
    p2 <- p.rtp.qbeta.simpa(K, p)
    p7 <- p.rtp.dbeta.integrate(K, p)
    p8 <- p.rtp.dbeta.simpa(K, p)
    p4 <- p.rtp.dbeta.riema(K, p)
    p3 <- p.rpt.dbeta.cuba(K, p)
    p10 <- p.art(K, p)
    p9 <- p.tfisher.soft(K / L, p)
    p6 <- ranktruncated(K, p)

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
    tab <- "  \t"
    writeLines("                                                   # of correct")
    writeLines("----Function-----------P-value--------Abs error----non zero digits")
    writeLines(paste("p.rtp.qbeta.integrate ", f1, tab, e1, tab, d1))
    writeLines(paste("p.rtp.qbeta.simpa     ", f2, tab, e2, tab, d2))
    writeLines(paste("p.rtp.dbeta.integrate ", f7, tab, e7, tab, d7))
    writeLines(paste("p.rpt.dbeta.cuba      ", f3, tab, e3, tab, d3))
    writeLines(paste("p.rtp.dbeta.simpa     ", f8, tab, e8, tab, d8))
    writeLines(paste("p.rtp.dbeta.riema     ", f4, tab, e4, tab, d4))

    writeLines(paste("mutoss/ranktruncated  ", f6, tab, e6, tab, d6))
    writeLines(paste("p.tfisher.soft        ", f9, tab, e9, tab, d9))
    writeLines(paste("p.art                 ", f10, tab, e10, tab, d10))
}

# plotQuantile(K=10, L=100, xmax = 1, small=1e-1, seed=0)
plotQuantile <- function(K, L, xmax = 1, small, seed) {
    p <- p.gen(L, small, seed)
    lw <- sum(log(p[1:K]))
    pval <- p.rtp.qbeta.integrate(K, p)

    f <- function(u) fQuantile(u, lw, K, L)
    plot.new()
    plot(f,
        type = "l", ylab = "", font.main = 1,
        col = "blue", cex.main = 1,
        xaxp = c(0, xmax, 2), xlim = c(0, xmax),
        yaxp = c(0, 1, 2), ylim = c(0, 1),
        xlab = paste(
            "Area under the blue integrand is ",
            "p = ", format(pval, digits = 2)
        ),
        main = "Integrand by Beta quantile function"
    )
    abline(v = 0, lty = 1, col = "#9e9d9d")
    abline(h = 0, lty = 1, col = "#9e9d9d")
    abline(v = 0.5, lty = 1, col = "#9e9d9d")
    abline(h = 0.5, lty = 1, col = "#9e9d9d")
}

# plotDensity(K=10, L=100, xmax=0.25, small=1e-4, seed=0)
plotDensity <- function(K, L, xmax, small, seed) {
    p <- p.gen(L, small, seed)
    pval <- p.rpt.dbeta.cuba(K, p)
    lw <- sum(log(p[1:K]))
    top <- fDenTop(lw, K, L)
    left <- exp(lw / K) * 1.5
    hight <- dbeta(K / (L - 1), K + 1, L - K)
    end <- min(1.0, top + 6 * betaSD(K + 1, L - K))

    f1 <- function(u) 1 - pgamma(K * log(u) - lw, K)
    f2 <- function(u) fDensity(u, lw, K, L)
    f3 <- function(u) dbeta(u, K + 1, L - K)

    plot.new()
    plot(f3,
        type = "l", lwd = 1, col = "#088308",
        xaxp = c(0, xmax, 10), xlim = c(0, xmax), ylim = c(0, hight),
        ylab = "", font.main = 1, cex.main = 1,
        xlab = c(
            paste(
                "Area under the blue integrand is p = ",
                format(pval, digits = 2)
            ),
            "Area under the green Beta density is 1"
        ),
        main = c(
            paste(
                "Beta density (green), ",
                "Gamma survival (red), ",
                "Beta x Gamma integrand (blue)."
            ),
            paste(
                " K = ", format(K, digits = 5),
                " L = ", format(L, digits = 5)
            )
        )
    )
    abline(v = left, lty = 2, col = "#f57575", lwd = 0.5)
    abline(v = K / (L - 1), lty = 2, col = "#0aa50a", lwd = 0.5)
    abline(v = top, lty = 2, col = "#5454f7", lwd = 0.5)
    abline(h = 0, lty = 1, col = "#9e9d9d", lwd = 1)
    abline(h = 1, lty = 1, col = "#9e9d9d", lwd = 1)

    par(new = TRUE)
    plot(f2,
        type = "l", lwd = 1, col = "blue",
        xlim = c(0, xmax), ylim = c(0, hight),
        yaxt = "n", xaxt = "n", ylab = "", xlab = ""
    )
    par(new = TRUE)
    plot(f1,
        type = "l", lwd = 1, col = "red",
        xlim = c(0, xmax), ylim = c(0, hight),
        xaxt = "n", xlab = "", ylab = "",
        if (L > 15) yaxp <- c(0, 1, 1) # stopped working?
    )
}

# plotLocation(K=10, L=100, xmax=0.3, small=1e-1, seed=0)
plotLocation <- function(K, L, xmax = 1, small, seed) {
    p <- p.gen(L, small, seed)
    lw <- sum(log(p[1:K]))
    pval <- p.rpt.dbeta.cuba(K, p)

    top <- fDenTop(lw, K, L)
    left <- exp(lw / K) * 1.5
    right <- K / (L - 1)
    end <- min(1.0, top + 6 * betaSD(K + 1, L - K))
    hwid <- 3 * betaSD(K + 1, L - K)

    f <- function(u) 1 - pgamma(K * log(u) - lw, K)
    f2 <- function(u) fDensity(u, lw, K, L)
    f3 <- function(u) dbeta(u, K + 1, L - K)

    plot.new()
    plot(f,
        type = "l", font.main = 1, lwd = 1, cex.main = 1,
        xaxp = c(0, xmax, 10), xlim = c(0, xmax), col = "red",
        yaxp = c(0, 1, 1),
        xlab = "Blue and green curves are vertically not in actual size",
        ylab = "",
        main = c(
            paste(
                "Gamma survival (red), ", "Beta density (green), ",
                "Beta x Gamma (blue)"
            ),
            paste(
                "K = ", format(K, digits = 5),
                " L = ", format(L, digits = 5),
                " p = ", format(pval, digits = 2)
            ),
            paste(
                "Left = ", format(left, digits = 2),
                " Top = ", format(top, digits = 2),
                " Right = ", format(right, digits = 2),
                " End = ", format(end, digits = 2)
            )
        )
    )
    points(
        x = c(max(0, top - hwid), top + hwid, end),
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
        type = "l", font.main = 1, lwd = 1,
        xaxp = c(0, xmax, 20), xlim = c(0, xmax),
        yaxt = "n", xaxt = "n", xlab = "", ylab = "",
        col = "#088308",
    )
    par(new = TRUE)
    plot(f2,
        type = "l", font.main = 1, lwd = 1,
        xaxp = c(0, xmax, 20), xlim = c(0, xmax),
        yaxt = "n", xaxt = "n", ylab = "", xlab = "",
        col = "blue",
    )
}

# bench(K=10, L=100, small=1e-8, seed=0)
bench <- function(K, L, small, seed) {
    library(microbenchmark)
    p <- p.gen(L, small, seed)
    Qinteg <- function() p.rtp.qbeta.integrate(K, p)
    Dinteg <- function() p.rtp.dbeta.integrate(K, p)
    Qsimpa <- function() p.rtp.qbeta.simpa(K, p)
    Dsimpa <- function() p.rtp.dbeta.simpa(K, p)
    Driema <- function() p.rtp.dbeta.riema(K, p)
    Dcuba <- function() p.rpt.dbeta.cuba(K, p)
    TFish <- function() p.tfisher.soft(K / L, p)
    mutoss <- function() ranktruncated(K, p)
    Art <- function() p.art(K, p)

    res <- microbenchmark(
        Qinteg(),
        Qsimpa(),
        Dinteg(),
        Dcuba(),
        Dsimpa(),
        Driema(),
        # sum(log(p[1:K])),
        # p <- p.gen(L, small, seed),
        # TFish(),
        # mutoss(),
        # Art(),
        times = 500
    )
    boxplot(res,
        font.main = 1, cex.main = 1.15,
        xlab = "Xyyyyy = Integrand and function",
        font.xlab = 1,
        # names = c(
        #     "Qinteg", "Qsimpa", "Dinteg", "Dcuba",
        #     "Dsimpa", "Driema", "TFish"
        # ),
        main = paste(
            "P-value integrating time (microseconds)\n",
            "(! logarithmic time)"
        )
    )
    res
}