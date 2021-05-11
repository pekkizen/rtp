
# plot.quantile.integrands(K=10, L=100, small=1e-4, seed=0)
plot.quantile.integrands <- function(K, L, small = 1e-1, seed = 0) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    if (seed > 0) set.seed(seed)
    p <- c(small, runif(L - 1))
    lw <- stat.rtp(K, p)
    pval <- p.rtp.qbeta(K, p)

    f1 <- function(u) fBetaQ.R(u, lw, K, L)
    f2 <- function(u) fGammaQ.R(u, lw, K, L)

    plot.new()
    plot(f1,
        type = "l", ylab = "", font.main = 1,
        col = "blue", cex.main = 1,
        xlim = c(0, 1), ylim = c(0, 1), yaxp = c(0, 1, 2),
        axis(1,
            at = pretty(c(0, 1)),
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
        xlim = c(0, 1), ylim = c(0, 1),
        xlab = "",
        axes = FALSE,
    )
    abline(v = 0, lty = 1, col = "#757474", lwd = 0.25)
    abline(v = 1, lty = 1, col = "#757474", lwd = 0.25)
    abline(v = 0.5, lty = 1, col = "#757474", lwd = 0.25)
    abline(h = 0, lty = 1, col = "#757474", lwd = 0.25)
    abline(h = 0.5, lty = 1, col = "#757474", lwd = 0.25)
}

# plot.BxG.integrand(K=0, L=100, small=1e-4, seed=0)
plot.BxG.integrand <- function(K, L, small = 1e-1, seed = 0) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    if (seed > 0) set.seed(seed)
    p <- c(small, runif(L - 1))
    lw <- stat.rtp(K, p)
    pval <- p.rtp.dbeta.cuba(K, p)
    init(K, p)
    top <- fBetaDtop()
    left <- exp(lw / K) * 1.5
    right <- K / (L - 1) # beta mode
    hight <- dbeta(right, K + 1, L - K) # function value at mode

    sd <- betaSD(K + 1, L - K)
    xmax <- min(1.0, top + 8 * sd)
    xmin <- max(0, top - 6 * sd)

    f1 <- function(u) 1 - pgamma(K * log(u) - lw, K)
    f2 <- function(u) fBetaD.R(u, lw, K, L)
    f3 <- function(u) dbeta(u, K + 1, L - K)

    plot.new()
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
            "BxG integrand (blue) = Beta PDF (green) x Gamma survival (red)",
            paste(
                " K = ", format(K, digits = 5),
                " L = ", format(L, digits = 5)
            ),
            paste(
                " Approx. integrand top (blue vertical line) = ",
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
    abline(v = right, lty = 2, col = "darkgreen", lwd = 0.5)
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

plotIntegrandLocation <- function(K, L, small = 1e-1, seed = 0) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    if (seed > 0) set.seed(seed)
    p <- c(small, runif(L - 1))
    lw <- stat.rtp(K, p)
    init(K, p)
    pval <- p.rtp.dbeta.cuba(K, p)
    iTop <- fBetaDtop()
    bTop <- K / (L - 1)
    bSD <- betaSD(K + 1, L - K)
    gTop <- exp((K - 1 + lw) / K)
    gSD <- exp((sqrt(K) + lw) / K)

    xmax <- min(1.0, iTop + 6 * bSD)
    xmin <- max(0.0, iTop - 6 * bSD)

    f1 <- function(u) 1 - pgamma(K * log(u) - lw, K)
    f2 <- function(u) fBetaD.R(u, lw, K, L)
    f3 <- function(u) dbeta(u, K + 1, L - K)

    plot.new()
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
                "Gamma survival (red), ", "Beta PDF (green), ",
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

# plot.GxB.integrand(K=10, L=100, small=1e-4, seed=0)
plot.GxB.integrand <- function(K, L, small = 1e-1, seed = 0) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    if (seed > 0) set.seed(seed)
    p <- c(small, runif(L - 1))
    lw <- stat.rtp(K, p)
    pval <- p.rtp.dbeta.cuba(K, p)

    init(K, p)
    hight <- dgamma(K - 1, K)

    f3 <- function(g) dgamma(g, K)
    f2 <- function(g) fGammaD.R(g, lw, K, L)
    f1 <- function(g) {
        b <- exp((g + lw) / K)
        pbeta(b, K + 1, L - K)
    }
    bTop <- K * log(K / (L - 1)) - lw
    gTop <- K - 1
    xmax <- floor(max(bTop, gTop) + 4 * sqrt(K))
    xmin <- max(0, floor(K - 1 - 3 * sqrt(K)))

    plot.new()
    plot(f2,
        type = "l", lwd = 1, font.main = 1, cex.main = 1, col = "blue",
        xlim = c(xmin, xmax), ylim = c(0, hight),
        yaxt = "n", xaxt = "n",
        ylab = "Integrand  and  Gamma  PDF",

        xlab = c(
            paste(
                "Area under the blue integrand is p = ",
                format(pval, digits = 2)
            ),
            "Area under the green Gamma PDF is 1"
        ),
        main = c(
            "GxB integrand (blue) = Gamma PDF (green) x Beta CDF (red)",
            paste(
                " K = ", format(K, digits = 5),
                " L = ", format(L, digits = 5)
            )
        )
    )
    axis(1,
        col.axis = "#3a3939",
        at = pretty(c(xmin, xmax))
    )
    axis(2,
        at = pretty(c(0, hight / 2, hight)),
        col.axis = "blue",
    )
    abline(v = gTop, lty = 2, col = "darkgreen", lwd = 0.5)
    abline(v = bTop, lty = 2, col = "red", lwd = 0.5)
    abline(h = 0, lty = 1, col = "#757474", lwd = 0.25)
    par(new = TRUE)
    plot(f1,
        type = "l", lwd = 1, col = "red",
        xlim = c(xmin, xmax), ylim = c(0, 1),
        xlab = "", ylab = "", axes = FALSE,
    )
    axis(4,
        at = c(0, 0.5, 1),
        col.axis = "red",
    )
    par(new = TRUE)
    plot(f3,
        type = "l", lwd = 1, col = "#035203",
        xlim = c(xmin, xmax), ylim = c(0, hight),
        xlab = "", ylab = "", axes = FALSE,
    )
}

# This compares Gamma and Normal distributions
# plot.GammaNorm(K=10)
plot.GammaNorm <- function(K) {
    f3 <- function(g) dgamma(g, K)
    f2 <- function(g) dnorm(g, K, sqrt(K))

    hight <- dgamma(K - 1, K)
    gTop <- K - 1
    xmax <- floor(gTop + 8 * sqrt(K))
    xmin <- max(0, floor(K - 1 - 3 * sqrt(K)))
    rpoint <- xmax

    plot.new()
    plot(f2,
        type = "l", lwd = 1, font.main = 1, cex.main = 1, col = "red",
        xlim = c(xmin, xmax), ylim = c(0, hight),
        ylab = "", yaxt = "n", xaxt = "n",
        xlab = c(
            paste(
                "Gamma density at ", format(rpoint, digits = 2), " :",
                format(dgamma(rpoint, K), digits = 2)
            ),
            paste(
                "Normal density at ", format(rpoint, digits = 2), " :",
                format(dnorm(rpoint, K, sqrt(K)), digits = 2)
            )
        ),
        main = c(
            "Gamma(K, 1) PDF (green) and Normal(K, sd=sqrt(K)) PDF (red)",
            paste(
                " K = ", format(K, digits = 5)
            )
        )
    )
    axis(1,
        at = pretty(c(xmin, xmax))
    )
    axis(2,
        at = pretty(c(0, hight / 2, hight)),
    )
    par(new = TRUE)
    plot(f3,
        type = "l", lwd = 1, col = "darkgreen",
        xlim = c(xmin, xmax), ylim = c(0, hight),
        xlab = "", ylab = "", axes = FALSE,
    )
}

# plot.BetaPDF(10, 200)
plot.BetaPDF <- function(K, L) {
    f1 <- function(p) pbeta(p, K + 1, L - K)
    f2 <- function(p) dbeta(p, K + 1, L - K)

    plim <- betaCutPoint(K, L)
    xmax <- min(1, plim + 5 * betaSD(K + 1, L - K))
    prob <- pbeta(plim, K + 1, L - K, lower.tail = FALSE)
    hight <- dbeta(K / (L - 1), K + 1, L - K)
    plot.new()
    plot(f1,
        type = "l", lwd = 1, font.main = 1, cex.main = 1, col = "#e91b1b",
        xlim = c(0, xmax), ylim = c(0, 1),
        ylab = "", yaxt = "n", xaxt = "n",
        main = c(
            paste("Beta CDF (red),", " Beta PDF (green)"),
            paste(
                " K = ", format(K, digits = 5),
                " L = ", format(L, digits = 5)
            )
        ),
        xlab = c(
            paste(
                "Beta right tail probability from ", format(plim, digits = 2), "=",
                format(prob, digits = 2)
            )
        ),
    )
    axis(1,
        at = pretty(c(0, xmax))
    )
    axis(4,
        at = c(0, 0.5, 1),
        col.axis = "red",
    )
    abline(h = 1, lty = 2, col = "#615959", lwd = 0.25)
    abline(v = plim, lty = 2, col = "#615959", lwd = 0.5)
    par(new = TRUE)
    plot(f2,
        type = "l", lwd = 1, col = "darkgreen",
        xlab = "", ylab = "",
        axes = FALSE,
        xlim = c(0, xmax),
    )
    axis(2,
        at = pretty(c(0, hight / 2, hight)),
        col.axis = "darkgreen",
    )
}