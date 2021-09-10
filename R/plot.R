
# plot.quantile.integrands(K=10, L=100, small=1e-4, seed=0)
plot.quantile.integrands <- function(K, L, small = 1e-1, seed = 0) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    if (seed > 0) set.seed(seed)
    p <- c(small, runif(L - 1))
    lw <- stat.rtp(K, p)
    pval <- p.rtp.dbeta.cuba(K, p, 1e-10)

    f1 <- function(u) fBetaQ.R(u, lw, K, L)
    f2 <- function(u) fGammaQ.R(u, lw, K, L)

    plot.new()
    plot(f1,
        type = "l", font.main = 1,
        col = "blue", cex.main = 1,
        xlim = c(0, 1), ylim = c(0, 1),
        ylab = "", xaxt = "n", yaxt = "n",
        xlab = paste(
            "Area under the integrands is ",
            "p = ", format(pval, digits = 2)
        ),
        main = c(
            "Integrand by inverse Beta CDF (quantile function)"
        )
    )
    axis(1, at = c(0, 0.5, 1))
    axis(2, at = c(0, 0.5, 1))


    abline(v = 0, lty = 1, col = "#757474", lwd = 0.25)
    abline(v = 1, lty = 1, col = "#757474", lwd = 0.25)
    abline(v = 0.5, lty = 1, col = "#757474", lwd = 0.25)
    abline(h = 0, lty = 1, col = "#757474", lwd = 0.25)
    abline(h = 0.5, lty = 1, col = "#757474", lwd = 0.25)
}

# plot.BxG.integrand(K=10, L=100, small=1e-4, seed=0)
plot.BxG.integrand <- function(K, L, small = 1e-1, seed = 0) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    if (seed > 0) set.seed(seed)
    p <- c(small, runif(L - 1))
    lw <- stat.rtp(K, p)
    pval <- p.rtp.dbeta.cuba(K, p) # calls init
    top <- fBetaDtop()
    left <- exp(lw / K) * 1.5
    right <- K / (L - 1) # beta mode
    hight <- dbeta(right, K + 1, L - K) # function value at mode

    sd <- betaSD(K + 1, L - K)
    xmax <- min(1.0, top + 8 * sd)
    xmin <- max(0, top - 6 * sd)

    f1 <- function(u) 1 - pgamma(K * log(u) + lw, K)
    f2 <- function(u) fBetaD.R(u, lw, K, L)
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
    axis(1, at = pretty(c(xmin, xmax)))
    axis(2, at = pretty(0:floor(hight)), col.axis = "blue")

    par(new = TRUE)
    plot(f2,
        type = "l", lwd = 1, col = "blue",
        xlim = c(xmin, xmax), ylim = c(0, hight),
        ylab = "", xlab = "", axes = FALSE, yaxt = "n",
    )
    par(new = TRUE)
    plot(f1,
        type = "l", lwd = 1, col = "red",
        xlim = c(xmin, xmax), ylim = c(0, 1),
        ylab = "", xlab = "", axes = FALSE, yaxt = "n",
    )
    axis(4, col.axis = "red", at = c(0, 0.5, 1))

    abline(v = left, lty = 2, col = "red", lwd = 0.5)
    abline(v = right, lty = 2, col = "darkgreen", lwd = 0.5)
    abline(v = top, lty = 2, col = "blue", lwd = 0.5)
    abline(h = 0, lty = 1, col = "#757474", lwd = 0.25)
    abline(v = 0, lty = 1, col = "#757474", lwd = 0.25)
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
    pval <- p.rtp.dbeta.cuba(K, p) # calls init
    hight <- dgamma(K - 1, K)


    f3 <- function(g) dgamma(g, K)
    f2 <- function(g) fGammaD.R(g, lw, K, L)
    f1 <- function(g) {
        b <- exp((g - lw) / K)
        pbeta(b, K + 1, L - K)
    }
    bTop <- K * log(K / (L - 1)) - lw
    gTop <- K - 1
    xmax <- floor(max(bTop, gTop) + 4 * sqrt(K))
    xmin <- max(0, floor(K - 1 - 3 * sqrt(K)))
    plim <- K / L + 5 * betaSD(K, L - K + 1)
    plim <- log(plim) * K - lw
    xmax <- max(plim, xmax)

    plot.new()
    plot(f2,
        type = "l", lwd = 1, font.main = 1, cex.main = 1, col = "blue",
        xlim = c(xmin, xmax),
        ylim = c(0, hight),
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
                ", L = ", format(L, digits = 5)
            )
        )
    )
    axis(1, at = pretty(xmin:xmax))
    axis(2, at = pretty(c(0, hight)), col.axis = "blue")

    abline(v = gTop, lty = 2, col = "darkgreen", lwd = 0.5)
    abline(v = bTop, lty = 2, col = "red", lwd = 0.5)
    abline(h = 0, lty = 1, col = "#757474", lwd = 0.25)
    abline(v = plim, lty = 2, col = "#615959", lwd = 0.5)

    par(new = TRUE)
    plot(f1,
        type = "l", lwd = 1, col = "red",
        xlim = c(xmin, xmax), ylim = c(0, 1),
        xlab = "", ylab = "", axes = FALSE,
        yaxt = "n", xaxt = "n",
    )
    axis(4, col.axis = "red", at = c(0, 0.5, 1))
    par(new = TRUE)
    plot(f3,
        type = "l", lwd = 1, col = "#035203",
        xlim = c(xmin, xmax), ylim = c(0, hight),
        xlab = "", ylab = "", axes = FALSE,
        yaxt = "n", xaxt = "n",
    )
}

# This compares Gamma and Normal distributions
# plot.GammaDist(K=10)
plot.GammaDist <- function(K) {
    hight <- dgamma(K - 1, K)
    gTop <- K - 1
    xmax <- floor(gTop + 8 * sqrt(K))
    xmin <- max(0, floor(K - 1 - 3 * sqrt(K)))
    rpoint <- xmax

    f3 <- function(g) dgamma(g, K)
    f2 <- function(g) dnorm(g, K, sqrt(K))

    plot.new()
    plot(f2,
        type = "l", lwd = 1, font.main = 1, cex.main = 1, col = "blue",
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
            "Gamma(K, 1) PDF (green) and Norm(K, sd=sqrt(K)) PDF (blue)",
            paste(
                " K = ", format(K, digits = 5)
            )
        )
    )
    axis(1, at = pretty(c(xmin, xmax)))
    axis(2, at = pretty(c(0, hight)))

    par(new = TRUE)
    plot(f3,
        type = "l", lwd = 1, col = "darkgreen",
        xlim = c(xmin, xmax), ylim = c(0, hight),
        xlab = "", ylab = "", axes = FALSE,
        yaxt = "n", xaxt = "n",
    )
}

# plot.BetaDist(5, 10)
plot.BetaDist <- function(K, L) {
    betaSkewness <- function(a, b) {
        2 * (b - a) * sqrt(a + b + 1) / ((a + b + 2) * sqrt(a * b))
    }

    sd <- betaSD(K + 1, L - K)
    plim <- (K + 1) / (L + 1) + 3 * sd
    xmax <- min(1, plim + 5 * sd)
    bprob <- pbeta(plim, K + 1, L - K, lower.tail = F)
    mean <- (K + 1) / (L + 1)
    skew <- betaSkewness(K + 1, L - K)
    nprob <- pnorm(plim, mean, sd, lower.tail = F)
    bhight <- dbeta(K / (L - 1), K + 1, L - K)
    nhight <- dnorm(mean, mean, sd)
    hight <- max(bhight, nhight)

    f1 <- function(p) pbeta(p, K + 1, L - K)
    f2 <- function(p) dbeta(p, K + 1, L - K)
    f3 <- function(p) dnorm(p, mean, sd)


    plot.new()
    plot(f2,
        type = "l", lwd = 1, font.main = 1, cex.main = 1, col = "darkgreen",
        xlim = c(0, xmax), ylim = c(0, hight),
        ylab = "", yaxt = "n", xaxt = "n",
        main = c(
            paste("Beta CDF (red),", " Beta PDF (green)"),
            "Norm(betaMean, betaSD) PDF (blue)",
            paste(
                " K = ", format(K, digits = 5),
                ",  L = ", format(L, digits = 5)
                # ",  mean = ", format(mean, digits = 3)
            )
        ),
        xlab = c(
            paste(
                "Beta/Norm right tail probability from ", format(plim, digits = 3), "=",
                format(bprob, digits = 2), "/", format(nprob, digits = 2)
            ),
            paste(
                "mean =", format(mean, digits = 2),
                " SD =", format(sd, digits = 2),
                " skewness =", format(skew, digits = 2)
            )
        ),
    )
    axis(1, at = pretty(c(0, xmax)))
    axis(2, at = pretty(c(0, hight)), col.axis = "darkgreen")

    par(new = TRUE)
    plot(f3,
        type = "l", lwd = 1, col = "blue",
        xlim = c(0, xmax), ylim = c(0, hight),
        xlab = "", ylab = "",
        yaxt = "n", xaxt = "n", axes = FALSE,
    )
    par(new = TRUE)
    plot(f1,
        type = "l", lwd = 1, col = "red",
        xlim = c(0, xmax), ylim = c(0, 1),
        xlab = "", ylab = "",
        yaxt = "n", xaxt = "n", axes = FALSE,
    )
    axis(4, at = c(0, 0.5, 1), col.axis = "red")

    abline(h = 1, lty = 2, col = "#615959", lwd = 0.25)
    abline(v = plim, lty = 2, col = "#615959", lwd = 0.5)
    abline(v = mean, lty = 2, col = "blue", lwd = 0.5)
}