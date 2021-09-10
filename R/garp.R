

p.rtp.garpR <- function(K, p) {
    L <- length(p)
    lw <- stat.rtp(K, p)
    mean <- meanRTP(K, L)
    var <- varRTP(K, L)
    theta <- var / mean
    shape <- mean * mean / var
    pgarp(lw, shape, theta, lower.tail = F)
}

# Garp distribution function.
pgarp <- function(x, shape, theta, lower.tail = T) {
    x <- gfun(x, shape, theta)
    pgamma(x, shape, scale = theta, lower.tail = lower.tail)
}

# Garp density function.
dgarp <- function(x, shape, theta) {
    g <- gfun(x, shape, theta)
    dgamma(g, shape, scale = theta) * gder(x, shape, theta)
}

# Garp quantile function
qgarp <- function(p, shape, theta) {
    if (p <= 0) {
        return(0)
    }
    if (p >= 1) {
        return(Inf)
    }
    q <- qgamma(p, shape, scale = theta, lower.tail = T)
    gfunInv(q, shape, theta)
}

gfun <- function(lw, shape, theta) {
    b <- c(
        2.06193e-01, 1.00015e+00, -2.35081e-01, -5.24942e-05,
        8.73881e-02, -2.27766e-01, 1.64767e-02, -7.98814e-02,
        2.53813e-01, -1.85481e-02, -2.87821e-03, 9.53178e-04,
        -3.93185e-04
    )
    mean <- shape * theta
    var <- mean * theta
    sqsh <- sqrt(shape)
    d <- (lw - mean) / sqrt(var)

    d <- b[1] + b[2] * lw + b[3] * theta + b[4] * shape +
        (b[5] + b[8] * theta + b[11] * sqsh) * d +
        (b[6] + b[9] * theta + b[12] * sqsh) * d * d +
        (b[7] + b[10] * theta + b[13] * sqsh) * d * d * d
    return(ifelse(d <= 0, 0, d))
}

gder <- function(x, shape, theta) {
    b <- c(
        2.06193e-01, 1.00015e+00, -2.35081e-01, -5.24942e-05,
        8.73881e-02, -2.27766e-01, 1.64767e-02, -7.98814e-02,
        2.53813e-01, -1.85481e-02, -2.87821e-03, 9.53178e-04,
        -3.93185e-04
    )
    mean <- shape * theta
    var <- mean * theta
    sqsh <- sqrt(shape)
    sd <- sqrt(var)
    d <- (x - mean) / sd
    d <- b[2] +
        (b[5] + b[8] * theta + b[11] * sqsh) / sd +
        (b[6] + b[9] * theta + b[12] * sqsh) / sd * 2 * d +
        (b[7] + b[10] * theta + b[13] * sqsh) / sd * 3 * d * d
    return(ifelse(d < 1e-300, 1e-300, d))
}

# gfun inverse function by Newton-Raphson iteration.
gfunInv <- function(y, shape, theta, tol = 1e-5) {
    x <- y
    for (i in 1:10) {
        delta <- (y - gfun(x, shape, theta)) / gder(x, shape, theta)
        x <- x + delta
        if (abs(delta) < tol) {
            return(x)
        }
    }
    NaN
}