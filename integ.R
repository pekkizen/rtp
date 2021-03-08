
riemannR <- function(lw, k, l, tol, upsteps) {
    top <- k / (l - 1)
    h <- top / upsteps
    I <- 0
    top <- 2 * top
    a <- h / 2
    while (a < 1) {
        d <- h * fDensity(a, lw, k, l)
        I <- I + d
        if (a > top & d < tol) {
            return(I)
        }
        a <- a + h
    }
    return(I)
}

adaSimpR <- function(f, a, b, fa, fm, fb, sPrev, tol, depth) {
    h <- (b - a) / 4
    # lm <- a + h
    # rm <- b - h
    fam <- f(a + h)
    fbm <- f(b - h)
    sl <- (fa + 4 * fam + fm) * (h / 3)
    sr <- (fm + 4 * fbm + fb) * (h / 3)

    sNew <- sl + sr
    delta <- (sNew - sPrev) / 15
    if (depth <= 0 || abs(delta) < tol) {
        return(sNew + delta)
    }
    m <- (a + b) / 2
    sl <- adaSimpR(f, a, m, fa, fam, fm, sl, tol / 2, depth - 1)
    sr <- adaSimpR(f, m, b, fm, fbm, fb, sr, tol / 2, depth - 1)
    return(sl + sr)
}