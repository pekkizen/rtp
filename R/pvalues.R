
# Used in plot.R and bench.R also.
checkPar <- function(K, L, small) {
    if (small < 1e-307) {
        return("Invalid: small < 1e-307")
    }
    if (K > L) {
        return("Invalid: K > L")
    }
    if (L <= 1) {
        return("Invalid: L <= 1")
    }
    if (K > 5000) {
        return("Invalid: K > 5000")
    }
    if (L > 1E+7) {
        return("Invalid: L > 1E+7")
    }
    if (K * L > 1E+8) {
        return("Invalid: K * L > 1E+8")
    }
    return("")
}

# pvalues.rtp(K=10, L=200, small=1e-6, seed=0)
pvalues.rtp <- function(K, L, small = 1, seed = 0, rounds = 5e5) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    if (seed > 0) set.seed(seed)
    if (small == 1) small <- 1 / L * 1e-4
    p <- c(small, runif(L - 1))

    p1 <- p.rtp.qbeta(K, p)
    p4 <- p.rtp.dbeta(K, p)
    p6 <- p.rtp.mutoss(K, p)

    p8 <- p.rtp.garp(K, p)
    p11 <- p.rtp.dgamma(K, p)
    # p13 <- p.rtp(K, p)
    p13 <- p.rtp.simulated(K, p, rounds = rounds)
    p12 <- p.rtp.dbeta.asimp(K, p, tol = 1e-4)
    # p12 <- 1
    pe <- p.rtp.dbeta.cuba(K, p, tol = 1e-14) # "exact" reference

    w <- function(s) writeLines(s)
    wl <- function(s, f, e, d) writeLines(paste(s, f, "  ", e, "  ", d))
    fpv <- function(p) {
        if (is.na(p)) {
            return("NaN")
        }
        f <- "%1.10f"
        if (p < 0.00001) f <- "%1.6e"
        sprintf(f, p)
    }
    fdig <- function(p) sprintf("%1.2f", (-log10(abs(p - pe) / p)))
    ferr <- function(p) {
        c <- abs(pe - p) / pe
        sprintf("%1.0e", c)
    }

    f1 <- fpv(p1)
    f4 <- fpv(p4)
    f6 <- fpv(p6)
    f8 <- fpv(p8)
    f11 <- fpv(p11)
    f12 <- fpv(p12)
    f13 <- fpv(p13)

    fe <- fpv(pe)

    e1 <- ferr(p1)
    e4 <- ferr(p4)
    e6 <- ferr(p6)
    e8 <- ferr(p8)
    e11 <- ferr(p11)
    e12 <- ferr(p12)
    e13 <- ferr(p13)

    d1 <- fdig(p1)
    d4 <- fdig(p4)
    d6 <- fdig(p6)
    d8 <- fdig(p8)
    d11 <- fdig(p11)
    d12 <- fdig(p12)
    d13 <- fdig(p13)

    w("\n                                           #nonzero")
    w("                                  relative   digits")
    w("RTP function----------P-value-------error----right")
    wl("p.rtp.qbeta        ", f1, e1, d1)
    wl("mutoss/ranktrunc   ", f6, e6, d6)
    wl("\np.rpt.dbeta.cuba   ", fe, "0 (ref)   ~15", "\n")
    wl("p.rtp.garp         ", f8, e8, d8)
    wl("p.rtp.dbeta        ", f4, e4, d4)
    wl("p.rtp.dgamma       ", f11, e11, d11)
    wl("p.rtp.dbeta.asimp  ", f12, e12, d12)
    wl("p.rtp.simulated    ", f13, e13, d13)
}

# pvalues.methods(K=10, L=100, small=1e-6, seed=0)
pvalues.methods <- function(K, L, small = 1e-1, seed = 0) {
    err <- checkPar(K, L, small)
    if (err != "") {
        return(err)
    }
    if (seed > 0) set.seed(seed)
    p <- c(small, runif(L - 1))
    tau <- (K + 1) / (L + 1) # beta mean

    art <- p.art(K, p)
    rtp <- p.rtp.dbeta.cuba(K, p, tol = 1e-10)

    tsoft <- p.tfisher.soft(tau, p)
    ttpm <- p.tfisher.tpm(tau, p)
    fish <- p.fisher(p)
    garp <- p.rtp.garp(K, p)

    w <- function(s, p, e) writeLines(paste(s, p, e))
    fd <- function(p) paste("   ", format(p / rtp, digits = 4))
    fpv <- function(p) {
        f <- "%1.6f"
        if (p < 0.00001) f <- "%1.4e"
        sprintf(f, p)
    }
    w("\nMethod----------P-value-----P-value/RTP", "", "")
    w("RTP            ", fpv(rtp), "    1")
    w("Garp           ", fpv(garp), fd(garp))
    w("TFisher soft   ", fpv(tsoft), fd(tsoft))
    w("TFisher tpm    ", fpv(ttpm), fd(ttpm))
    w("Fisher         ", fpv(fish), fd(fish))
    w("Art            ", fpv(art), fd(art))
}