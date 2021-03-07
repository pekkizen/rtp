

check <- function(K, L, small) {
    if (small < 1e-300) {
        return("Invalid: small < 1e-300")
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

# pvaluesRTP(K=10, L=100, small=1e-6, seed=0)
pvaluesRTP <- function(K, L, small, seed) {
    err <- check(K, L, small)
    if (err != "") {
        return(err)
    }
    if (seed > 0) set.seed(seed)
    p <- c(small, runif(L - 1))
    q <- sort(p)
    tau <- (K + 1) / (L + 1)

    p1 <- p.rtp.qbeta.integrate(K, q)
    p4 <- p.rtp.dbeta.riema(K, p, stepscale = 1)
    p6 <- ranktruncated(K, q)
    p8 <- p.rtp.dbeta.simp.a(K, p, abstol = 1e-7, reltol = 1e-3)
    p9 <- p.tfisher.soft(tau, p)
    p2 <- p.tfisher.tpm(tau, p)
    p10 <- p.art(K, p)
    p11 <- p.rtp.dgamma.riema(K, p, tol = 1e-10, stepscale = 1)
    p13 <- p.rtp.dgamma.simp(K, p, tol = 1e-10, stepscale = 1)

    pe <- p.rtp.dbeta.cuba(K, p, tol = 1e-15) # "exact" reference

    w <- function(s) writeLines(s)
    wl <- function(s, f, e, d) writeLines(paste(s, f, "  ", e, "   ", d))
    fpv <- function(p) {
        if (is.na(p)) {
            return("NaN")
        }
        f <- "%1.10f"
        if (p < 0.00001) f <- "%1.6e"
        sprintf(f, p)
    }
    fdig <- function(p) format(-log10(abs(p - pe) / pe), digits = 2)
    ferr <- function(p) format(abs(p - pe), digits = 1)

    f1 <- fpv(p1)
    f2 <- fpv(p2)
    f4 <- fpv(p4)
    f6 <- fpv(p6)
    f8 <- fpv(p8)
    f9 <- fpv(p9)
    f10 <- fpv(p10)
    f11 <- fpv(p11)
    f13 <- fpv(p13)

    fe <- fpv(pe)

    e1 <- ferr(p1)
    e2 <- ferr(p2)
    e4 <- ferr(p4)
    e6 <- ferr(p6)
    e8 <- ferr(p8)
    e9 <- ferr(p9)
    e10 <- ferr(p10)
    e11 <- ferr(p11)
    e13 <- ferr(p13)

    d1 <- fdig(p1)
    d2 <- fdig(p2)
    d4 <- fdig(p4)
    d6 <- fdig(p6)
    d8 <- fdig(p8)
    d9 <- fdig(p9)
    d10 <- fdig(p10)
    d11 <- fdig(p11)
    d13 <- fdig(p13)

    w("\n                                           # nonzero")
    w("                                     Abs     digits")
    w("RTP function----------P-value-------error----right")
    wl("p.rtp.qbeta.integ  ", f1, e1, d1)
    wl("mutoss/ranktrunc   ", f6, e6, d6)
    wl("\np.rpt.dbeta.cuba   ", fe, "0 (ref)   ~14", "\n")
    wl("p.rtp.dbeta.simp.a ", f8, e8, d8)
    wl("p.rtp.dbeta.riema  ", f4, e4, d4)
    wl("p.rtp.dgamma.simp  ", f13, e13, d13)
    wl("p.rtp.dgamma.riema ", f11, e11, d11)
    # wl("\np.tfisher.soft     ", f9, e9, d9)
    # wl("p.tfisher.tpm      ", f2, e2, d2)
    # wl("p.art              ", f10, e10, d10)
    # w(paste("p.fisher           ", fpv(p.fisher(p))))
}

# pvaluesMethods(K=10, L=100, small=1e-6, seed=0)
pvaluesMethods <- function(K, L, small, seed) {
    err <- check(K, L, small)
    if (err != "") {
        return(err)
    }
    if (seed > 0) set.seed(seed)
    p <- c(small, runif(L - 1))
    q <- sort(p)
    tau <- (K + 1) / (L + 1)

    rtp <- p.rtp.dbeta.cuba(K, p)
    art <- p.art(K, p)
    tsoft <- p.tfisher.soft(tau, p)
    ttpm <- p.tfisher.tpm(tau, p)
    fish <- fisher(p)

    w <- function(s, p, e) writeLines(paste(s, p, e))
    # fd <- function(p) paste("   ", format((p - rtp) / rtp, digits = 4))
    fd <- function(p) paste("   ", format(p / rtp, digits = 4))
    fpv <- function(p) {
        f <- "%1.6f"
        if (p < 0.00001) f <- "%1.4e"
        sprintf(f, p)
    }

    w("\nFunction--------P-value-----P-value/RTP", "", "")
    w("RTP            ", fpv(rtp), "    1")
    w("TFisher soft   ", fpv(tsoft), fd(tsoft))
    w("Art            ", fpv(art), fd(art))
    w("TFisher tpm    ", fpv(ttpm), fd(ttpm))
    w("Fisher         ", fpv(fish), fd(fish))
}