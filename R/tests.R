

# p <- p.gen(K=10, L = 100, small = 1e-4, seed = 0)
p.gen <- function(K, L, small, seed = 0) {
    if (seed > 0) set.seed(seed)
    sort(c(small, runif(L - 1)), partial = c(1:K))
}

# testSelect(100, 1000, 1, 123)
# This tests uniSelect(K, p) vs. R sort(p, partial = c(1:K))
testSelect <- function(kmax, lmax, kstep, seed) {
    if (seed > 0) set.seed(seed)
    K <- 2
    rounds <- 0
    while (K <= kmax) {
        L <- 2 * K
        while (L <= lmax) {
            rounds <- rounds + 1
            p <- c(runif(L))
            psp <- sort(p, partial = c(1:K))[1:K]
            # p[1] <- 1e-10
            init(K, p) # calls uniSelect(K, p)
            pus <- sort(p[1:K])
            if (!identical(psp, pus)) {
                return(writeLines(paste("Round = ", format(rounds), " not same")))
            }
            L <- L + 2
        }
        K <- K + kstep
    }
    writeLines(paste("Rounds = ", format(rounds), " all same"))
}

terminalCommands <- function() {

    # Vsevolozhskaya et al (2019) S-3
    # "The resulting combined P-values are ART=0.045 and PRTP=0.047".
    p <- sort(c(0.7, 0.07, 0.15, 0.12, 0.08, 0.09))
    K <- 4
    L <- 7
    p.art(K, p) # 0.04487285
    p.rtp.dbeta.cuba(K, p) # 0.04741096

    # Test Fisher/TFisher R and C++ functions equivalence.
    K <- 10
    L <- 100
    p <- p.gen(K, L, small = 1e-4, seed = 0)
    tau <- (K + 1) / (L + 1)

    sprintf("%1.16e", p.tfisher.soft(tau, p))
    sprintf("%1.16e", p.tfisher.softR(tau, p))

    sprintf("%1.16e", p.tfisher.soft(tau, p, tol = 0))
    sprintf("%1.16e", p.tfisher.soft(tau, p))

    sprintf("%1.16f", p.tfisher.tpm(tau, p))
    sprintf("%1.16f", p.tfisher.tpmR(tau, p))

    sprintf("%1.16f", fisher(p))
    sprintf("%1.16f", p.fisher(p))

    p <- p.gen(K = 10, L = 100, small = 1e-5, seed = 0)
    K <- 10
    L <- 100
    p.rtp.dgamma.simp(K, p)
    p.rtp.dgamma.riema(K, p)
    p.rtp.qbeta.integrate(K, p)
    p.rtp.dbeta.cuba(K, p)
    p.rtp.dbeta.simp.a(K, p)
    p.rtp.dbeta.riema(K, p)
    ranktruncated(K, p)
    fisher(p)


    pvaluesRTP(K = 5, L = 100, small = 1e-6, seed = 0)
    pvaluesMethods(K = 10, L = 100, small = 1e-6, seed = 0)

    plotQuantile(K = 10, L = 100, small = 1e-1, seed = 0)
    plotBxGintegrand(K = 10, L = 100, small = 1e-3, seed = 0)
    plotGxBintegrand(K = 10, L = 100, small = 1e-3, seed = 0)
    plotIntegrandLocation(K = 10, L = 100, small = 1e-4, seed = 0)

    benchIntegrals(K = 5, L = 100, small = 1e-5)
    benchSelect(K = 25, L = 1000, times = 2000)
    benchIntegrands(K = 5, L = 100, small = 1e-3)
}