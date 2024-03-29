# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

p.rtp.garp.lw <- function(lw, K, L) {
    .Call(`_rtp_rtpGarpLW`, lw, K, L)
}

p.rtp.garp <- function(K, p) {
    .Call(`_rtp_rtpGarp`, K, p)
}

rgarp <- function(samplesize, k, l) {
    .Call(`_rtp_randGarp`, samplesize, k, l)
}

garpError <- function(r1, r2, samples, plim, K, L) {
    .Call(`_rtp_garpError`, r1, r2, samples, plim, K, L)
}

baseNull <- function(x) {
    .Call(`_rtp_baseNull`, x)
}

stat.rtp <- function(k, q) {
    .Call(`_rtp_statRTP`, k, q)
}

init <- function(k, p, integrand = 0L) {
    .Call(`_rtp_init`, k, p, integrand)
}

betaSD <- function(a, b) {
    .Call(`_rtp_betaSD`, a, b)
}

survbinom <- function(k, n, p) {
    .Call(`_rtp_survbinom`, k, n, p)
}

survgamma <- function(g, k) {
    .Call(`_rtp_survgamma`, g, k)
}

fBetaD <- function(b) {
    .Call(`_rtp_fBetaD`, b)
}

fGammaD <- function(g) {
    .Call(`_rtp_fGammaD`, g)
}

fBetaDtop <- function() {
    .Call(`_rtp_fBetaDtop`)
}

p.rtp.dbeta.lw <- function(lw, k, l, tol = 1e-12, stepscale = 1) {
    .Call(`_rtp_rtpDbetaLW`, lw, k, l, tol, stepscale)
}

p.rtp.dbeta <- function(k, p, tol = 1e-10, stepscale = 1) {
    .Call(`_rtp_rtpDbeta`, k, p, tol, stepscale)
}

p.rtp.dbeta.asimp <- function(k, p, tol = 1e-4) {
    .Call(`_rtp_rtpDbetaAsimp`, k, p, tol)
}

p.rtp.dgamma.lw <- function(lw, k, l, tol = 1e-10, stepscale = 1) {
    .Call(`_rtp_rtpDgammaLW`, lw, k, l, tol, stepscale)
}

p.rtp.dgamma <- function(k, p, tol = 1e-10, stepscale = 1) {
    .Call(`_rtp_rtpDgamma`, k, p, tol, stepscale)
}

p.rtp.lw <- function(lw, k, l, tol = 1e-14, stepscale = 0.5) {
    .Call(`_rtp_rtpLW`, lw, k, l, tol, stepscale)
}

p.rtp <- function(k, p, tol = 1e-14, stepscale = 0.5) {
    .Call(`_rtp_rtp`, k, p, tol, stepscale)
}

p.rtp.vector <- function(w, k, l, tol = 1e-10, stepscale = 0.5) {
    .Call(`_rtp_rtpLWvec`, w, k, l, tol, stepscale)
}

drtp <- function(x, k, l, tol = 1e-16, stepscale = 1) {
    .Call(`_rtp_drtp`, x, k, l, tol, stepscale)
}

rrtp <- function(n, k, l) {
    .Call(`_rtp_rrtp`, n, k, l)
}

rtpSimulated <- function(k, p, samples) {
    .Call(`_rtp_rtpSimulated`, k, p, samples)
}

uniSel <- function(k, p) {
    .Call(`_rtp_uniSel`, k, p)
}

simpleSel <- function(k, p) {
    .Call(`_rtp_simpleSel`, k, p)
}

nth_element <- function(k, p) {
    .Call(`_rtp_nth_element`, k, p)
}

p.tfisher <- function(lw, L, tau1, tau2, tol = 1e-16) {
    .Call(`_rtp_tfisher`, lw, L, tau1, tau2, tol)
}

