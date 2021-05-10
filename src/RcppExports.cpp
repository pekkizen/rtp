// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// baseNull
double baseNull(double x);
RcppExport SEXP _rtp_baseNull(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(baseNull(x));
    return rcpp_result_gen;
END_RCPP
}
// init
double init(int k, NumericVector p, int denfunc);
RcppExport SEXP _rtp_init(SEXP kSEXP, SEXP pSEXP, SEXP denfuncSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type denfunc(denfuncSEXP);
    rcpp_result_gen = Rcpp::wrap(init(k, p, denfunc));
    return rcpp_result_gen;
END_RCPP
}
// betaSD
double betaSD(double a, double b);
RcppExport SEXP _rtp_betaSD(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(betaSD(a, b));
    return rcpp_result_gen;
END_RCPP
}
// survbinom
double survbinom(double k, double n, double p);
RcppExport SEXP _rtp_survbinom(SEXP kSEXP, SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(survbinom(k, n, p));
    return rcpp_result_gen;
END_RCPP
}
// survgamma
double survgamma(double g, double k);
RcppExport SEXP _rtp_survgamma(SEXP gSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type g(gSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(survgamma(g, k));
    return rcpp_result_gen;
END_RCPP
}
// fBetaD
double fBetaD(double b);
RcppExport SEXP _rtp_fBetaD(SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(fBetaD(b));
    return rcpp_result_gen;
END_RCPP
}
// fGammaD
double fGammaD(double g);
RcppExport SEXP _rtp_fGammaD(SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(fGammaD(g));
    return rcpp_result_gen;
END_RCPP
}
// fBetaQ
double fBetaQ(double p);
RcppExport SEXP _rtp_fBetaQ(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(fBetaQ(p));
    return rcpp_result_gen;
END_RCPP
}
// fGammaQ
double fGammaQ(double p);
RcppExport SEXP _rtp_fGammaQ(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(fGammaQ(p));
    return rcpp_result_gen;
END_RCPP
}
// fBetaDtop
double fBetaDtop();
RcppExport SEXP _rtp_fBetaDtop() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(fBetaDtop());
    return rcpp_result_gen;
END_RCPP
}
// rtpDbetaRiema
double rtpDbetaRiema(int k, NumericVector p, double tol, double stepscale);
RcppExport SEXP _rtp_rtpDbetaRiema(SEXP kSEXP, SEXP pSEXP, SEXP tolSEXP, SEXP stepscaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type stepscale(stepscaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rtpDbetaRiema(k, p, tol, stepscale));
    return rcpp_result_gen;
END_RCPP
}
// rtpDbetaAsimp
double rtpDbetaAsimp(int k, NumericVector p, double abstol, double reltol);
RcppExport SEXP _rtp_rtpDbetaAsimp(SEXP kSEXP, SEXP pSEXP, SEXP abstolSEXP, SEXP reltolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type abstol(abstolSEXP);
    Rcpp::traits::input_parameter< double >::type reltol(reltolSEXP);
    rcpp_result_gen = Rcpp::wrap(rtpDbetaAsimp(k, p, abstol, reltol));
    return rcpp_result_gen;
END_RCPP
}
// rtpDgammaRiema
double rtpDgammaRiema(int k, NumericVector p, double tol, double stepscale);
RcppExport SEXP _rtp_rtpDgammaRiema(SEXP kSEXP, SEXP pSEXP, SEXP tolSEXP, SEXP stepscaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type stepscale(stepscaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rtpDgammaRiema(k, p, tol, stepscale));
    return rcpp_result_gen;
END_RCPP
}
// rtpDgammaSimp
double rtpDgammaSimp(int k, NumericVector p, double tol, double stepscale);
RcppExport SEXP _rtp_rtpDgammaSimp(SEXP kSEXP, SEXP pSEXP, SEXP tolSEXP, SEXP stepscaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type stepscale(stepscaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rtpDgammaSimp(k, p, tol, stepscale));
    return rcpp_result_gen;
END_RCPP
}
// rtpRiema
double rtpRiema(int k, NumericVector p, double tol, double stepscale);
RcppExport SEXP _rtp_rtpRiema(SEXP kSEXP, SEXP pSEXP, SEXP tolSEXP, SEXP stepscaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type stepscale(stepscaleSEXP);
    rcpp_result_gen = Rcpp::wrap(rtpRiema(k, p, tol, stepscale));
    return rcpp_result_gen;
END_RCPP
}
// uniSel
int uniSel(int k, NumericVector p);
RcppExport SEXP _rtp_uniSel(SEXP kSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(uniSel(k, p));
    return rcpp_result_gen;
END_RCPP
}
// simpleSel
int simpleSel(int k, NumericVector p);
RcppExport SEXP _rtp_simpleSel(SEXP kSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(simpleSel(k, p));
    return rcpp_result_gen;
END_RCPP
}
// nth_element
int nth_element(int k, NumericVector p);
RcppExport SEXP _rtp_nth_element(SEXP kSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(nth_element(k, p));
    return rcpp_result_gen;
END_RCPP
}
// tfisher
double tfisher(double lw, double L, double tau1, double tau2, double tol);
RcppExport SEXP _rtp_tfisher(SEXP lwSEXP, SEXP LSEXP, SEXP tau1SEXP, SEXP tau2SEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lw(lwSEXP);
    Rcpp::traits::input_parameter< double >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type tau1(tau1SEXP);
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(tfisher(lw, L, tau1, tau2, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rtp_baseNull", (DL_FUNC) &_rtp_baseNull, 1},
    {"_rtp_init", (DL_FUNC) &_rtp_init, 3},
    {"_rtp_betaSD", (DL_FUNC) &_rtp_betaSD, 2},
    {"_rtp_survbinom", (DL_FUNC) &_rtp_survbinom, 3},
    {"_rtp_survgamma", (DL_FUNC) &_rtp_survgamma, 2},
    {"_rtp_fBetaD", (DL_FUNC) &_rtp_fBetaD, 1},
    {"_rtp_fGammaD", (DL_FUNC) &_rtp_fGammaD, 1},
    {"_rtp_fBetaQ", (DL_FUNC) &_rtp_fBetaQ, 1},
    {"_rtp_fGammaQ", (DL_FUNC) &_rtp_fGammaQ, 1},
    {"_rtp_fBetaDtop", (DL_FUNC) &_rtp_fBetaDtop, 0},
    {"_rtp_rtpDbetaRiema", (DL_FUNC) &_rtp_rtpDbetaRiema, 4},
    {"_rtp_rtpDbetaAsimp", (DL_FUNC) &_rtp_rtpDbetaAsimp, 4},
    {"_rtp_rtpDgammaRiema", (DL_FUNC) &_rtp_rtpDgammaRiema, 4},
    {"_rtp_rtpDgammaSimp", (DL_FUNC) &_rtp_rtpDgammaSimp, 4},
    {"_rtp_rtpRiema", (DL_FUNC) &_rtp_rtpRiema, 4},
    {"_rtp_uniSel", (DL_FUNC) &_rtp_uniSel, 2},
    {"_rtp_simpleSel", (DL_FUNC) &_rtp_simpleSel, 2},
    {"_rtp_nth_element", (DL_FUNC) &_rtp_nth_element, 2},
    {"_rtp_tfisher", (DL_FUNC) &_rtp_tfisher, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_rtp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
