
#include <Rcpp.h>
using namespace Rcpp;

/* See referencies in the README.md.

The K+1'th smallest of L unif(0, 1) numbers is distributed Beta(K+1, L-K).
The density function of order statistic x ~ Beta(K+1, L-K) is
    dbeta(x, K+1, L-K) = L! / (K! * L-K-1!) * x^K * (1-x)^(L-K-1)
                       = (L choose K+1) * (K+1) * x^K * (1-x)^(L-K-1)

Beta distribution's relation to Binomial distribution.
https://dlmf.nist.gov/8.17#E5. Formulas 8.17.4/5.

    pbeta(x, K+1, L-K)         = 1 - pbinom(K, L, x)
    pbeta(x, K+1, L-K)         = survbinom(K, L, x)
    pbeta(x, K+1, L-K)         = 1 - pbeta(1 - x, L - K, K + 1)
    pbeta(1 - x, L - K, K + 1) = pbinom(K, L, x)
    pbeta(x, L - K, K + 1)     = pbinom(K, L, 1 - x)
    dbeta(x, K + 1, L - K)     = (L - K) / (1 - x) * dbinom(K, L, x)
    dbeta(x, K + 1, L- K + 1)  = (L + 1) * dbinom(K, L, x)

    Binomial mass probabilities have recurrence relation:
    dbinom(K+1, L, p) = dbinom(K, L, p) * p / (1 - p) * (L - K) / (K + 1)

    Check in R:
    K <- 10
    L <- 100
    x <- 0.1
    p <- 0.1

    dbeta(x, K + 1, L - K)
    choose(L, K + 1) * (K + 1) * x^K * (1 - x)^(L - K - 1)
    (L - K) / (1 - x) * dbinom(K, L, x)

    pbeta(x, K + 1, L - K)
    1 - pbeta(1 - x, L - K, K + 1)
    1 - pbinom(K, L, x)
    survbinom(K, L, x)

    pbeta(x, L - K, K + 1)
    pbinom(K, L, 1 - x)

    dbinom(K+1, L, p)
    dbinom(K, L, p) * p / (1 - p) * (L - K) / (K+1)

    For larger K and L Beta and Gamma distributions converge to
    normal distribution with Beta and Gamma mean and SD. Try
    functions plot.BetaDist(K, L) and plot.GammaDist(K).

*/

double riemann(double (*f)(double), double a, double h, double tol);
double simpson(double (*f)(double), double a, double h, double tol);
double adaSimpson(double (*f)(double), double a, double b, double fa,
                  double fm, double fb, double Iprev, double abstol,
                  double reltol, int depth);
void quickUniSelect(int k, NumericVector p);
static double fisher(NumericVector p);
static double sidak(NumericVector p);
double betaMean(double a, double b);
double betaSD(double a, double b);
double betaSkewness(double a, double b);

// Benchmark baseline function
// [[Rcpp::export]]
double baseNull(double x) {
    return x;
}

// Global "constants" in integration
static double K;     // rank, number of smallest
static double L;     // number of p-values
static double LW;    // log(p1 x ... x pK), test statistic
static double LBETA; // log(L! / (K! * L-K-1!))
static double LKF;   // log((K-1)!)
static double PCUT;
static int ERR = 0; // yet not used
static int HIT = 0;
static int MISS = 0;

#define OK 2

// betaCutPoint returns limit for which
// 1 - pbeta(x, K+1, L-K) < ~1e-12, when x > limit.
// Positive left skewness lifts the right tail off zero.
//
// [[Rcpp::export]]
double
betaCutPoint(double k, double l) {
    double pcut, dist;
    dist = 7 + fmax(0, 9 * betaSkewness(k + 1, l - k));
    pcut = betaMean(k + 1, l - k) + dist * betaSD(k + 1, l - k);
    return fmin(1, pcut);
}

// rtpStat calculates RPT method test statistic from p-value vector p.
// U1 x U2 x ... x Un ~ e^-n (not 2^-n) for big n and Ui ~ Unif(0, 1).
static double rtpStatistic(int k, NumericVector p) {
    quickUniSelect(k, p);
    double lw = 0;
    for (int i = 0; i < k; i++)
        lw += log(p[i]);
    return lw;
}

// [[Rcpp::export]]
double init(double k, NumericVector p, int integrand = 0) {
    L = p.size();
    K = k;
    ERR = 0;
    if (L < K || L < 1) return ERR = -1;
    if (K == 1) return sidak(p);
    if (K == L) return fisher(p);

    LW = rtpStatistic(k, p);
    if (integrand == 1)
        LBETA = R::lbeta(K + 1, L - K);
    else {
        LKF = lgamma(K);
        PCUT = betaCutPoint(K, L);
    }
    return OK;
}

// [[Rcpp::export]]
double hitMiss(int reset) {
    double r = (double)HIT / (HIT + MISS);
    if (reset > 0) {
        HIT = 0;
        MISS = 0;
    }
    return r;
}

// Beta standard deviation.
// betaSD(K+1, L-K) ~ sqrt(K) / L, for small K/L.
// [[Rcpp::export]]
double betaSD(double a, double b) {
    double c = a + b;
    return sqrt(a * b / (c * c * (c + 1)));
}

// betaMode(K+1, L-K) = K / (L-1).
inline static double betaMode(double a, double b) {
    return (a - 1) / (a + b - 2);
}

// betaMean(K+1, L-K) = (K+1) / (L+1).
double betaMean(double a, double b) {
    return a / (a + b);
}

// [[Rcpp::export]]
double betaSkewness(double a, double b) {
    return 2 * (b - a) * sqrt(a + b + 1) / ((a + b + 2) * sqrt(a * b));
}

// sumRight sums bin(n, p) probabilities from a to b. prob must be dbinom(a, n, p).
inline static double sumRight(double a, double b, double n, double p, double prob) {
    const double reltol = 1e-12;
    double sum = prob;
    for (double j = a + 1; j <= b; j++) {
        prob *= p / (1 - p) * (n - j + 1) / j;
        sum += prob;
        if (prob <= sum * reltol) break;
    }
    return sum;
}

// sumLeft sums from right to left, a > b.
inline static double sumLeft(double a, double b, double n, double p, double prob) {
    const double reltol = 1e-12;
    double sum = prob;
    for (double j = a; j > b; j--) {
        prob *= (1 - p) / p * j / (n - j + 1);
        sum += prob;
        if (prob <= sum * reltol) break;
    }
    return sum;
}

// Binomial density / mass probability function (for survbinom).
inline static double dbinom(double k, double n, double p) {
    double lg = lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
    return exp(lg + k * log(p) + (n - k) * log(1 - p));
}

// Faster binomial survival / beta CDF function.
// survbinom(K, L, b) ~ R::pbeta(b, K + 1, L - K, 1, 0).
// Absolute difference of the functions is < 1e-12 near 1
// and < 1e-13 elsewhere. Relative difference is always < 1e-11.
// [[Rcpp::export]]
double survbinom(double k, double n, double p, double pcut) {
    const double minNormal = 0x1p-1022; // min normal double
    double prob;

    if (p >= pcut) return 1;
    if (p <= 0) return 0;

    if (k < n * p + 3 * sqrt(n * p * (1 - p))) {
        prob = exp(n * log(1 - p));
        if (prob > minNormal) {
            return 1 - sumRight(0, k, n, p, prob);
        }
        prob = dbinom(k, n, p);
        if (prob < minNormal)
            return 1;
        return 1 - sumLeft(k, 0, n, p, prob);
    }
    prob = exp(n * log(p));
    if (prob > minNormal)
        return sumLeft(n, k + 1, n, p, prob);

    prob = dbinom(k + 1, n, p);
    if (prob < minNormal)
        return prob;
    return sumRight(k + 1, n, n, p, prob);
}

// Beta density function.
inline static double dbeta(double x, double a, double b, double lbeta) {
    if (x <= 0 || x >= 1) return 0;

    return exp(-lbeta + (a - 1) * log(x) + (b - 1) * log(1 - x));
}

inline static double gammaSD(double k) {
    return sqrt(k);
}

inline static double gammaMode(double k) {
    return k - 1;
}

inline static double gammaMean(double k) {
    return k;
}

// Fast gamma survival function for (small) positive integers k = shape and rate 1.
// survgamma(g, K) = 1 - pgamma(g, K)
// [[Rcpp::export]]
double survgamma(double g, double k) {
    if (g <= 0) return 1;

    double p = exp(-g);
    if (p < 3e-308)
        return R::pgamma(g, k, 1, 0, 0); // right tail

    double s = p;
    for (double j = 1; j < k; j++) {
        p *= g / j;
        s += p;
    }
    return s;
}

// Gamma density function (g,k) = e^-g * g^(k-1) / (k-1)!
inline static double dgamma(double g, double k, double lkf) {
    if (g <= 0) return 0;

    return exp((k - 1) * log(g) - g - lkf);
}

// fBetaD is rtp integrand Beta PDF x (1 - Gamma CDF) over [0, 1].
// This implements the equations in Dudbridge and Koeleman.
// [[Rcpp::export]]
double fBetaD(double b) {
    if (b <= 0 || b >= 1) return 0;

    double g = K * log(b) - LW;
    return dbeta(b, K + 1, L - K, LBETA) * survgamma(g, K);
    // return R::dbeta(b, K + 1, L - K, 0) * R::pgamma(g, K, 1, 0, 0);
}

// fGammaD is rtp integrand Gamma PDF x Beta CDF over [0, inf).
// [[Rcpp::export]]
double fGammaD(double g) {
    if (g <= 0) return 0;

    double b = exp((g + LW) / K);
    return dgamma(g, K, LKF) * survbinom(K, L, b, PCUT);
    // return R::dgamma(g, K, 1, 0) * R::pbeta(b, K + 1, L - K, 1, 0);
}

// fBetaQ is rtp integrand over Beta probabilities in [0, 1].
// This is the integrand in Vsevolozhskaya et al.
// [[Rcpp::export]]
double fBetaQ(double p) {
    if (p <= 0) return 1;

    double b = R::qbeta(p, K + 1, L - K, 1, 0); // Inverse beta CDF method
    double g = K * log(b) - LW;
    return survgamma(g, K);
    // return R::pgamma(g, K, 1, 0, 0);
}

// fGammaQ is rtp integrand over Gamma probabilities in [0, 1].
// fGammaQ(1 - u) is very close to fBetaQ(u). ~2 x faster than fBetaQ.
// [[Rcpp::export]]
double fGammaQ(double p) {
    if (p >= 1) return 1;

    double g = R::qgamma(p, K, 1, 1, 0); // Inverse gamma CDF method
    double b = exp((g + LW) / K);
    return survbinom(K, L, b, PCUT);
    // return R::pbeta(b, K + 1, L - K, 1, 0);
}

// fBetaDtop approximates the location of highest point of fBetaD.
// This is a manually fitted from hat model, good enough anyway.
// Maximum from equations has no closed form solution?
// [[Rcpp::export]]
double fBetaDtop() {
    double right = betaMode(K + 1, L - K);
    double left = exp(LW / K) * 1.5;
    if (left > right) return right;

    double weight = 0.2 + 3 * left / right + 2 * right;
    return (left + weight * right) / (1 + weight);
}

// rtpDbetaRiema integrates fBetaD from 0 to 1 by Riemann sum integral.
// [[Rcpp::export]]
double rtpDbetaRiema(double k, NumericVector p, double tol = 1e-10, double stepscale = 1) {
    const double cStep = 0.75;
    double h, s, l;

    if ((s = init(k, p, 1)) != OK) return s;

    l = p.size();
    h = cStep * betaSD(k + 1, l - k) * stepscale;
    if (k < 8) h *= k / 8;

    return riemann(&fBetaD, 0, h, tol); // x >= 1 -> fBetaD(x) = 0.
}

// rtpDbetaAsimp integrates fBetaD from 0 to 1 by adaptive Simpson's 1/3 rule.
// [[Rcpp::export]]
double rtpDbetaAsimp(double k, NumericVector p, double abstol = 1e-7, double reltol = 1e-3) {
    const double depth = 25;
    double top, fa, fm, fb, I, s;

    if ((s = init(k, p, 1)) != OK) return s;

    top = fBetaDtop();
    fa = 0;
    fm = fBetaD(top / 2);
    fb = fBetaD(top);
    I = adaSimpson(&fBetaD, 0, top, fa, fm, fb, 2, abstol / 2, reltol, depth);

    fa = fb;
    fm = fBetaD((top + 1) / 2);
    fb = 0;
    I += adaSimpson(&fBetaD, top, 1, fa, fm, fb, 2, abstol / 2, reltol, depth);
    return I;
}

// rtpDgammaRiema integrates fGammaD from 0 to inf by Riemann sum integral.
// For k = 10 this needs < 10 integrand function evaluations.
// [[Rcpp::export]]
double rtpDgammaRiema(double k, NumericVector p, double tol = 1e-10, double stepscale = 1) {
    const double cStep = 1.25;
    double h, s;
    if ((s = init(k, p, 0)) != OK) return s;

    h = cStep * gammaSD(k) * stepscale;

    return riemann(&fGammaD, 0, h, tol);
}

// rtpDgammaSimp integrates fGammaD from 0 to inf by fixed step Simpson's 1/3 rule.
// [[Rcpp::export]]
double rtpDgammaSimp(double k, NumericVector p, double tol = 1e-10, double stepscale = 1) {
    const double cStep = 1.25;
    double h, s;
    if ((s = init(k, p, 0)) != OK) return s;

    h = cStep * gammaSD(k) * stepscale;

    return simpson(&fGammaD, 0, h, tol);
}

// rtpRiema uses rtpDgammaRiema for K < L / 6 and
// rtpDbetaRiema for bigger K's. This gives good
// all over accuracy.
// [[Rcpp::export]]
double rtpRiema(double k, NumericVector p, double tol = 1e-10, double stepscale = 1) {
    double l = p.size();
    if (k / l > 1.0 / 3.0)
        return rtpDbetaRiema(k, p, tol, stepscale);
    return rtpDgammaRiema(k, p, tol, stepscale);
}

// sidak returns the probability of getting one or
// more p-values = minimum observed p-value.
// RTP for K == 1.
static double sidak(NumericVector p) {
    int l = p.size();
    double pmin = p[0];

    for (int i = 1; i < l; i++)
        if (pmin > p[i]) pmin = p[i];

    return R::pbinom(0, l, pmin, 0, 0); // rigth tail
    // return 1 - pow(1 - pmin, l);
}

// Standard Fisher's method using all p-values.
// Solved by Gamma distribution. RTP for K == L.
static double fisher(NumericVector p) {
    int l = p.size();
    double lw = 0;
    for (int i = 0; i < l; i++)
        lw += log(p[i]);
    return R::pgamma(-lw, l, 1, 0, 0);
}