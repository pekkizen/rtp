
#include <Rcpp.h>
using namespace Rcpp;

// "rtp.h"
double riemann(double (*f)(double), double a, double h, double tol, double SD);
double simpson(double (*f)(double), double a, double h, double tol, double SD,
               double hlim, double hmul);
double adaSimpson(double (*f)(double), double a, double b, double fa,
                  double fm, double fb, double Iprev, double abstol,
                  double reltol, int depth);
void uniSelect(int k, NumericVector p);

/*
https://en.wikipedia.org/wiki/Gamma_distribution
Characterization using shape α and rate β
https://en.wikipedia.org/wiki/Beta_distribution
Probability density function
https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function
https://en.wikipedia.org/wiki/Order_statistic
Order statistics sampled from a uniform distribution

Density function of the K + 1'th Uniform(0, 1) order statistic is
dbeta(x, K+1, L-K) = L! / (K! * L-K-1!) * x^K * (1-x)^(L-K-1)
                   = (L choose K+1) * (K+1) * x^K * (1-x)^(L-K-1)
Beta distribution's relation to Binomial distribution.
pbeta(x, K+1, L-K) = 1 - pbinom(K, L, x)
dbeta(x, K+1, L-K) = (L-K) / (1-x) * dbinom(K, L, x)
*/

// [[Rcpp::export]]
double baseNull(double x) {
    return x;
}

// Global "constants" in integration
static double K;     // rank, number of smallest
static double L;     // number of p-values
static double LBETA; // lbeta(K + 1, L - K) = log(K! * L-K-1! / L!)
static double LKF;   // lgamma(K) = log((K - 1)!)
static double LW;    // log(p1 x ... x pK), test statistic
static int ERR = -1;

// [[Rcpp::export]]
double init(int k, NumericVector p) {
    L = p.size();
    K = k;
    if (L < K || L < 1) return ERR;
    uniSelect(K, p);
    LW = 0;
    for (int i = 0; i < K; i++)
        LW += log(p[i]);

    LBETA = R::lbeta(K + 1, L - K);
    LKF = lgamma(K);
    return L;
}

// Beta standard deviation.
// betaSD(K+1, L-K) = ~sqrt(K) / L, for small K/L.
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
inline static double betaMean(double a, double b) {
    return a / (a + b);
}

// Fast binomial survival function for small K.
// sbinom(K, L, b) = pbeta(b, K+1, L-K) = 1 - pbinom(K, L, b)
static double sbinom(double k, double n, double p) {
    if (p <= 0) return 0;
    if (p >= 1) return 1;

    double prob = R::dbinom(0, n, p, 0);
    double cdf = prob;

    for (double j = 1.0; j <= k; j++) {
        prob *= p * (n + 1 - j) / ((1 - p) * j);
        cdf += prob;
    }
    if (cdf > (1.0 - 1e-10) || cdf == 0)
        return R::pbinom(k, n, p, 0, 0);

    return 1 - cdf;
}

// Beta density function.
// lbeta is lbeta(a, b).
inline static double dbeta(double x, double lbeta, double a, double b) {
    if (x <= 0 || x >= 1) return 0;

    return exp(-lbeta + (a - 1) * log(x) + (b - 1) * log(1 - x));
}

// Hight of Beta density function. For plot. Function value at mode.
// [[Rcpp::export]]
double dbetaHight(double a, double b) {
    return dbeta(betaMode(a, b), R::lbeta(a, b), a, b);
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

// Fast gamma survival function for small positive
// integers k=shape and rate 1.
static double sgamma(double g, double k) {
    if (g <= 0) return 1;

    if (g > 700 || k > 100)
        return R::pgamma(g, k, 1, 0, 0); // right tail

    double p = exp(-g); // p > 0
    double s = p;
    for (double j = 1.0; j < k; j++) {
        p *= g / j;
        s += p;
    }
    return fmin(1, s);
}

// Gamma density function (g,k) = e^-g * g^(k-1) / (k-1)!
// lkf is logarithm of (k-1)!.
inline static double dgamma(double g, double lkf, double k) {
    if (g <= 0) return 0;

    return exp((k - 1) * log(g) - g - lkf);
}

// fBetaD is rtp integrand over Beta density in [0, 1].
// This implements the equations in Dudbridge and Koeleman.
// [[Rcpp::export]]
double fBetaD(double b) {
    if (b <= 0 || b >= 1) return 0;

    double g = K * log(b) - LW;
    return dbeta(b, LBETA, K + 1, L - K) * sgamma(g, K);
}

// fGammaD is rtp integrand over Gamma density in [0, inf).
// [[Rcpp::export]]
double fGammaD(double g) {
    if (g <= 0) return 0;

    double b = exp((g + LW) / K);
    return dgamma(g, LKF, K) * sbinom(K, L, b);
    // return dgamma(g, LKF, K) * R::pbeta(b, K + 1, L - K, 1, 0);
}

// fBetaQ is rtp integrand over Beta probabilities in [0, 1].
// This is the integrand in Vsevolozhskaya et al.
// [[Rcpp::export]]
double fBetaQ(double p) {
    if (p <= 0) return 1;

    double b = R::qbeta(p, K + 1, L - K, 1, 0); // Inverse beta CDF method
    double g = K * log(b) - LW;
    return sgamma(g, K);
}

// fGammaQ is rtp integrand over Gamma probabilities in [0, 1].
// fGammaQ(1 - u) is very close to fBetaQ(u). ~2 x faster than fBetaQ.
// [[Rcpp::export]]
double fGammaQ(double p) {
    if (p >= 1) return 1;

    double g = R::qgamma(p, K, 1, 1, 0); // Inverse gamma CDF method
    double b = exp((g + LW) / K);
    return sbinom(K, L, b);
    // return R::pbeta(b, K + 1, L - K, 1, 0);
}

// fBetaDtop approximates the location of highest point of fBetaD.
// This is manually fitted from hat model.
// Maximum from equations has no closed form solution?
// [[Rcpp::export]]
double fBetaDtop() {

    double right = betaMode(K + 1, L - K);
    double left = exp(LW / K) * 1.5;
    if (left > right) return right;

    double weight = 0.2 + 3 * left / right + 2 * right;
    return (left + weight * right) / (1 + weight);
}

// Smallest p-value "method". K == 1.
static double probSmallest(NumericVector p) {
    double l = p.size();
    double pmin = p[0];

    for (int i = 1; i < l; i++)
        if (pmin > p[i]) pmin = p[i];
    return R::pbinom(0, l, pmin, 0, 0); // rigth tail
}

// Standard Fisher's method using all p-values.
static double fisher(NumericVector p) {
    double l = p.size();
    double lw = 0;
    for (int i = 0; i < l; i++)
        lw += log(p[i]);
    return R::pgamma(-lw, l, 1, 0, 0);
}

// rtpDbetaRiema integrates fBetaD from 0 to 1 by Riemann sum integral.
// [[Rcpp::export]]
double rtpDbetaRiema(int k, NumericVector p, double tol = 1e-10, double stepscale = 1) {
    const double steps = 8.0, maxstep = 0.05;
    double betaTop, h, SD, l;

    if (k == 1) return probSmallest(p);
    if (k == p.size()) return fisher(p);
    if (init(k, p) < 0) return ERR;

    l = p.size();
    SD = betaSD(k + 1, l - k);
    betaTop = betaMode(k + 1, l - k);
    h = fmin(maxstep, fmin(6 * SD, betaTop) / steps);
    h *= stepscale;

    return riemann(&fBetaD, 0, h, tol, SD); // x >= 1 -> fBetaD(x) = 0.
}

// rtpDbetaAsimp integrates fBetaD from 0 to 1.
// [[Rcpp::export]]
double rtpDbetaAsimp(int k, NumericVector p, double abstol = 1e-7, double reltol = 1e-3) {
    double top, fa, fm, fb, I;

    if (k == 1) return probSmallest(p);
    if (k == p.size()) return fisher(p);
    if (init(k, p) < 0) return ERR;

    top = fBetaDtop();
    fa = 0;
    fm = fBetaD(top / 2);
    fb = fBetaD(top);
    I = adaSimpson(&fBetaD, 0, top, fa, fm, fb, 2, abstol / 2, reltol, 25);

    fa = fb;
    fm = fBetaD((top + 1) / 2);
    fb = 0;
    I += adaSimpson(&fBetaD, top, 1, fa, fm, fb, 2, abstol / 2, reltol, 25);
    return I;
}

// rtpDgammaRiema integrates fGammaD from 0 to inf by Rieman sum integral.
// [[Rcpp::export]]
double rtpDgammaRiema(int k, NumericVector p, double tol = 1e-10, double stepscale = 1) {
    const double cStep = 1.25;
    double h, SD;

    if (k == 1) return probSmallest(p);
    if (k == p.size()) return fisher(p);
    if (init(k, p) < 0) return ERR;

    SD = gammaSD(k);
    h = cStep * SD;
    if (k < 6) h *= (double)k / 6;
    h *= stepscale;

    return riemann(&fGammaD, 0, h, tol, SD);
}

// rtpDgammaSimp integrates fGammaD from 0 to inf by fixed step Simpson's 1/3 rule.
// [[Rcpp::export]]
double rtpDgammaSimp(int k, NumericVector p, double tol = 1e-10, double stepscale = 1) {
    const double cStep = 1.5;
    double h, SD, hlim = 0.15, hmul = 1.5;

    if (k == 1) return probSmallest(p);
    if (k == p.size()) return fisher(p);
    if (init(k, p) < 0) return ERR;

    if (stepscale != 1) hlim = 0; // h, hlim and hmul are dependant
    SD = gammaSD(k);
    h = cStep * SD;
    if (k < 6) h *= (double)k / 6;
    h *= stepscale;

    return simpson(&fGammaD, 0, h, tol, SD, hlim, hmul);
}
