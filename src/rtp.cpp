

/* See referencies in the README.md.
https://mran.microsoft.com/snapshot/2020-05-03/web/packages/metap/metap.pdf

The K+1'th smallest of L unif(0, 1) numbers is distributed Beta(K+1, L-K).
The density function of order statistic x ~ Beta(K+1, L-K) is
    dbeta(x, K+1, L-K) = L! / (K! * L-K-1!) * x^K * (1-x)^(L-K-1)
                       = (L choose K+1) * (K+1) * x^K * (1-x)^(L-K-1)

Beta distribution's relation to Binomial distribution.
https://dlmf.nist.gov/8.17#E5. Formulas 8.17.4/5.

    pbeta(x, K+1, L-K)         = 1 - pbinom(K, L, x)
    pbeta(x, K+1, L-K)         = survbinom(K, L, x)
    pbeta(x, K, L)             = survbinom(K-1, L+K-1, x)
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

*/

#include <Rcpp.h>
using namespace Rcpp;

#include "fundec.h"

// Global "constants" in integration
static double K;     // rank, number of smallest numbers
static double L;     // number of p-values
static double LW;    // log(p1 x ... x pK), test statistic
static double LBETA; // log(L! / (K! * L-K-1!))
static double LGAM;  // log((K-1)!)

// Benchmark baseline function
// [[Rcpp::export]]
double
baseNull(double x) {
    return x;
}

// betaCutPoint returns limit for which
// 1 - pbeta(x, K+1, L-K) < ~1e-12, when x > limit.
// Leaning left (positive left skewness) lifts the right leg/tail off the
// ground and we must go further right to get near zero.
//
// // [[Rcpp::export]]
// double betaCutPoint(double k, double l) {
//     double dist = 7 + fmax(0, 9 * betaSkewness(k + 1, l - k));
//     double pcut = betaMean(k + 1, l - k) + dist * betaSD(k + 1, l - k);
//     return fmin(1, pcut);
// }

// statRTP returns RPT method test statistic from p-value vector p.
// [[Rcpp::export(name = "stat.rtp")]]
double statRTP(long k, NumericVector q) {

    NumericVector p = clone(q); // we don't want to modify p
    kSelect(k, p);
    double lw = 0;
    for (long i = 0; i < k; i++)
        lw += -log(p[i]);
    return lw;
}

// [[Rcpp::export]]
double init(double k, NumericVector p, int integrand = 0) {
    L = p.size();
    K = k;
    if (L < K || L < 1) return ERR;
    if (K == 1 && L == 1) return p[1];
    if (K == 1) return single(p);
    if (K == L) return fisher(p);
    if (integrand == 0) return OK;

    LW = statRTP(K, p);
    if (integrand == 1) {
        LBETA = R::lbeta(K + 1, L - K);
    } else if (integrand == 2) {
        LGAM = lgamma(K);
    }
    return OK;
}

// Beta standard deviation.
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

inline double betaSkewness(double a, double b) {
    return 2 * (b - a) * sqrt(a + b + 1) / ((a + b + 2) * sqrt(a * b));
}

// sumRight sums bin(n, p) probabilities from a to b. prob is dbinom(a, n, p).
static double sumRight(double a, double b, double n, double p, double prob) {
    const double reltol = 1e-12;
    double sum = prob;
    for (double j = a + 1; j <= b; j++) {
        prob *= p / (1 - p) * (n - j + 1) / j;
        sum += prob;
        if (prob < sum * reltol) break;
    }
    return sum;
}

// sumLeft sums from right to left, a > b.
static double sumLeft(double a, double b, double n, double p, double prob) {
    const double reltol = 1e-12;
    double sum = prob;
    for (double j = a; j > b; j--) {
        prob *= (1 - p) / p * j / (n - j + 1);
        sum += prob;
        if (prob < sum * reltol) break;
    }
    return sum;
}

// Binomial density / mass probability function.
inline static double dbinom(double k, double n, double p) {
    double lg = lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
    return exp(lg + k * log(p) + (n - k) * log(1 - p));
}

// Faster Binomial survival function.
// survbinom(K, L, p) ~ pbinom(K, L, p, lower.tail = F).
// Absolute difference of the functions is < 5e-14 and
// relative difference < 5e-12 for k <= 100, n <= 2000,
// n > 1.5 x k, values > 1e-300.
// Arithmetic with subnormal doubles is avoided, because of slow
// calculation and large errors produced.
// Summing cumulative probability is done from small to big values,
// if possible (smallest > minNormal).
// [[Rcpp::export]]
double survbinom(double k, double n, double p) {
    double prob;

    if (p >= 1) return 1;
    if (p <= 0) return 0;

    double mean = n * p;
    double sd = sqrt(n * p * (1 - p));

    if (k < mean + 2.5 * sd) { // left tail sum, sum << 1
        prob = exp(n * log(1 - p));
        if (prob > DBL_MIN) { // normal exit
            return 1 - sumRight(0, k, n, p, prob);
        }
        prob = dbinom(k, n, p);
        if (prob < DBL_MIN) {
            return 1;
        }
        return 1 - sumLeft(k, 0, n, p, prob);
    }
    // right tail sum
    prob = exp(n * log(p));
    if (prob > DBL_MIN) {
        return sumLeft(n, k + 1, n, p, prob);
    }
    prob = dbinom(k + 1, n, p);
    if (prob < DBL_MIN) {
        return prob;
    }
    return sumRight(k + 1, n, n, p, prob);
}

// Faster Beta distribution function by survbinom.
inline static double pbeta(double x, double a, double b) {

    return survbinom(a - 1, a + b - 1, x);
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

// Fast gamma survival function for Gamma(K, 1) and not very big k.
// survgamma(g, K) = 1 - pgamma(g, K)
// [[Rcpp::export]]
double survgamma(double g, double k) {
    if (g <= 0) return 1;
    double p = exp(-g);
    if (p < DBL_MIN)
        return R::pgamma(g, k, 1, 0, 0); // right tail
    double s = p;
    for (double j = 1; j < k; j++) {
        p *= g / j;
        s += p;
    }
    return s;
}

// Gamma density function for Gamma(k, 1).
inline static double dgamma(double g, double k, double lkf) {
    if (g <= 0) return 0;

    return exp((k - 1) * log(g) - g - lkf);
}

// fBetaD is rtp integrand Beta PDF x (1 - Gamma CDF) over [0, 1].
// This implements the equations in Dudbridge and Koeleman.
// [[Rcpp::export]]
double fBetaD(double b) {
    if (b <= 0 || b >= 1) return 0;

    double g = K * log(b) + LW;
    return dbeta(b, K + 1, L - K, LBETA) * survgamma(g, K);
    // return R::dbeta(b, K + 1, L - K, 0) * R::pgamma(g, K, 1, 0, 0);
}

// fGammaD is rtp integrand Gamma PDF x Beta CDF over [0, inf).
// [[Rcpp::export]]
double fGammaD(double g) {
    if (g <= 0) return 0;
    double x = dgamma(g, K, LGAM);
    if (g - LW < 0) {
        x *= pbeta(exp((g - LW) / K), K + 1, L - K);
    }
    return x;
    // return R::dgamma(g, K, 1, 0) * R::pbeta((g - LW) / K), K + 1, L - K, 1, 0);
}

// fBetaDtop approximates the location of highest point of fBetaD.
// [[Rcpp::export]]
double fBetaDtop() {
    double right = betaMode(K + 1, L - K);
    double left = exp(-LW / K) * 1.5;
    if (left > right) return right;

    double weight = 0.2 + 3 * left / right + 2 * right;
    return (left + weight * right) / (1 + weight);
}

// [[Rcpp::export(name = "p.rtp.dbeta.lw")]]
double rtpDbetaLW(double lw, double k, double l, double tol = 1e-12,
                  double stepscale = 1) {
    const double cStep = 0.5;
    K = k;
    L = l;
    LBETA = R::lbeta(k + 1, l - k);
    LW = lw;
    double h = cStep * betaSD(k + 1, l - k) * stepscale;
    if (k < 8) h *= k / 8;
    return riemann(&fBetaD, 0, h, tol);
}

// rtpDbeta integrates fBetaD from 0 to 1 by Riemann sum integral.
// [[Rcpp::export(name = "p.rtp.dbeta")]]
double rtpDbeta(double k, NumericVector p, double tol = 1e-10,
                double stepscale = 1) {
    double s = init(k, p, 0);
    if (s != OK) return s;
    double lw = statRTP(K, p);
    return rtpDbetaLW(lw, K, L, tol, stepscale);
}

// rtpDbetaAsimp integrates fBetaD from 0 to 1 by adaptive Simpson's 1/3 rule.
// [[Rcpp::export(name = "p.rtp.dbeta.asimp")]]
double rtpDbetaAsimp(double k, NumericVector p, double tol = 1e-4) {
    const double depth = 25;
    double top, fa, fm, fb, I, s;

    if (tol < 1e-12) tol = 1e-12;
    if (tol > 1e-2) tol = 1e-2;

    if ((s = init(k, p, 1)) != OK) return s;

    top = K / L * 0.5; // very approx.
    fa = 0;
    fm = fBetaD(top / 2);
    fb = fBetaD(top);
    I = adaSimpson(&fBetaD, 0, top, fa, fm, fb, 2, tol, depth);

    fa = fb;
    fm = fBetaD((top + 1) / 2);
    fb = 0;
    I += adaSimpson(&fBetaD, top, 1, fa, fm, fb, 2, tol, depth);
    return fmin(1, I);
}

// [[Rcpp::export(name = "p.rtp.dgamma.lw")]]
double rtpDgammaLW(double lw, double k, double l,
                   double tol = 1e-10, double stepscale = 1) {
    const double cStep = 1.25;
    K = k;
    L = l;
    LGAM = lgamma(K);
    LW = lw;

    double h = cStep * gammaSD(k) * stepscale;
    return riemann(&fGammaD, 0, h, tol);
}

// rtpDgamma integrates fGammaD from 0 to inf by Riemann sum integral.
// [[Rcpp::export(name = "p.rtp.dgamma")]]
double rtpDgamma(double k, NumericVector p, double tol = 1e-10,
                 double stepscale = 1) {
    double s = init(k, p, 0);
    if (s != OK) return s;
    double l = p.size();

    double lw = statRTP(k, p);
    return rtpDgammaLW(lw, k, l, tol, stepscale);
}

// rtpLW uses rtpDgammaLW for K/L <= 1/3 and rtpDbetaLW for K/L > 1/3.
// [[Rcpp::export(name = "p.rtp.lw")]]
double rtpLW(double lw, double k, double l, double tol = 1e-14,
             double stepscale = 0.5) {
    if (k / l > 1.0 / 3.0)
        return rtpDbetaLW(lw, k, l, tol, stepscale);
    return rtpDgammaLW(lw, k, l, tol, stepscale);
}

// rtp uses rtpDgamma for K/L <= 1/3 and rtpDbeta for K/L > 1/3.
// [[Rcpp::export(name = "p.rtp")]]
double rtp(double k, NumericVector p, double tol = 1e-14,
           double stepscale = 0.5) {
    double s = init(k, p, 0);
    if (s != OK) return s;

    double lw = statRTP(K, p);
    return rtpLW(lw, K, L, tol, stepscale);
}

// [[Rcpp::export(name = "p.rtp.vector")]]
NumericVector rtpLWvec(NumericVector w, double k, double l,
                       double tol = 1e-10, double stepscale = 0.5) {
    long s = w.size();
    NumericVector v(s);
    for (long i = 0; i < s; i++) {
        v[i] = rtpLW(w[i], k, l, tol, stepscale);
    }
    return v;
}

// single returns the probability of getting one or
// more p-values <= minimum observed p-value.
double single(NumericVector p) {
    long l = p.size();
    if (l == 1) return p[0];
    double pmin = p[0];

    for (long i = 1; i < l; i++)
        if (pmin > p[i]) pmin = p[i];

    return R::pbinom(0, l, pmin, 0, 0); // rigth tail
    // return 1 - pow(1 - pmin, l); // loses accuracy
}

// Standard Fisher's method using all p-values.
// Solved by Gamma distribution. RTP for K == L.
double fisher(NumericVector p) {
    const double e = 2.71828;
    long n = p.size();
    double lw = 0, w = 1;
    for (long i = 0; i < n; i++)
        w *= p[i] * e; // mean(w) = 1
    lw = log(w) - n;
    if (w == 0 || w > DBL_MAX) {
        lw = 0;
        for (long i = 0; i < n; i++)
            lw += log(p[i]);
    }
    return R::pgamma(-lw, n, 1, 0, 0);
}

static double densityIntegrand(double b) {
    if (b <= 0 || b >= 1) return 0;
    double g = K * log(b) + LW;
    return dbeta(b, K + 1, L - K, LBETA) * dgamma(g, K, LGAM);
}

// Vectorized RTP density function.
// [[Rcpp::export(name = "drtp")]]
NumericVector drtp(NumericVector x, double k, double l, double tol = 1e-16,
                   double stepscale = 1) {
    const double cStep = 0.05;
    long n = x.size();
    NumericVector v(n);
    K = k;
    L = l;
    LGAM = lgamma(k);
    LBETA = R::lbeta(k + 1, l - k);
    double h = cStep * k / l * stepscale;
    for (long i = 0; i < n; i++) {
        LW = x[i];
        v[i] = riemann(&densityIntegrand, 0, h, tol);
    }
    return v;
}

// Vectorized logarithmic RTP random variable generator.
// [[Rcpp::export(name = "rrtp")]]
NumericVector rrtp(long n, double k, double l) {
    NumericVector v(n);

    for (long i = 0; i < n; i++)
        v[i] = k * log(R::rbeta(k + 1, l - k)) - R::rgamma(k, 1);

    return v;
}

// rtpSimulated solves rtp p-value by Monte Carlo simulation.
// lw = -log(p1 x ... pk) is the test statistic from the given p-values.
// The k+1's smallest of L Unif(0,1) numbers is simulated by
// one Beta(K+1, l-K) random number.
// The logarithm of product U1 x ... Uk is simulated by a single
// Gamma(K, 1) random number.
//
// [[Rcpp::export]]
double rtpSimulated(double k, NumericVector p, uint64_t samples) {
    double l = p.size(), more = 0, lw, lwr;
    if (k < 1 || k >= l || l < 2) return ERR;
    lw = statRTP(k, p);
    for (uint64_t i = 0; i < samples; i++) {
        lwr = -k * log(R::rbeta(k + 1, l - k)); // Simulated k+1'th smallest of l numbers
        lwr += R::rgamma(k, 1);                 // Simulated random product
        if (lwr > lw) more++;
    }
    return more / samples;
}
