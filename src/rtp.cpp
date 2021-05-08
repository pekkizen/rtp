
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
    pbeta(x, L - K, K + 1)     = pbinom(K, L, 1-x)
    
    dbeta(x, K + 1, L - K) = (L - K) / (1 - x) * dbinom(K, L, x)
    Check in R:
    K <- 10
    L <- 100
    x <- 0.05
    dbeta(x, K + 1, L - K)
    choose(L, K + 1) * (K + 1) * x^K * (1 - x)^(L - K - 1)
    (L - K) / (1 - x) * dbinom(K, L, x)
    pbeta(x, K + 1, L - K)
    1 - pbeta(1 - x, L - K, K + 1)
    1 - pbinom(K, L, x)
    survbinom(K, L, x)

    

    For large k the gamma(k,1) distribution converges to 
    normal distribution with mean k and SD sqrt(k).
    For k > 10 distributions start to look quite alike.
    Test by plot.GammaNorm(K) function.
*/
double riemann(double (*f)(double), double a, double h, double tol);
double simpson(double (*f)(double), double a, double h, double tol);
double adaSimpson(double (*f)(double), double a, double b, double fa,
                  double fm, double fb, double Iprev, double abstol,
                  double reltol, int depth);
void selectUnif(int k, NumericVector p);
double fisher(NumericVector p);
static double sidak(NumericVector p);

// Benchmark baseline function
// [[Rcpp::export]]
double baseNull(double x) {
    return x;
}

// Global "constants" in integration
static double K;  // rank, number of smallest
static double L;  // number of p-values
static double LW; // log(p1 x ... x pK), test statistic
static int ERR = 0;

#define OK 2

// [[Rcpp::export]]
double init(int k, NumericVector p, int density = 0) {
    L = p.size();
    K = k;
    ERR = 0;
    if (L < K || L < 1) return ERR = -1;
    if (K == 1) return sidak(p);
    if (K == L) return fisher(p);

    selectUnif(K, p);
    LW = 0;
    for (int i = 0; i < K; i++)
        LW += log(p[i]);
    return OK;
}

// Beta standard deviation.
// betaSD(K+1, L-K) = ~sqrt(K) / L, for small K/L.
//
// [[Rcpp::export]]
double betaSD(double a, double b) {
    double c = a + b;
    return sqrt(a * b / (c * c * (c + 1)));
}

// betaMode(K+1, L-K) = K / (L-1).
double betaMode(double a, double b) {
    return (a - 1) / (a + b - 2);
}

// betaMean(K+1, L-K) = (K+1) / (L+1).
double betaMean(double a, double b) {
    return a / (a + b);
}

inline static double ldbinom(double k, double n, double p) {
    static double lbc = 0;
    if (lbc == 0) lbc = R::lchoose(n, k);

    return lbc + k * log(p) + (n - k) * log(1 - p);
}

// Binomial distribution from k+1 to n.
// Right tail, 1-CDF, survival function.
//
inline static double pbinomRT(double k, double n, double p) {
    return R::pbinom(k, n, p, 0, 0); // right tail
}

// Fast binomial survival function for not very big k.
// survbinom(K, L, b) = 1 - pbinom(K, L, b) = pbeta(b, K+1, L-K)
// https://dlmf.nist.gov/8.17#E5. Formula 8.17.5
//
// [[Rcpp::export]]
double survbinom(double k, double n, double p) {
    if (p >= 1) return 1;
    if (p <= 0) return 0;
    if (k > 100 || n > 10000)
        return pbinomRT(k, n, p);

    double prob = exp(n * log1p(-p)); // (1-p)^n
    if (prob == 0) {
        if (ldbinom(k, n, p) < -744) // dbinom(k, n, p) = 0
            return 1;
        return pbinomRT(k, n, p);
    }
    double cdf = prob;
    for (double j = 1; j <= k; j++) {
        prob *= p / (1 - p) * (n + 1 - j) / j;
        cdf += prob;
    }
    if (cdf > (1 - 1e-14))
        // gives also very small (<1e-16) probabilities
        return pbinomRT(k, n, p);

    return 1 - cdf;
}

// Beta density function.
//
inline static double dbeta(double x, double a, double b) {
    static double lb = 0;
    if (x <= 0 || x >= 1) return 0;
    if (lb == 0) lb = R::lbeta(a, b);

    return exp(-lb + (a - 1) * log(x) + (b - 1) * log(1 - x));
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
//
// [[Rcpp::export]]
double survgamma(double g, double k) {
    if (g <= 0) return 1;

    double p = exp(-g);
    if (p == 0)
        return R::pgamma(g, k, 1, 0, 0); // right tail
    double s = p;
    for (double j = 1; j < k; j++) {
        p *= g / j;
        s += p;
    }
    return s;
}

// Gamma density function (g,k) = e^-g * g^(k-1) / (k-1)!
//
inline static double dgamma(double g, double k) {
    static double lkf = 0;
    if (g <= 0) return 0;
    if (lkf == 0) lkf = lgamma(k);

    return exp((k - 1) * log(g) - g - lkf);
}

// fBetaD is rtp integrand Beta PDF x (1 - Gamma CDF) over [0, 1].
// This implements the equations in Dudbridge and Koeleman.
//
// [[Rcpp::export]]
double fBetaD(double b) {
    if (b <= 0 || b >= 1) return 0;

    double g = K * log(b) - LW;
    return dbeta(b, K + 1, L - K) * survgamma(g, K);
    // return dbeta(b, K + 1, L - K) * R::pgamma(g, K, 1, 0, 0);
}

// fGammaD is rtp integrand Gamma PDF x Beta CDF over [0, inf).
//
// [[Rcpp::export]]
double fGammaD(double g) {
    if (g <= 0) return 0;

    double b = exp((g + LW) / K);
    return dgamma(g, K) * survbinom(K, L, b);
    // return dgamma(g, K) * R::pbeta(b, K + 1, L - K, 1, 0);
}

// fBetaQ is rtp integrand over Beta probabilities in [0, 1].
// This is the integrand in Vsevolozhskaya et al.
//
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
//
// [[Rcpp::export]]
double fGammaQ(double p) {
    if (p >= 1) return 1;

    double g = R::qgamma(p, K, 1, 1, 0); // Inverse gamma CDF method
    double b = exp((g + LW) / K);
    return survbinom(K, L, b);
    // return R::pbeta(b, K + 1, L - K, 1, 0);
}

// fBetaDtop approximates the location of highest point of fBetaD.
// This is a manually fitted from hat model, quite good anyway.
// Maximum from equations has no closed form solution?
//
// [[Rcpp::export]]
double fBetaDtop() {
    double right = betaMode(K + 1, L - K);
    double left = exp(LW / K) * 1.5;
    if (left > right) return right;

    double weight = 0.2 + 3 * left / right + 2 * right;
    return (left + weight * right) / (1 + weight);
}

// rtpDbetaRiema integrates fBetaD from 0 to 1 by Riemann sum integral.
//
// [[Rcpp::export]]
double rtpDbetaRiema(int k, NumericVector p, double tol = 1e-12, double stepscale = 1) {
    const double cStep = 0.5;
    double h, s, l;

    if ((s = init(k, p, 1)) != OK) return s;

    l = p.size();
    h = cStep * betaSD(k + 1, l - k) * stepscale;
    if (k < 8) h *= (double)k / 8;

    return riemann(&fBetaD, 0, h, tol); // x >= 1 -> fBetaD(x) = 0.
}

// rtpDbetaAsimp integrates fBetaD from 0 to 1 by adaptive Simpson's 1/3 rule.
//
// [[Rcpp::export]]
double rtpDbetaAsimp(int k, NumericVector p, double abstol = 1e-7, double reltol = 1e-3) {
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

// rtpDgammaRiema integrates fGammaD from 0 to inf by Rieman sum integral.
// For k = 10 this needs < 10 integrand function evaluations.
//
// [[Rcpp::export]]
double rtpDgammaRiema(int k, NumericVector p, double tol = 1e-12, double stepscale = 1) {
    const double cStep = 1.25;
    double h, s;
    if ((s = init(k, p, 0)) != OK) return s;

    h = cStep * gammaSD(k) * stepscale;

    return riemann(&fGammaD, 0, h, tol);
}

// rtpDgammaSimp integrates fGammaD from 0 to inf by fixed step Simpson's 1/3 rule.
// [[Rcpp::export]]
double rtpDgammaSimp(int k, NumericVector p, double tol = 1e-10, double stepscale = 1) {
    const double cStep = 1.25;
    double h, s;
    if ((s = init(k, p, 0)) != OK) return s;

    h = cStep * gammaSD(k) * stepscale;

    return simpson(&fGammaD, 0, h, tol);
}
// double normalPDF(double x) {
//     return R::dnorm(x, 0, 1, 0);
// }
// // [[Rcpp::export]]
// double integnormal(double a, double h, double tol) {

//     // double h = 1;

//     return riemann(&normalPDF, a, h, 1e-10);
// }

// sidak returns the probability of getting one or
// more p-values = minimum observed p-value.
// RTP for K == 1.
//
static double sidak(NumericVector p) {
    int l = p.size();
    double pmin = p[0];

    for (int i = 1; i < l; i++)
        if (pmin > p[i]) pmin = p[i];
    // return R::pbinom(0, l, pmin, 0, 0); // rigth tail
    return 1 - pow(1 - pmin, l);
}

// Standard Fisher's method using all p-values.
// Solved by Gamma distribution. RTP for K == L.
//
// [[Rcpp::export]]
double fisher(NumericVector p) {
    int l = p.size();
    double lw = 0;
    for (int i = 0; i < l; i++)
        lw += log(p[i]);
    return R::pgamma(-lw, l, 1, 0, 0);
}