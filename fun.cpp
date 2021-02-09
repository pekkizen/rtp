#include <Rcpp.h>
using namespace Rcpp;

// Constants in integration
static double LK;
static double LG;
static double LW;
static double L;
static double K;

// For benchmark.
// [[Rcpp::export]]
static double fNull(double a) {
    return a;
}

// lbetaC returns logarithm of the Beta(a, b) factorial coefficient.
// lbetaC(k + 1, l - k) gives lgamma(l + 1) - lgamma(k + 1) - lgamma(l - k),
// which exponentiates to L! / (K! * L-K-1!) = (L choose K+1)*(K+1).
// See Dudbridge and Koeleman (2003).
static double lbetaC(double a, double b) {
    return lgamma(a + b) - lgamma(a) - lgamma(b);
    // return -R::lbeta(a, b);
}

// [[Rcpp::export]]
double init(double k, NumericVector p) {
    LW = 0;
    for (int i = 0; i < k; i++)
        LW += log(p[i]);

    K = k;
    L = p.length();
    LG = lbetaC(K + 1, L - K);
    LK = lgamma(K);
    return L;
}

// Beta standard deviation.
// [[Rcpp::export]]
static double betaSD(double a, double b) {
    double c = a + b;
    return sqrt(a * b / (c * c * (c + 1)));
}

inline static double betaMode(double a, double b) {
    return (a - 1) / (a + b - 2);
}

inline static double betaMean(double a, double b) {
    return a / (a + b);
}

// Beta quantile/inverse CDF function.
inline static double qbeta(double p, double a, double b) {
    return R::qbeta(p, a, b, 1, 0);
}

// Beta distribution function.
inline static double pbeta(double x, double a, double b) {
    if (x >= 1) return 1;
    return R::pbeta(x, a, b, 1, 0);
}

// Beta distribution function as binomial cdf.
// pbeta(x, k + 1, l - k) = 1 - pbinom(k, l, x)
// e.g. equation 8.17.5 in https://dlmf.nist.gov/8.17.

// Binomial distribution function (right tail).
// 1 - left tail loses much accuracy.
inline static double pbinomRT(double k, double l, double x) {
    if (x >= 1) return 1; // NaNs and inf loops otherwise
    return R::pbinom(k, l, x, 0, 0);
}

// Density of the k+1 th order statistic is
// l! / (k! * l-k-1!) * x^k * (1-x)^(l-k-1) =
// Beta density(x, k+1, l-k)

// Beta density function.
// lg is precalculated logarithm of l! / (k! * l-k-1!).
inline static double dbeta(double x, double lg, double a, double b) {
    if (x <= 0 || x >= 1) return 0;

    return exp(lg + (a - 1) * log(x) + (b - 1) * log(1 - x));
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

// Gamma distribution function (right tail).
inline static double pgammaRT(double g, double k) {
    return R::pgamma(g, k, 1, 0, 0);
}

// Gamma survival function.
// [[Rcpp::export]]
double gammaSurv(double g, double k) {
    double p, q, j;
    if (g <= 0) return 1;

    if (g > 700 || k > 100)
        return pgammaRT(g, k);

    q = exp(-g); // q > 0
    p = q;
    for (j = 1.0; j < k; j++) {
        q *= g / j;
        p += q;
    }
    return p;
}

// Gamma density function.
// lk is precalculated logarithm of (k-1)!.
inline static double dgamma(double g, double lk, double k) {
    if (g <= 0) return 0;

    return exp((k - 1) * log(g) - g - lk);
}

// Gamma quantile/inverse CDF function.
inline static double qgamma(double p, double k) {
    return R::qgamma(p, k, 1, 1, 0);
}

// fBeta is integrand over Beta density in [0, 1]. This is equivalent
// to the integrand equation in Dudbridge and Koeleman (2003).
// [[Rcpp::export]]
double fBeta(double b) {
    if (b <= 0 || b >= 1) return 0;

    double g = K * log(b) - LW;
    return dbeta(b, LG, K + 1, L - K) * gammaSurv(g, K);
}

// fGamma is integrand over Gamma density in [0, inf).
// [[Rcpp::export]]
static double fGamma(double g) {
    if (g <= 0) return 0;

    double b = exp((g + LW) / K);
    double dg = dgamma(g, LK, K);
    if (b >= 1) return dg;
    return dg * pbeta(b, K + 1, L - K);

    // Same results and speed:
    // return dg * pbinomRT(K, L, b);
}

// fBetaQuantile is integrand over Beta probabilities in [0, 1].
// This is the integrand in Vsevolozhskaya et al. (2019).
// [[Rcpp::export]]
static double fBetaQuantile(double p) {
    if (p <= 0) return 1;

    double b = qbeta(p, K + 1, L - K); // Inverse CDF method
    double g = K * log(b) - LW;
    return gammaSurv(g, K);
}

// fGammaQuantile is integrand over Gamma probabilities in[0, 1].
// fGammaQuantile(1 - u) is very close to fBetaQuantile(u).
// It is ~2 x faster than fBetaQuantile.
// [[Rcpp::export]]
double fGammaQuantile(double p) {
    if (p >= 1) return 1;

    double g = qgamma(p, K);
    double b = exp((g + LW) / K);
    return pbeta(b, K + 1, L - K);
}

// riemann integrates f from a to inf by Riemann sum. Integration
// stops when relative integral value on SD distance is < tol.
// Assumed unimodal f(x) -> 0 when x -> inf.
// Functions in f(x) must give a value, 0 or 1, for x > 1, not NaNs.
static double riemann(double (*f)(double), double a,
                      double h, double tol, double SD) {
    double fa, fsum = 0;

    tol /= SD;
    tol *= h;
    a += h / 2;
    while (true) {
        fa = f(a);
        fsum += fa;
        a += h;

        if (fa < fsum * tol) return h * fsum;
    }
}

// simpson integrates f from a to inf by Simpson's 1/3 rule.
// When f gets small, stepsize h is increased by hlim and hmul.
// Integration stops when relative integral value on SD distance is < tol.
// Assumed unimodal f(x) -> 0 when x -> inf.
static double simpson(double (*f)(double), double a, double h,
                      double tol, double SD, double hlim, double hmul) {
    double fa, fm, fb, fab, I = 0;

    tol /= SD;
    hlim /= SD;
    fb = f(a);
    while (true) {
        fa = fb;
        fm = f(a + h / 2);
        fb = f(a + h);
        fab = (fa + 4 * fm + fb) / 6;
        I += h * fab;
        a += h;

        if (fab < I * tol) return I;

        if (fab < I * hlim) h *= hmul;
    }
}

// adaSimpson integrates f from a to b by adaptive Simpson's 1/3 rule.
// Non standard exit condition uses also relative tolerance.
static double adaSimpson(double (*f)(double), double a, double b,
                         double fa, double fm, double fb, double Iprev,
                         double abstol, double reltol, int depth) {
    double h, fam, fmb, m, Ia, Ib, Iab, error;

    h = (b - a) / 4;
    fam = f(a + h);
    fmb = f(b - h);
    Ia = (fa + 4 * fam + fm) * (h / 3);
    Ib = (fm + 4 * fmb + fb) * (h / 3);

    Iab = Ia + Ib;
    error = (Iab - Iprev) / 15;
    if (depth <= 0 || fabs(error) < abstol || fabs(error) < Iab * reltol)
        return Iab + error;

    m = (a + b) / 2;
    Ia = adaSimpson(f, a, m, fa, fam, fm, Ia, abstol / 2, reltol, depth - 1);
    Ib = adaSimpson(f, m, b, fm, fmb, fb, Ib, abstol / 2, reltol, depth - 1);
    return Ia + Ib;
}

// riemannBeta integrates fBeta from 0 to 1 by Riemann sum integral.
// [[Rcpp::export]]
double riemannBeta(double k, NumericVector p, double tol = 1e-8, double stepscale = 1) {
    const double steps = 8.0, minstep = 0.05;
    double bTop, h, a, SD, l;

    l = init(k, p);
    SD = betaSD(k + 1, l - k);
    bTop = betaMode(k + 1, l - k);

    h = fmin(minstep, fmin(6 * SD, bTop) / steps);
    h *= stepscale;
    a = fmax(0, bTop - 6 * SD);
    SD = fmax(h, SD);

    return riemann(&fBeta, a, h, tol, SD); // x >= 1 -> fBeta(x) = 0.
}

// riemannGamma integrates fGamma from 0 to inf by Rieman sum integral.
// [[Rcpp::export]]
double riemannGamma(double k, NumericVector p, double tol = 1e-8, double stepscale = 1) {
    const double cSD = 1.25;
    double h, SD;
    init(k, p);

    SD = gammaSD(k);
    h = cSD * SD;
    if (k < 6) h *= k / 6;
    h *= stepscale;
    SD = fmax(h, SD);

    return riemann(&fGamma, 0, h, tol, SD);
}

// simpsonGamma integrates fGamma from 0 to inf by fixed steps Simpson's 1/3 rule.
// Tolerance tol is adjusted for Gamma density standard deviation.
// [[Rcpp::export]]
double simpsonGamma(double k, NumericVector p, double tol = 1e-8, double stepscale = 1) {
    const double cSD = 1.4, hmul = 2.0;
    double h, SD, hlim = 0.15;
    init(k, p);

    if (stepscale != 1) hlim = 0; // h, hlim and hmul are dependant

    SD = gammaSD(k);
    h = cSD * SD;
    if (k < 6) h *= k / 6;
    h *= stepscale;
    SD = fmax(h, SD);

    return simpson(&fGamma, 0, h, tol, SD, hlim, hmul);
}

// fBetaTop approximates the location of highest point of fBeta.
// [[Rcpp::export]]
double fBetaTop() {
    double left, right, weight;

    right = betaMode(K + 1, L - K);
    left = exp(LW / K) * 1.5;
    if (left > right)
        return right;

    //This manually fitted model from hat works quite well.
    //Solving maximun from equations goes complicated.
    weight = 0.2 + 3 * left / right + 2 * right;
    return (left + weight * right) / (1 + weight);
}

// simpsonAdaBeta integrates fBeta from 0 to 1.
// [[Rcpp::export]]
double simpsonAdaBeta(double k, NumericVector p, double abstol = 1e-7,
                      double reltol = 1e-3, int depth = 25) {
    double top, right, fa, fm, fb, I, l;

    l = init(k, p);
    if (l > 1e3) abstol *= 1e6 / (l * l);
    if (abstol < 1e-14) abstol = 1e-14;
    abstol /= 2;

    top = fBetaTop();
    fa = 0;
    fm = fBeta(top / 2);
    fb = fBeta(top);
    I = adaSimpson(&fBeta, 0, top, fa, fm, fb, 2, abstol, reltol, depth);

    right = top + 6 * betaSD(k + 1, l - k);
    right = fmin(1, right);
    fa = fb;
    fm = fBeta((top + right) / 2);
    fb = fBeta(right);
    I += adaSimpson(&fBeta, top, right, fa, fm, fb, 2, abstol, reltol, depth);
    return I;
}

// simpsonAdaGammaQuantile integrates fGammaQuantile from 0 to 1.
// [[Rcpp::export]]
double simpsonAdaGammaQuantile(double k, NumericVector p, double abstol = 1e-7,
                               double reltol = 1e-3, int depth = 25) {
    double fa, fm, fb;
    init(k, p);

    fa = fGammaQuantile(0);
    fm = fGammaQuantile(0.5);
    fb = 1;
    return adaSimpson(&fGammaQuantile, 0, 1, fa, fm, fb, 2, abstol, reltol, depth);
}

// simpsonAdaBetaQuantile integrates fBetaQuantile from 0 to 1.
// [[Rcpp::export]]
double simpsonAdaBetaQuantile(double k, NumericVector p, double abstol = 1e-7,
                              double reltol = 1e-3, int depth = 25) {
    double fa, fm, fb;
    init(k, p);

    fa = 1;
    fm = fBetaQuantile(0.5);
    fb = fBetaQuantile(1);
    return adaSimpson(&fBetaQuantile, 0, 1, fa, fm, fb, 2, abstol, reltol, depth);
}

// pTFisherCpp implements R function p.tfisher. 10-50 x faster.
// [[Rcpp::export]]
double pTFisherCpp(double lw, double l, double tau1, double tau2, double tol) {
    double lqTau, ldeltaB, ldbinom, deltaP, sum, gammaS;
    double prod = 0, cumP = 0;

    lqTau = log(tau1 / tau2);
    ldeltaB = log(tau1) - log(1 - tau1);
    ldbinom = l * log(1 - tau1);
    lw /= 2;

    if (tau1 == tau2) prod = exp(-lw); // soft TFisher
    sum = prod;

    for (double k = 1.0; k <= l; k++) {

        ldbinom += ldeltaB + log((l + 1 - k) / k);
        // ldbinom = log(R::dbinom(k, l, tau1, 0));

        if (prod > 0) {
            gammaS = sum;
            prod *= lw / k;
            sum += prod;
        } else {
            double wk = lw + k * lqTau; // lqTau=0 for soft TFisher
            if (wk < 0)
                break; // TPM main exit
            gammaS = gammaSurv(wk, k);
        }
        deltaP = (1 - gammaS) * exp(ldbinom);
        cumP += deltaP;

        if ((k > tau1 * l) && (deltaP < k * tol * (1 - cumP)))
            break;
    }
    cumP += pow(1 - tau1, l);
    return fmax(0, 1 - cumP);
}
