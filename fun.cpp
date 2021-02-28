#include <Rcpp.h>
using namespace Rcpp;

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

// Global "constants" in integration
static double K;
static double L;
static double LBeta;
static double LF;
static double LW;
static int ERR = 0;

// Benchmark baseline reference.
// [[Rcpp::export]]
static double fNull(double a) {
    return a;
}

// lbeta returns logarithm of Beta(a, b) = Gamma(a) * Gamma(b) / Gamma(a+b).
// exp(-lbeta(K+1, L-K)) = L! / (K! * L-K-1!).
inline static double lbeta(double a, double b) {
    return lgamma(a) + lgamma(b) - lgamma(a + b);
}

static void selectSmalls(int k, NumericVector p);

// [[Rcpp::export]]
double init(double k, NumericVector p) {
    L = p.size();
    K = k;
    ERR = 0;
    if (L < K) {
        ERR = -1;
        return ERR;
    }
    LW = 0;
    selectSmalls(K, p);
    for (int i = 0; i < k; i++)
        LW += log(p[i]);
    LBeta = -lbeta(K + 1, L - K);
    LF = lgamma(K);
    return L;
}

// Beta standard deviation.
// betaSD(K+1, L-K) = ~sqrt(K) / L, for small K and big L.
// [[Rcpp::export]]
static double betaSD(double a, double b) {
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

// Beta quantile/inverse CDF function.
inline static double qbeta(double p, double a, double b) {
    if (p >= 1) return 1;
    if (p <= 0) return 0;
    return R::qbeta(p, a, b, 1, 0);
}

// Beta distribution function.
inline static double pbeta(double x, double a, double b) {
    return R::pbeta(x, a, b, 1, 0); // (x >= 1 || x <= 0) returns 0
}

// Fast binomial survival function.
// sbinom(K, L, b) = pbeta(b, K+1, L-K) = 1 - pbinom(K, L, b)
static double sbinom(double k, double n, double p) {
    double prob, cdf;
    if (p <= 0) return 0;
    if (p >= 1) return 1;

    prob = pow(1 - p, n); // dbinom(0, n, p)
    cdf = prob;

    for (double j = 1.0; j <= k; j++) {
        prob *= p * (n + 1 - j) / ((1 - p) * j);
        cdf += prob;
    }
    if (cdf > (1.0 - 1e-10) || cdf == 0)
        // these give also values in [0, 1e-16]
        // return  R::pbinom(k, n, p, 0, 0)
        return pbeta(p, k + 1, n - k);

    return 1 - cdf;
}

// Beta density function.
// lg is -log((a+b-1)! / ((a-1)! * (b-1)!)) = -lbeta(a, b).
inline static double dbeta(double x, double lg, double a, double b) {
    if (x <= 0 || x >= 1) return 0;

    return exp(lg + (a - 1) * log(x) + (b - 1) * log(1 - x));
}

// Hight of Beta density function. Function value at mode.
// [[Rcpp::export]]
double dbetaHight(double a, double b) {
    return dbeta(betaMode(a, b), -lbeta(a, b), a, b);
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

// Fast gamma survival function for positive integer shape k and rate 1.
static double sgamma(double g, double k) {
    double s, p;
    if (g <= 0) return 1;

    if (g > 700 || k > 100)
        return R::pgamma(g, k, 1, 0, 0); // right tail

    p = exp(-g); // p > 0
    s = p;
    for (double j = 1.0; j < k; j++) {
        p *= g / j;
        s += p;
    }
    return fmin(1, s);
}

// Gamma density function (g,k) = e^-g * g^(k-1) / (k-1)!
// lf is precalculated logarithm of (k-1)!.
inline static double dgamma(double g, double lf, double k) {
    if (g <= 0) return 0;

    return exp((k - 1) * log(g) - g - lf);
}

// Gamma quantile/inverse CDF function.
inline static double qgamma(double p, double k) {
    if (p <= 0) return 0;
    if (p >= 1) return R_PosInf;
    return R::qgamma(p, k, 1, 1, 0);
}

// fBetaD is integrand over Beta density in [0, 1]. This is equivalent
// to the integrand equation in Dudbridge and Koeleman (2003).
// [[Rcpp::export]]
double fBetaD(double b) {
    if (b <= 0 || b >= 1) return 0;

    double g = K * log(b) - LW;
    return dbeta(b, LBeta, K + 1, L - K) * sgamma(g, K);
}

// fGammaD is integrand over Gamma density in [0, inf).
// [[Rcpp::export]]
double fGammaD(double g) {
    if (g <= 0) return 0;

    double b = exp((g + LW) / K);
    return dgamma(g, LF, K) * sbinom(K, L, b);
    // return dgamma(g, LF, K) * pbeta(b, K + 1, L - K);
}

// fBetaQ is integrand over Beta probabilities in [0, 1].
// This is the integrand in Vsevolozhskaya et al. (2019).
// [[Rcpp::export]]
double fBetaQ(double p) {
    if (p <= 0) return 1;

    double b = qbeta(p, K + 1, L - K); // Inverse beta CDF method
    double g = K * log(b) - LW;
    return sgamma(g, K);
}

// fGammaQ is integrand over Gamma probabilities in[0, 1].
// fGammaQ(1 - u) is very close to fBetaQ(u). ~2 x faster than fBetaQ.
// [[Rcpp::export]]
double fGammaQ(double p) {
    if (p >= 1) return 1;

    double g = qgamma(p, K); // Inverse gamma CDF method
    double b = exp((g + LW) / K);
    return sbinom(K, L, b);
    // return pbeta(b, K + 1, L - K);
}

// riemann integrates f from a to inf by Riemann sum. Integration
// stops when relative integral value on SD distance is < tol.
// Assumed unimodal f(x) -> 0 when x -> inf.
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
        if (fa < fsum * tol) return fmin(1, h * fsum);
    }
}

// simpson integrates f from a to inf by Simpson's 1/3 rule.
// When f gets small, stepsize h is increased by hlim and hmul.
// Integration stops when relative integral value on SD distance is < tol.
// Assumed unimodal f(x) -> 0 when x -> inf.
static double simpson(double (*f)(double), double a, double h, double tol,
                      double SD, double hlim, double hmul) {
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
        if (fab < I * tol) return fmin(1, I);

        if (fab * h < I * hlim) h *= hmul;
    }
}

// adaSimpson integrates f from a to b by adaptive Simpson's 1/3 rule.
// Non standard exit condition uses also relative tolerance, || not &&.
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
    abstol /= 2;

    // bisect only if Ia differs enough from 0.
    if (Ia > abstol && Ia > Iab * reltol)
        Ia = adaSimpson(f, a, m, fa, fam, fm, Ia, abstol, reltol, depth - 1);

    if (Ib > abstol && Ib > Iab * reltol)
        Ib = adaSimpson(f, m, b, fm, fmb, fb, Ib, abstol, reltol, depth - 1);
    return Ia + Ib;
}

// fBetaDtop approximates the location of highest point of fBetaD.
// This is manually fitted from hat model.
// Maximum from equations has no closed form solution?
// [[Rcpp::export]]
double fBetaDtop() {
    double left, right, weight;

    right = betaMode(K + 1, L - K);
    left = exp(LW / K) * 1.5;
    if (left > right) return right;

    weight = 0.2 + 3 * left / right + 2 * right;
    return (left + weight * right) / (1 + weight);
}

// riemannBeta integrates fBetaD from 0 to 1 by Riemann sum integral.
// [[Rcpp::export]]
double riemannBeta(double k, NumericVector p, double tol = 1e-10, double stepscale = 1) {
    const double steps = 8.0, maxstep = 0.05;
    double betaTop, h, SD, l;
    l = init(k, p);
    if (l < 0) return ERR;

    SD = betaSD(k + 1, l - k);
    betaTop = betaMode(k + 1, l - k);
    h = fmin(maxstep, fmin(6 * SD, betaTop) / steps);
    h *= stepscale;

    return riemann(&fBetaD, 0, h, tol, SD); // x >= 1 -> fBetaD(x) = 0.
}

// simpsonAdaBeta integrates fBetaD from 0 to 1.
// [[Rcpp::export]]
double simpsonAdaBeta(double k, NumericVector p, double abstol = 1e-7, double reltol = 1e-3) {
    double top, fa, fm, fb, I;
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

// riemannGamma integrates fGammaD from 0 to inf by Rieman sum integral.
// [[Rcpp::export]]
double riemannGamma(double k, NumericVector p, double tol = 1e-10, double stepscale = 1) {
    const double cSD = 1.25;
    double h, SD;
    if (init(k, p) < 0) return ERR;

    SD = gammaSD(k);
    h = cSD * SD;
    if (k < 6) h *= k / 6;
    h *= stepscale;

    return riemann(&fGammaD, 0, h, tol, SD);
}

// simpsonGamma integrates fGammaD from 0 to inf by fixed step Simpson's 1/3 rule.
// [[Rcpp::export]]
double simpsonGamma(double k, NumericVector p, double tol = 1e-10, double stepscale = 1) {
    const double cSD = 1.5;
    double h, SD, hlim = 0.15, hmul = 1.5;
    if (init(k, p) < 0) return ERR;

    if (stepscale != 1) hlim = 0; // h, hlim and hmul are dependant

    SD = gammaSD(k);
    h = cSD * SD;
    if (k < 6) h *= k / 6;
    h *= stepscale;

    return simpson(&fGammaD, 0, h, tol, SD, hlim, hmul);
}

// simpsonAdaGammaQ integrates fGammaQ from 0 to 1.
// [[Rcpp::export]]
double simpsonAdaGammaQ(double k, NumericVector p, double abstol = 1e-7,
                        double reltol = 1e-3, int depth = 25) {
    double fa, fm, fb;
    if (init(k, p) < 0) return ERR;

    fa = fGammaQ(0);
    fm = fGammaQ(0.5);
    fb = 1;
    return adaSimpson(&fGammaQ, 0, 1, fa, fm, fb, 2, abstol, reltol, depth);
}

// simpsonAdaBetaQ integrates fBetaQ from 0 to 1.
// [[Rcpp::export]]
double simpsonAdaBetaQ(double k, NumericVector p, double abstol = 1e-7,
                       double reltol = 1e-3, int depth = 25) {
    double fa, fm, fb;
    if (init(k, p) < 0) return ERR;

    fa = 1;
    fm = fBetaQ(0.5);
    fb = fBetaQ(1);
    return adaSimpson(&fBetaQ, 0, 1, fa, fm, fb, 2, abstol, reltol, depth);
}

// pTFisher implements R function p.tfisher from Zhang et al (2020).
// [[Rcpp::export]]
double pTFisher(double lw, double n, double tau1, double tau2, double tol = 1e-14) {
    double lqTau, p, probBin, deltaP, survG;
    double prod = 0, cumP = 0, fastBin = false;

    lqTau = log(tau1 / tau2); // lqTau=0 for soft TFisher
    p = tau1;
    lw /= 2; // Chisq stat to Gamma stat.

    probBin = pow(1 - p, n); // dbinom(0, n, p)
    cumP = probBin;
    if (probBin > 0) fastBin = true;

    if (tau1 == tau2) prod = exp(-lw); // soft TFisher
    survG = prod;

    // For testing. Use only R::functions.
    if (tol == 0) {
        fastBin = false;
        prod = 0;
        tol = 1e-14;
    }
    for (double k = 1.0; k <= n; k++) {

        if (fastBin)
            probBin *= p * (n + 1 - k) / ((1 - p) * k);
        else
            probBin = R::dbinom(k, n, p, 0);

        if (prod > 0) {
            deltaP = 1 - survG;
            prod *= lw / k;
            survG += prod;
        } else
            deltaP = R::pgamma(lw + k * lqTau, k, 1, 1, 0);

        deltaP *= probBin;
        cumP += deltaP;

        if ((k > p * n) && (deltaP < tol * ((1 + 1e-10) - cumP)))
            break;
    }
    return fmax(0, 1 - cumP);
}

// Below are functions for efficient selection of k smallest
// of n Uniform(0, 1) numbers.

// select swaps k smallest values in range p[lo], ..., p[hi]
// to p[lo], ..., p[lo+k-1]. For small k and large n this is O(n) operation.
static void select(int k, int lo, int hi, NumericVector p) {
    if (k <= 0 || k >= hi - lo) return;

    int max = 0;
    k += lo;
    for (int j = lo; j < k; j++)
        if (p[max] < p[j]) max = j;
    double pmax = p[max];

    for (int i = k; i <= hi; i++) {
        if (pmax > p[i]) {
            p[max] = p[i];
            p[i] = pmax;

            max = 0;
            for (int j = lo; j < k; j++)
                if (p[max] < p[j]) max = j;
            pmax = p[max];
        }
    }
}

// This is Hoare's quicksort partition with parameter pivot.
static int partition(int lo, int hi, double pivot, NumericVector p) {

    int i = lo - 1, j = hi + 1;

    while (true) {
        do {
            i++;
        } while (p[i] < pivot && i < hi);

        do {
            j--;
        } while (p[j] > pivot && j > lo);

        if (i >= j) return j;

        double s = p[i];
        p[i] = p[j];
        p[j] = s;
    }
}

// [[Rcpp::export]]
int selectBench(int k, NumericVector p) {
    selectSmalls(k, p);
    return k;
}

// selectSmalls picks k smallest numbers in p
// and swaps them to p[0], ... , p[k-1], unordered.
// This is very efficient if numbers are near uniform(0, 1) distributed and
// not very bad in other cases.
// Equivalent std function: std::nth_element(p.begin(), p.begin() + k, p.end());
// R function: p <- sort(p, partial = c(1:k))
static void selectSmalls(int k, NumericVector p) {
    int n = p.size();
    if (k <= 0 || k >= n) return;

    if (n <= 20 || k <= 3) {
        select(k, 0, n - 1, p);
        return;
    }

    // The k'th smallest of n unif(0, 1) numbers ~Beta(k, n-k+1).
    double K = k, N = n + 1;
    double SD = sqrt(K * (N - K) / (N * N * (N + 1)));
    double mean = K / N; // k / (n+1)
    double pivot = mean + 2 * SD;

    int pi, hi = n - 1;

    pi = partition(0, hi, pivot, p);

    if (pi <= 0 || pi >= hi) {
        select(k, 0, hi, p);
        return;
    }
    if (pi + 1 >= k) hi = pi;

    if (pi > k) {
        pivot *= K / (pi + 1);
        pi = partition(0, hi, pivot, p);
    }
    if (pi + 1 < k) {
        select(k - pi - 1, pi + 1, hi, p);
        return;
    }
    select(k, 0, pi + 1, p); // why pi+1, not pi?
}
