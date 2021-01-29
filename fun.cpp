#include <Rcpp.h>
using namespace Rcpp;

// Finding maximum of dbeta(u, K + 1, L - K)
// f(u) = u ^ K x (1 - u) ^ (L - K - 1) (x constant)
// g(u) = log(f(u)) = K log(u) + (L - K - 1) log(1 - u)
// g'(u) = K/u - (L-K-1)/(1-u)
// =((1 - u)K - u(L - K - 1)) / (u(1 - u))
// g'(u) = 0 ->
// (1 - u)K - u(L - K - 1) = 0
// K - uK - uL + uK + u = 0
// u(-L + 1) = -K
// u = K / (L - 1) = mode of Beta(K + 1, L - K)

// Constants in integration
static double LF;
static double LC;
static double LW;
static double L;
static double K;

// [[Rcpp::export]]
void initCpp(double lw, double k, double l)
{
    LW = lw;
    K = k;
    L = l;
    LC = lgamma(l + 1) - lgamma(k + 1) - lgamma(l - k);
    LF = lgamma(k);

    // exp(LC) = L! / (K! * L-K-1!) = (L choose K+1)*(K+1)
    // See Dudbridge and Koeleman (2003).
}

inline static double pgamma(double g, double k)
{
    return R::pgamma(g, k, 1, 1, 0);
}

inline static double dgammaR(double g, double k)
{
    return R::dgamma(g, k, 1, 0);
}

inline static double qbeta(double p, double k, double l)
{
    return R::qbeta(p, k, l, 1, 0);
}

inline static double pbeta(double u, double k, double l)
{
    return R::pbeta(u, k, l, 1, 0);
}

// dgamma is density function of g ~ Gamma(k, 1).
// lf is precalculated logarithm of (k-1)!.
inline static double dgamma(double g, double lf, double k)
{
    return (g <= 0) ? 0 : exp((k - 1) * log(g) - g - lf);
}

// dbeta is density function of u ~ Beta(lc, k, l).
// lc is precalculated logarithm of l! / (k! x l-k-1!).
inline static double dbeta(double u, double lc, double k, double l)
{
    return (u <= 0 || u >= 1) ? 0 : exp(lc + (k - 1) * log(u) + (l - 1) * log(1 - u));
}

// gammaSurv = 1 - pgamma
// [[Rcpp::export]]
double gammaSurv(double g, double k)
{
    double p, q, j;
    if (g <= 0)
        return 1;

    if (g > 700 || k > 100)
        return 1.0 - pgamma(g, k);

    q = exp(-g); // q > 0
    p = q;
    for (j = 1; j < k; j++)
    {
        q *= g / j;
        p += q;
    }
    return p;
}

// betaSD returns standard deviation of Beta(a, b).
static double betaSD(double a, double b)
{
    double c = a + b;
    return sqrt(a * b / (c * c * (c + 1)));
}

// fBetaCpp is integrand over Beta density in [0, 1]. This is
// equivalent to the integrand equation in Dudbridge and Koeleman (2003).
// [[Rcpp::export]]
double fBetaCpp(double u)
{
    double g;
    if (u <= 0 || u >= 1)
        return 0;

    g = K * log(u) - LW; // exp(-g) = w/u^k
    return dbeta(u, LC, K + 1, L - K) * gammaSurv(g, K);
}

// fGamma is integrand over Gamma density in [0, inf).
static double fGamma(double g)
{
    double u, dg;
    u = exp((g + LW) / K);
    dg = dgamma(g, LF, K);
    return (u >= 1) ? dg : dg * pbeta(u, K + 1, L - K);
}

// fBetaQuantile is integrand over Beta probabilities in [0, 1].
// This is the integrand in Vsevolozhskaya et al. (2019).
static double fBetaQuantile(double p)
{
    double u, g;
    if (p <= 0)
        return 1;

    u = qbeta(p, K + 1, L - K);
    g = K * log(u) - LW;
    return gammaSurv(g, K);
}

// riemannSum integrates f from a to b by Riemann sum.
// When f -> 0, stepsize h is increased.
// Integration stops when relative increment is < tol.
static double riemannSum(double (*f)(double), double a, double b, double h,
                         double tol, double hlim, double hmul)
{
    double I, d;

    I = 0;
    a += h / 2;

    while (a < b)
    {
        d = h * f(a);
        I += d;
        a += h;

        if (d < I * tol)
            return I;
        if (d < I * hlim)
            h *= hmul;
    }
    return I;
}

// simpson integrates f from a to inf by Simpson's 1/3 rule.
// When f -> 0, stepsize h is increased.
// Integration stops when relative increment is < tol.
static double simpson(double (*f)(double), double a, double h,
                      double tol, double hlim, double hmul)
{
    double fa, fm, fb, I, d;

    fb = f(a);
    I = 0;
    while (true)
    {
        fa = fb;
        fm = f(a + h / 2);
        fb = f(a + h);
        d = (fa + 4 * fm + fb) * (h / 6);
        I += d;
        a += h;

        if (d < I * tol)
            return I;
        if (d < I * hlim)
            h *= hmul;
    }
}

// riemannBetaCpp integrates fBetaCpp from 0 to 1 by Rieman sum integral.
// [[Rcpp::export]]
double riemannBetaCpp(double lw, double k, double l,
                      double tol, double stepscale)
{
    double width, betaTop, h, a;

    initCpp(lw, k, l);

    width = 6 * betaSD(k + 1, l - k);
    betaTop = k / (l - 1);
    h = fmin(0.05, fmin(width, betaTop) / 8);
    h *= stepscale;
    a = fmax(0, betaTop - 1.5 * width);

    return riemannSum(&fBetaCpp, a, 1, h, tol, 0, 1);
}

// riemannGammaCpp integrates fGamma from 0 to inf by Rieman sum integral.
// [[Rcpp::export]]
double riemannGammaCpp(double lw, double k, double l, double tol, double stepscale)
{
    double a, h, hlim = 0.01, hmul = 1.35;

    initCpp(lw, k, l);

    h = fmax(1, (k - 1) / sqrt(k)); // sqrt(k) is Gamma SD
    h *= stepscale;
    hlim = (stepscale < 1) ? 0 : hlim;

    a = fmax(0, (k - 1) - 5 * sqrt(k)); // Gamma mode - 5 * SD

    return riemannSum(&fGamma, a, 100 * k, h, tol, hlim, hmul);
}

// simpsonGammaCpp integrates fGamma from 0 to inf by fixed steps Simpson's 1/3 rule.
// [[Rcpp::export]]
double simpsonGammaCpp(double lw, double k, double l, double tol, double stepscale)
{
    double a, h, hlim = 0.15, hmul = 2;

    initCpp(lw, k, l);

    h = fmax(1, 1.5 * (k - 1) / sqrt(k)); // sqrt(k) is Gamma SD
    h *= stepscale;
    hlim = (stepscale < 1) ? 0 : hlim; // h, hlim and hmul are dependant

    a = fmax(0, (k - 1) - 5 * sqrt(k)); // Gamma mode - 5 * SD

    return simpson(&fGamma, a, h, tol, hlim, hmul);
}

// Simpson's 1/3 rule adaptive integrator.
static double adaSimpson(double (*f)(double), double a, double b,
                         double fa, double fm, double fb, double Iprev,
                         double abstol, double reltol, int depth)
{
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

// fBetaTop approximates the location of highest point of fBetaCpp.
// [[Rcpp::export]]
double fBetaTop(double lw, double k, double l)
{
    double left, right, weight;

    right = k / (l - 1);
    left = exp(lw / k) * 1.5;
    if (left > right)
    {
        return right;
    }
    //This works, but better formulas surely exist.
    weight = 0.2 + 3 * left / right + 2 * right;
    return (left + weight * right) / (1 + weight);
}

// simpsonAdaBetaCpp integrates fBetaCpp from 0 to 1.
// [[Rcpp::export]]
double simpsonAdaBetaCpp(double lw, double k, double l,
                         double tol, double rtol, int depth)
{
    double top, right, fa, fm, fb, I;

    initCpp(lw, k, l);

    if (rtol < 1e-10)
        rtol = 1e-10;
    if (depth > 25)
        depth = 25;
    if (l > 1e3)
        tol *= 1e6 / (l * l);
    if (tol < 1e-14)
        tol = 1e-14;

    top = fBetaTop(lw, k, l);
    fa = 0;
    fm = fBetaCpp(top / 2);
    fb = fBetaCpp(top);
    I = adaSimpson(&fBetaCpp, 0, top, fa, fm, fb, 2, tol / 2, rtol, depth - 1);

    right = top + 6 * betaSD(k + 1, l - k);
    right = fmin(1, right);
    fa = fb;
    fm = fBetaCpp((top + right) / 2);
    fb = fBetaCpp(right);
    I += adaSimpson(&fBetaCpp, top, right, fa, fm, fb, 2, tol / 2, rtol, depth - 1);
    return I;
}

// simpsonAdaBetaQuantileCpp integrates fBetaQuantile from 0 to 1.
// [[Rcpp::export]]
double simpsonAdaBetaQuantileCpp(double lw, double k, double l,
                                 double tol, double rtol, int depth)
{
    double fa, fm, fb;

    initCpp(lw, k, l);

    if (tol < 1e-14)
        tol = 1e-14;
    if (rtol < 1e-10)
        rtol = 1e-10;
    if (depth > 25)
        depth = 25;

    fa = 1;
    fm = fBetaQuantile(0.5);
    fb = 0;
    return adaSimpson(&fBetaQuantile, 0, 1, fa, fm, fb, 2, tol, rtol, depth);
}

// pTFisherCpp implements R function p.tfisher. 10-50 x faster.
// [[Rcpp::export]]
double pTFisherCpp(double lw, double l, double tau1, double tau2, double tol)
{
    double qTauL, ldeltaB, lbinom, deltaP, prod, sum, gammaS, cdf;

    qTauL = log(tau1 / tau2);
    ldeltaB = log(tau1) - log(1 - tau1);
    lbinom = l * log(1 - tau1);

    lw /= 2;
    cdf = 0;
    prod = 0;
    if (tau1 == tau2) // soft TFisher
        prod = exp(-lw);
    sum = prod;

    for (double k = 1.0; k <= l; k++)
    {
        lbinom += ldeltaB + log((l + 1 - k) / k);

        if (prod > 0)
        {
            gammaS = sum;
            prod *= lw / k;
            sum += prod;
        }
        else
        {
            double wk = lw + k * qTauL;
            if (wk < 0)
                break; // TPM main exit

            gammaS = gammaSurv(wk, k);
        }
        deltaP = (1 - gammaS) * exp(lbinom);
        cdf += deltaP;

        if ((k > tau1 * l) && (deltaP < k * tol * (1 - cdf)))
            break;
    }
    cdf += pow(1 - tau1, l);
    return fmax(0, 1 - cdf);
}
