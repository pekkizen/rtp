#include <Rcpp.h>
using namespace Rcpp;

// Constants in integration
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

    // exp(LC) = L! / (K! * L-K-1!) = (L choose K+1)*(K+1)
    // See Dudbridge and Koeleman (2003).
}

inline static double pgamma(double g, double k)
{
    return R::pgamma(g, k, 1, 1, 0);
}

inline static double qbeta(double u, double k, double l)
{
    return R::qbeta(u, k, l, 1, 0);
}

// betaSD returns standard deviation of Beta(a, b).
// [[Rcpp::export]]
inline double betaSD(double a, double b)
{
    double c = a + b;
    return sqrt(a * b / (c * c * (c + 1)));
}

// dBeta is density function of u ~ Beta(lc, k, l), where lc is
// logarithm of precalculated factorial coefficient.
inline static double dbeta(double u, double lc, double k, double l)
{
    if (u <= 0 || u >= 1)
        return 0;
    return exp(lc + (k - 1) * log(u) + (l - 1) * log(1 - u));
}

// gammaSurv = 1 - pgamma
// [[Rcpp::export]]
inline static double gammaSurv(double g, double k)
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

// fDensityCpp is integrand by Beta density x Gamma survival. This is
// equivalent to the integrand equation in Dudbridge and Koeleman (2003).
// [[Rcpp::export]]
double fDensityCpp(double u)
{
    double g;
    if (u <= 0 || u >= 1)
        return 0;

    g = K * log(u) - LW; // exp(-g) = w/u^k
    return dbeta(u, LC, K + 1, L - K) * gammaSurv(g, K);
}

// fQuantileCpp is integrand by Beta quantile function qbeta.
// This is the integrand in Vsevolozhskaya et al. (2019)
double fQuantileCpp(double p)
{
    double b, g;
    if (p <= 0)
        return 1;

    b = qbeta(p, K + 1, L - K);
    g = K * log(b) - LW;
    return gammaSurv(g, K);
}

// riemannCpp integrates fDensityCpp from 0 to 1 by Rieman sum integral.
// [[Rcpp::export]]
double riemannCpp(double lw, double k, double l, double tol, double steps)
{
    double I, a, d, width, betaTop, h, start;

    initCpp(lw, k, l);

    width = 6 * betaSD(k + 1, l - k);
    betaTop = k / (l - 1);
    h = fmin(0.05, fmin(width, betaTop) / steps);

    start = fmax(0, betaTop - 1.5 * width);
    I = 0;
    a = start + h / 2;

    while (a < 1)
    {
        d = h * fDensityCpp(a);
        I += d;
        if (d < I * tol)
            return I;
        a += h;
    }
    return I;
}

// Simpson's 1/3 rule adaptive integrator.
double adaSimp(double (*f)(double), double a, double b,
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
    Ia = adaSimp(f, a, m, fa, fam, fm, Ia, abstol / 2, reltol, depth - 1);
    Ib = adaSimp(f, m, b, fm, fmb, fb, Ib, abstol / 2, reltol, depth - 1);
    return Ia + Ib;
}

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

// fDenTop approximates the location of highest point of fDensity.
// [[Rcpp::export]]
double fDenTop(double lw, double k, double l)
{
    double left, right, weight;

    right = k / (l - 1);
    left = exp(lw / k) * 1.5;
    if (left > right)
    {
        return right;
    }
    //This weighted mean works, but more accurate formulas surely exist.
    weight = 0.4 + right * 10;
    return (left + weight * right) / (weight + 1);
}

// adaSimpDensityCpp integrates fDensityCpp from 0 to 1.
// [[Rcpp::export]]
double adaSimpDensityCpp(double lw, double k, double l, double tol, double rtol, int depth)
{
    double top, end, fa, fm, fb, I;

    initCpp(lw, k, l);

    if (rtol < 1e-10)
        rtol = 1e-10;
    if (depth > 25)
        depth = 25;
    if (l > 1e3)
        tol *= 1e6 / (l * l);
    if (tol < 1e-14)
        tol = 1e-14;

    top = fDenTop(lw, k, l);
    fa = 0;
    fm = fDensityCpp(top / 2);
    fb = fDensityCpp(top);
    I = adaSimp(&fDensityCpp, 0, top, fa, fm, fb, 2, tol / 2, rtol, depth - 1);

    end = fmin(1, top + 6 * betaSD(k + 1, l - k));
    fa = fb;
    fm = fDensityCpp((top + end) / 2);
    fb = fDensityCpp(end);
    I += adaSimp(&fDensityCpp, top, end, fa, fm, fb, 2, tol / 2, rtol, depth - 1);
    return I;
}

// adaSimpQuantileCpp integrates fQuantileCpp from 0 to 1.
// [[Rcpp::export]]
double adaSimpQuantileCpp(double lw, double k, double l, double tol, double rtol, int depth)
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
    fm = fQuantileCpp(0.5);
    fb = 0;
    return adaSimp(&fQuantileCpp, 0, 1, fa, fm, fb, 2, tol, rtol, depth);
}

// pTFisherCpp implements R function p.tfisher for independant p-values.
// [[Rcpp::export]]
double pTFisherCpp(double w, double n, double tau1, double tau2, double tol)
{
    double qTauL, ldeltaB, lbinom, deltaP, prod, sum, gammaCDF, cdf;

    qTauL = log(tau1 / tau2);
    ldeltaB = log(tau1) - log(1 - tau1);
    lbinom = n * log(1 - tau1);

    w /= 2;
    cdf = 0;
    prod = 0;
    if (tau1 == tau2) // soft TFisher
        prod = exp(-w);
    sum = prod;

    for (double k = 1.0; k <= n; k++)
    {
        lbinom += ldeltaB + log((n + 1 - k) / k);

        if (prod > 0)
        {
            gammaCDF = 1 - sum;
            prod *= w / k;
            sum += prod;
        }
        else
        {
            double wk = w + k * qTauL;
            if (wk < 0)
                break; // TPM main exit

            gammaCDF = 1 - gammaSurv(wk, k);
        }
        deltaP = gammaCDF * exp(lbinom);
        cdf += deltaP;

        if ((k > tau1 * n) && (deltaP < k * tol * (1 - cdf)))
            break;
    }
    cdf += pow(1 - tau1, n);
    return fmax(0, 1 - cdf);
}
