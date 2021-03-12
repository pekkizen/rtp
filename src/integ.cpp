#include <Rcpp.h>

// riemann integrates f from a to inf by Riemann sum. Integration
// stops when relative integral value on SD distance is < tol.
// Assumed unimodal f(x) -> 0 when x -> inf.
double riemann(double (*f)(double), double a,
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
double simpson(double (*f)(double), double a, double h, double tol,
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
double adaSimpson(double (*f)(double), double a, double b,
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

    // bisect only if Ia differs enough from 0 ???
    if (Ia > abstol && Ia > Iab * reltol)
        Ia = adaSimpson(f, a, m, fa, fam, fm, Ia, abstol, reltol, depth - 1);

    if (Ib > abstol && Ib > Iab * reltol)
        Ib = adaSimpson(f, m, b, fm, fmb, fb, Ib, abstol, reltol, depth - 1);
    return Ia + Ib;
}
