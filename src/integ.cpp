#include <Rcpp.h>

// riemann integrates f from a to inf by Riemann sum.
// Integration stops when relative integral increment is < tol.
double riemann(double (*f)(double), double a, double h, double tol) {
    double fa, fsum = 0, maxdist = h * 10000;

    a += h / 2;
    while (a < maxdist) {
        fa = f(a);
        fsum += fa;
        a += h;
        if (fa < fsum * tol) return fmin(1, h * fsum);
    }
    return h * fsum;
}

// adaSimpson integrates f from a to b by adaptive Simpson's 1/3 rule.
// Non standard exit condition uses relative error.
double adaSimpson(double (*f)(double), double a, double b,
                  double fa, double fm, double fb, double Iprev,
                  double reltol, int depth) {
    double h, fam, fmb, m, Ia, Ib, Iab, error;

    h = (b - a) / 4;
    fam = f(a + h);
    fmb = f(b - h);
    Ia = (fa + 4 * fam + fm) * (h / 3);
    Ib = (fm + 4 * fmb + fb) * (h / 3);

    Iab = Ia + Ib;
    error = (Iab - Iprev) / 15;
    if (fabs(error) <= Iab * reltol || depth <= 0)
        return Iab + error;

    if (error != error) return error; //NAN
    m = (a + b) / 2;
    Ia = adaSimpson(f, a, m, fa, fam, fm, Ia, reltol, depth - 1);
    Ib = adaSimpson(f, m, b, fm, fmb, fb, Ib, reltol, depth - 1);
    return Ia + Ib;
}
