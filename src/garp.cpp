
#include <Rcpp.h>
using namespace Rcpp;
#include "fundec.h"

static inline double meanRTP(double K, double L) {
    return K - K * (R::digamma(K + 1) - R::digamma(L + 1));
}
static inline double varRTP(double K, double L) {
    return K + (K * K) * (R::trigamma(K + 1) - R::trigamma(L + 1));
}

static const double b[13] = {
    2.06193e-01, 1.00015e+00, -2.35081e-01, -5.24942e-05,
    8.73881e-02, -2.27766e-01, 1.64767e-02, -7.98814e-02,
    2.53813e-01, -1.85481e-02, -2.87821e-03, 9.53178e-04,
    -3.93185e-04};

static double gfun(double lw, double shape, double theta) {
    double mean = shape * theta;
    double var = shape * theta * theta;
    double sqsh = sqrt(shape);
    double d = (lw - mean) / sqrt(var);
    if (d > 10) return lw * 0.75; // TODO gfunInv

    d = b[0] + b[1] * lw + b[2] * theta + b[3] * shape +
        (b[4] + b[7] * theta + b[10] * sqsh) * d +
        (b[5] + b[8] * theta + b[11] * sqsh) * (d * d) +
        (b[6] + b[9] * theta + b[12] * sqsh) * (d * d) * d;
    return fmax(0, d);
}

static double gder(double lw, double shape, double theta) {
    double mean = shape * theta;
    double var = shape * theta * theta;
    double sqsh = sqrt(shape);
    double sd = sqrt(var);
    double d = (lw - mean) / sd;

    d = b[1] +
        (b[4] + b[7] * theta + b[10] * sqsh) / sd +
        (b[5] + b[8] * theta + b[11] * sqsh) / sd * 2 * d +
        (b[6] + b[9] * theta + b[12] * sqsh) / sd * 3 * d * d;
    return fmax(3e-308, d);
}

static double gfunInv(double y, double shape, double theta, double tol = 1e-5) {
    if (y <= 0) return 0;
    double x = y;
    for (int i = 0; i < 10; i++) {
        double delta = (y - gfun(x, shape, theta)) / gder(x, shape, theta);
        x += delta;
        if (fabs(delta) < tol) return fmax(0, x);
    }
    return R_NaN;
}

double qgarp(double p, double shape, double theta) {
    if (p <= 0) return 0;
    if (p >= 1) return R_PosInf;

    double q = R::qgamma(p, shape, theta, 1, 0);
    return gfunInv(q, shape, theta);
}

// [[Rcpp::export(name = "p.rtp.garp.lw")]]
double rtpGarpLW(double lw, double K, double L) {
    double mean = meanRTP(K, L);
    double var = varRTP(K, L);
    double shape = mean * mean / var;
    double theta = var / mean;
    lw = gfun(lw, shape, theta);
    return R::pgamma(lw, shape, theta, 0, 0);
}

// [[Rcpp::export(name = "p.rtp.garp")]]
double rtpGarp(double K, NumericVector p) {
    double L = p.size();
    if (L < K || L < 1) return ERR;
    if (K == 1 && L == 1) return p[0];
    if (K == 1) return single(p);
    if (K == L) return fisher(p);
    return rtpGarpLW(statRTP(K, p), K, L);
}

// [[Rcpp::export(name = "rgarp")]]
NumericVector randGarp(long samplesize, double k, double l) {
    NumericVector v(samplesize);
    double mean = meanRTP(k, l);
    double var = varRTP(k, l);
    double shape = mean * mean / var;
    double theta = var / mean;

    for (long i = 0; i < samplesize; i++) {
        v[i] = gfunInv(R::rgamma(shape, theta), shape, theta);
    }
    return v;
}

// [[Rcpp::export]]
double garpError(double r1, double r2, int samples, double plim,
                 double K, double L) {
    double pg, pr, lw = 0, err, maxerr = 0;
    if (r1 >= r2) return ERR;
    int j = 0;
    for (int i = 0; i < samples; i++) {
        pr = 0;
        while (pr < plim) {
            lw = r1 + (r2 - r1) * R::unif_rand();
            pr = rtpDgammaLW(lw, K, L, 1e-16, 0.1);
            if (++j > 10 * samples) return ERR;
        }
        pg = rtpGarpLW(lw, K, L);

        // err = fmax(pg, pr) / fmin(pg, pr);
        err = fabs(pg - pr) / pr;
        if (err > maxerr) {
            maxerr = err;
        }
        // maxerr += -log10(err);
    }
    // return maxerr / samples;
    return maxerr;
}
