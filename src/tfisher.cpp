#include <Rcpp.h>

// tfisher implements R function p.tfisher in Zhang et al.
// [[Rcpp::export]]
double tfisher(double lw, double n, double tau1, double tau2, double tol = 1e-16) {
    double lqTau, p, probBin, deltaP, survG;
    double prod = 0, cumP = 0, fastBin = false;

    lqTau = log(tau1 / tau2); // lqTau=0 for soft TFisher
    p = tau1;
    lw /= 2; // Chisq stat to Gamma stat

    probBin = pow(1 - p, n);    // dbinom(0, n, p)
    if (lw > 0) cumP = probBin; // soft TFisher & tpm always
    if (probBin > 0) fastBin = true;

    if (tau1 == tau2) prod = exp(-lw); // soft TFisher
    survG = prod;

    // For testing. Only R::functions used
    if (tol == 0) {
        fastBin = false;
        prod = 0;
        tol = 1e-16;
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

        if (cumP >= 1) return 0;
        if ((k > p * n) && (deltaP < tol * (1 - cumP))) break;
    }
    return 1 - cumP;
}
