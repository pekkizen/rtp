#include <Rcpp.h>

// tfisher implements R function p.tfisher in Zhang et al.
// [[Rcpp::export]]
double tfisher(double lw, double L, double tau1, double tau2, double tol = 1e-16) {
    double lqTau, p, probBin, deltaP, survG;
    double prod = 0, cumP = 0;

    lqTau = log(tau1 / tau2); // lqTau=0 for soft TFisher
    lw /= 2;                  // Chisq stat to Gamma stat
    p = tau1;

    probBin = pow(1 - p, L); // dbinom(0, L, p)

    if (lw > 0) cumP = probBin;        // soft TFisher & tpm always
    if (tau1 == tau2) prod = exp(-lw); // soft TFisher
    survG = prod;

    for (double k = 1; k <= L; k++) {

        if (probBin > 0)
            probBin *= p / (1 - p) * (L + 1 - k) / k;
        else
            probBin = R::dbinom(k, L, p, 0);

        if (prod > 0) {
            deltaP = 1 - survG;
            prod *= lw / k;
            survG += prod;
        } else
            deltaP = R::pgamma(lw + k * lqTau, k, 1, 1, 0);

        deltaP *= probBin;
        cumP += deltaP;

        if (cumP >= 1) return 0;
        if ((k > p * L) && (deltaP < tol * (1 - cumP))) break;
    }
    return 1 - cumP;
}
