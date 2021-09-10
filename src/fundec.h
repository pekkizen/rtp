
double riemann(double (*f)(double), double a, double h, double tol);
double adaSimpson(double (*f)(double), double a, double b, double fa,
                  double fm, double fb, double Iprev, double abstol,
                  int depth);
void kSelect(long k, NumericVector p);
double fisher(NumericVector p);
double single(NumericVector p);
double betaMean(double a, double b);
double betaSD(double a, double b);
double betaSkewness(double a, double b);
double statRTP(long k, NumericVector q);
double rtpDgammaLW(double lw, double k, double l, double tol, double stepscale);

#define ERR -1
#define OK 2
