
double riemann(double (*f)(double), double a, double h, double tol);
double simpson(double (*f)(double), double a, double h, double tol);
double adaSimpson(double (*f)(double), double a, double b, double fa,
                  double fm, double fb, double Iprev, double abstol,
                  double reltol, int depth);
void quickUniSelect(long k, NumericVector p);
double fisher(NumericVector p);
double sidak(NumericVector p);
double betaMean(double a, double b);
double betaSD(double a, double b);
inline double betaSkewness(double a, double b);
