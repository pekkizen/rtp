#include <Rcpp.h>
using namespace Rcpp;

double betaSD(double a, double b);
double betaMean(double a, double b);

// Here are functions for efficient selection of k smallest
// of n uniform(0, 1) numbers. Functions implements what is known
// as the Nth Element Algorith.

// selectSmall swaps k smallest values in range p[lo], ..., p[hi]
// to the beginning of the range: p[lo], ..., p[lo+k-1].
// For fixed (small) k and large n this is O(n) ~ n x compare.
static void selectSmall(int k, int lo, int hi, NumericVector p) {
    if (k <= 0 || k >= hi - lo + 1) return;

    int imax = lo;
    k += lo;
    for (int j = lo + 1; j < k; j++)
        if (p[imax] < p[j]) imax = j;
    double pmax = p[imax];

    for (int i = k; i <= hi; i++) {
        if (pmax > p[i]) {
            p[imax] = p[i];
            p[i] = pmax;

            imax = lo;
            for (int j = lo + 1; j < k; j++)
                if (p[imax] < p[j]) imax = j;
            pmax = p[imax];
        }
    }
}

// selectBig swaps k biggest values in range p[lo], ..., p[hi]
// to the end of the range: p[hi-k+1], ..., p[hi].
static void selectBig(int k, int lo, int hi, NumericVector p) {
    if (k <= 0 || k >= hi - lo + 1) return;

    int hilo = hi - k + 1;
    int imin = hilo;
    for (int j = hilo + 1; j <= hi; j++)
        if (p[imin] > p[j]) imin = j;
    double pmin = p[imin];

    for (int i = lo; i < hilo; i++) {
        if (pmin < p[i]) {
            p[imin] = p[i];
            p[i] = pmin;

            imin = hilo;
            for (int j = hilo + 1; j <= hi; j++)
                if (p[imin] > p[j]) imin = j;
            pmin = p[imin];
        }
    }
}

static void select(int k, int lo, int hi, NumericVector p) {
    int b = hi - lo + 1 - k;
    if (b < k)
        selectBig(b, lo, hi, p);
    else
        selectSmall(k, lo, hi, p);
}

// Hoare's quicksort partition with external pivot value.
// partition returns index j and permutation of p for which
// i <= j -> p[i] <= pivot  and i > j -> p[i] > pivot.
// If j = lo, p[lo] can be > pivot for too small pivot.
static int partition(int lo, int hi, double pivot, NumericVector p) {
    int i = lo - 1, j = hi + 1;

    while (true) {

        while (p[++i] <= pivot && i < hi)
            ;
        while (p[--j] > pivot && j > lo)
            ;
        if (i >= j) return j;

        double s = p[i];
        p[i] = p[j];
        p[j] = s;
    }
}

// quickUniSelect picks k smallest numbers in p and swaps
// them to p[0], ... , p[k-1], unordered.
// This is very efficient if numbers are near unif(0, 1) distributed.
// The k'th smallest of n unif(0, 1) numbers is distributedmBeta(k, n - k + 1).
// [[Rcpp::export]]
void quickUniSelect(int k, NumericVector p) {
    const double dist = 1;
    int n = p.size();
    if (k <= 0 || k >= n) return;

    if (n <= 25 || k <= 5 || n - k <= 5) {
        select(k, 0, n - 1, p);
        return;
    }
    double SD = betaSD(k, n + 1 - k);
    double mean = betaMean(k, n + 1 - k);
    double pivot = mean + dist * SD;

    int pi, lo = 0, hi = n - 1;

    pi = partition(0, hi, pivot, p);

    if (pi <= 0 || pi >= hi) {
        select(k, 0, hi, p);
        return;
    }
    int i = pi + 1;
    if (i >= k)
        hi = pi;
    else {
        lo = i;
        k -= i;
    }
    if (abs(i - lo - k) <= 3) { // partition missed 3 or less numbers
        select(k, lo, hi, p);
        return;
    }
    pivot *= (double)k / (pi + 1); // inter/extrapolate to k numbers
    pi = partition(lo, hi, pivot, p);

    if (p[pi] > pivot) { // pivot too small, very rare
        select(k, lo, hi, p);
        return;
    }
    i = pi - lo + 1;
    if (i == k)
        return;
    if (i > k)
        hi = pi;
    else {
        k -= i;
        lo = pi + 1;
    }
    select(k, lo, hi, p);
    return;

    // Make nth element algorithm's compatible
    // selectBig(1, 0, k - 1, p); // k'th smallest to p[k-1]
}

// For bechmarking only
// [[Rcpp::export]]
int uniSel(int k, NumericVector p) {
    quickUniSelect(k, p);
    return 1;
}
// [[Rcpp::export]]
int simpleSel(int k, NumericVector p) {
    select(k, 0, p.size() - 1, p);
    return 1;
}
// [[Rcpp::export]]
int nth_element(int k, NumericVector p) {
    std::nth_element(p.begin(), p.begin() + k, p.end());
    return 1;
}
