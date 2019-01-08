#include <Rcpp.h>
using namespace Rcpp;


NumericVector rank_tie_average(NumericVector x, bool decreasing)
{
    R_xlen_t sz = x.size();
    IntegerVector w = seq(0, sz - 1);
    std::sort(w.begin(), w.end(),
              [&](R_xlen_t i,R_xlen_t j){
                  return decreasing ? (x[i] > x[j]) : (x[i] < x[j]);
              }
    );

    NumericVector r = no_init_vector(sz);
    for (R_xlen_t n, i = 0; i < sz; i += n) {
        n = 1;
        while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;
        for (R_xlen_t k = 0; k < n; k++) {
            r[w[i + k]] = i + (n + 1) / 2.;
        }
    }

    return r;
}


IntegerVector order(NumericVector x, bool decreasing)
{
    R_xlen_t sz = x.size();
    IntegerVector w = seq(0, sz - 1);
    std::sort(w.begin(), w.end(),
              [&](R_xlen_t i,R_xlen_t j){
                  return decreasing ? (x[i] > x[j]) : (x[i] < x[j]);
              }
    );
    return w;
}


IntegerVector rank_tie_last(NumericVector x, bool decreasing)
{
    IntegerVector w = (!decreasing) ? order(x, true) : order(x, false);
    std::reverse(w.begin(), w.end());
    IntegerVector r = order(as<NumericVector>(w), false);
    return r;
}


NumericVector reorder_by_rank(NumericVector x, IntegerVector index)
{
    if (x.length() != index.length()) {
        stop("Unmatched size for two input vectors!");
    }

    NumericVector y = no_init_vector(x.length());
    for ( int i = 0; i < x.length(); i++ )
    {
        y[index[i]] = x[i];
    }

    return y;
}

//' p-value by weighted estimation
//'
//' @param tt a vector of test statistics.
//' @param ttperm a matrix of test statistics from permutaitons.
//' @param W a matrix containing weights for each spot in ttperm.
//' @export
// [[Rcpp::export]]
NumericVector pvalue_weighted_est(NumericVector tt, NumericVector ttperm,
                                  NumericVector w)
{
    R_xlen_t L1 = tt.length();
    R_xlen_t L2 = ttperm.length();

    if (L2 != w.length()) {
        stop("Unmatched size for permutation matrix and weight matrix!");
    }

    NumericVector ttall = no_init_vector(L1 + L2);
    NumericVector wall = no_init_vector(L1 + L2);

    std::copy(tt.begin(), tt.end(), ttall.begin());
    std::copy(ttperm.begin(), ttperm.end(), ttall.begin() + L1);
    std::copy(w.begin(), w.end(), wall.begin() + L1);
    for ( int j = 0; j < L1; j++ ) {
        wall[ j ] = 1.0;
    }

    IntegerVector r0 = rank_tie_last(ttall, true);

    NumericVector wallr = reorder_by_rank(wall, r0);
    NumericVector wrank = cumsum(wallr);
    NumericVector r2 = rank_tie_average(tt, true);

    NumericVector pv = no_init_vector(L1);
    double s = sum(w);
    for ( int j = 0; j < L1; j++ )
    {
        pv[ j ] = ( wrank[ r0[j] ] - r2[ j ] ) / s;

    }

    return pv;
}



