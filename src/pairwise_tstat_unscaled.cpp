#include <Rcpp.h>
using namespace Rcpp;

//' pairwise t-statistics (unscaled)
//'
//' @param ymean a numeric matrix containing group means.
//' @param stdevUnscaled a numeric matrix containing unscaled standard
//'     deviations of the group means.
// [[Rcpp::export]]
NumericMatrix pairwise_tstat_unscaled(NumericMatrix ymean, NumericMatrix stdevUnscaled)
{
    int K = ymean.ncol();
    int ngene = ymean.nrow();

    if (K != stdevUnscaled.ncol() || ngene != stdevUnscaled.nrow()) {
        stop("Unmatched size for two input matrices!");
    }

    NumericMatrix pairtUnscaled(ngene, K-1);
    for( int g = 0; g < ngene; g++ )
    {
        NumericMatrix::Row yrow = ymean( g , _ );
        NumericMatrix::Row srow = stdevUnscaled( g , _ );

        IntegerVector index(K);
        for ( int i = 0; i < K; i++ )
        {
            index[i] = i;
        }
        std::sort( index.begin(), index.end(),
              [&](int i,int j){
                  return (yrow[i] > yrow[j]);
              }
        );

        for ( int i = 1; i < K; i++ )
        {
            pairtUnscaled( g, i-1 ) = ( yrow[ index[0] ] - yrow[ index[i] ] )
            / sqrt( srow[ index[0] ] * srow[ index[0] ]
                        + srow[ index[i] ] * srow[ index[i] ] );
        }
    }

    return pairtUnscaled;

}
