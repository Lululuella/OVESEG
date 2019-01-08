#include <Rcpp.h>
using namespace Rcpp;

//' min value for each row
//'
//' @param Y a numeric matrix.
// [[Rcpp::export]]
NumericVector row_min(NumericMatrix Y)
{
    int K = Y.ncol();
    int N = Y.nrow();
    NumericVector rowMin(N);

    for( int j = 0; j < N; j++ )
    {
        NumericMatrix::Row yrow = Y( j , _ );
        double minvalue = yrow[0];
        for ( int i = 1; i < K; i++ )
        {
            if ( yrow[i] < minvalue )
            {
                minvalue = yrow[i];
            }
            rowMin( j ) = minvalue;
        }
    }

    return rowMin;

}
