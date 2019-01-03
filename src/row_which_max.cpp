#include <Rcpp.h>
using namespace Rcpp;

//' Multiply a number by two
//'
//' @param y A single integer.
//' @param group A single integer.
//' @export
// [[Rcpp::export]]
IntegerVector row_which_max(NumericMatrix Y)
{
    int K = Y.ncol();
    int N = Y.nrow();
    IntegerVector maxIdx(N);

    for( int j = 0; j < N; j++ )
    {
        int index = 0;
        NumericMatrix::Row yrow = Y( j , _ );
        for ( int i = 1; i < K; i++ )
        {
            if ( yrow[i] > yrow[index] )
            {
                index = i;
            }
            maxIdx( j ) = index + 1;
        }
    }

    return maxIdx;

}
