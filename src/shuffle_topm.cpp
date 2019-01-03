#include <Rcpp.h>
#include <random>
using namespace Rcpp;

//' Multiply a number by two
//'
//' @param y A single integer.
//' @param group A single integer.
//' @export
// [[Rcpp::export]]
List shuffle_topm(NumericMatrix y, IntegerVector group,
    Nullable<NumericMatrix> weights,
    IntegerMatrix combM, IntegerVector geneSubset, int seed)
{
    int ncombM = combM.nrow();
    int ngene = y.nrow();
    int nsample = y.ncol();
    NumericMatrix YMat(ngene, nsample);

    int shuffletype[ncombM][nsample];
    int a[nsample];
    bool b[nsample];
    int c[nsample];

    for( int j = 0; j < ncombM; j++ )
    {
        IntegerMatrix::Row row = combM( j , _ );

        int n = 0;
        for ( int i = 0; i < nsample; i++ )
        {
            a[i] = i;
            b[i] = std::find(row.begin(), row.end(), group(i))
                != row.end();
            if ( b[i] )
            {
                c[n] = i;
                n++;
            }

        }
        std::shuffle(c, c+n, std::default_random_engine(seed));
        int ni = 0;
        for ( int i = 0; i < nsample; i++ )
        {
            if ( b[i] )
            {
                a[i] = c[ni];
                ni++;
            }
        }

        std::copy( a, a + nsample, &shuffletype[j][0] );
    }

    for( int g = 0; g < ngene; g++ )
    {
        int* shuffle = &shuffletype[geneSubset(g) - 1][0];
        NumericMatrix::Row v = y( g , _ );
        for ( int i = 0; i < nsample; i++ )
        {
            YMat( g , i ) = v[ shuffle[i] ];
        }
    }

    checkUserInterrupt();

    if (weights.isNotNull())
    {
        NumericMatrix WMatrixOrig(weights);
        NumericMatrix WMat(ngene, nsample);
        for( int g = 0; g < ngene; g++ )
        {
            int* shuffle = &shuffletype[geneSubset(g) - 1][0];
            NumericMatrix::Row w = WMatrixOrig( g , _ );
            for ( int i = 0; i < nsample; i++ )
            {
                WMat( g , i ) = w[ shuffle[i] ];
            }
        }
        List L = List::create(Named("y") = YMat , Named("weights") = WMat);
        return L;
    } else {
        List L = List::create(Named("y") = YMat ,
                              Named("weights") = R_NilValue);
        return L;
    }

}
