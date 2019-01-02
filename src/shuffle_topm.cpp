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
    IntegerMatrix tmp(ncombM, nsample, &shuffletype[0][0]);
    Rcout << "shuffletype matrix is " << std::endl <<  tmp << std::endl;

    for( int g = 0; g < ngene; g++ )
    {
        int* shuffle = &shuffletype[geneSubset(g)][0];
        NumericMatrix::Row v = y( g , _ );
        for ( int i = 0; i < nsample; i++ )
        {
            YMat( g , i ) = v[ shuffle[i] ];
        }
    }
    Rcout << "YMat matrix is " << std::endl << YMat << std::endl;

    if (weights.isNotNull())
    {
        Rcout << "weights is not null " << std::endl;
        NumericMatrix WMatrixOrig(weights);
        NumericMatrix WMat(ngene, nsample);
        for( int g = 0; g < ngene; g++ )
        {
            int* shuffle = &shuffletype[geneSubset[g]][0];
            NumericMatrix::Row w = WMatrixOrig( g , _ );
            for ( int i = 0; i < nsample; i++ )
            {
                WMat( g , i ) = w[ shuffle[i] ];
            }
            Rcout << "WMat matrix is " << std::endl << WMat << std::endl;
        }
        List L = List::create(Named("y") = YMat , Named("weights") = WMat);
        return L;
    } else {
        List L = List::create(Named("y") = YMat ,
                              Named("weights") = R_NilValue);
        return L;
    }
}
