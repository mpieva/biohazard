#include <assert.h>

void c_loop( double *restrict w_, const double *restrict v_, const double *restrict p7_, int n7, const double *restrict p5_, int n5_ )
{
    int n5 = n5_ * 32 ;
    double        *restrict w = __builtin_assume_aligned(  w_, 16 ) ;
    const double *restrict  v = __builtin_assume_aligned(  v_, 16 ) ;
    const double *restrict p5 = __builtin_assume_aligned( p5_, 16 ) ;
    const double *restrict p7 = __builtin_assume_aligned( p7_, 16 ) ;

    double acc = 0 ;
    for( int i = 0, k = 0 ; i != n7 ; ++i, k += n5 ) {
        double p7i = p7[i] ;
        for( int j = 0 ; j != n5 ; ++j ) {
            acc += v[k+j] * p7i * p5[j] ;
        }
    }

    acc = 1.0 / acc ;
    for( int i = 0, k = 0 ; i != n7 ; ++i ) {
        double p7i = p7[i] ;
        for( int j = 0 ; j != n5 ; ++j, ++k ) {
            w[k] += acc * v[k] * p7i * p5[j] ;
        }
    }
}

