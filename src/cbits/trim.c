#include <stdint.h>

static const uint8_t compls[] =
    { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 } ;

int prim_match_reads( int i1
                    , int i2
                    , int r
                    , const uint8_t *rd1
                    , const uint8_t *qs1
                    , const uint8_t *rd2
                    , const uint8_t *qs2 )
{
    int acc = 0 ;
    while( r != 0 )
    {
        --i2 ;
        uint8_t n1 = rd1[ i1 ] ;
        uint8_t n2 = rd2[ i2 ] ;
        uint8_t q1 = qs1[ i1 ] ;
        uint8_t q2 = qs2[ i2 ] ;

        acc += (n1 & 0xF) == compls[ n2 & 0xF ] ? 0 : 5 + (q1 < q2 ? q1 : q2) ;

        ++i1 ;
        --r ;
    }
    return acc ;
}

int prim_match_ad( int off
                 , int i
                 , const uint8_t *rd
                 , const uint8_t *qs
                 , const uint8_t *ad )
{
    int acc = 0 ;
    while( i > 0 )
    {
        --i;
        acc += rd[ i+off ] == ad[ i ] ? 0 : 5 +
               (qs[ i+off ] < 25 ? qs[ i+off ] : 25) ;
    }
    return acc ;
}

