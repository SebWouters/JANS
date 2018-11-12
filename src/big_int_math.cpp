/*
   JANS: just another number sieve
   Copyright (C) 2018 Sebastian Wouters

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <assert.h>

#include "big_int.h"

int jans::big_int::__sum3set__( ubase_t * r, ubase_t * a, const int la, ubase_t * b, const int lb ){

   // r = a + b

   if (( r != a ) && ( r != b )){ __clear__( r ); }
   const int upper = ( ( la > lb ) ? la : lb );

   ucarry_t z = 0;

   for ( int i = 0; i < upper; i++ ){
      z = ( ( z + a[ i ] ) + b[ i ] );
      r[ i ] = z & __11111111__;
      z = z >> BLOCK_BIT;
   }

   if ( z == 0 ){ return upper; }
   assert( upper < NUM_BLOCK ); // If z != 0, z == 1, and we need to place z in upper: Overflow exception
   r[ upper ] = 1;
   return ( upper + 1 );

}

int jans::big_int::__sum1__( ubase_t * r, const int lr, const ubase_t b ){

   // r = r + b

   ucarry_t z = b;

   int ir = 0;
   while ( ( ir < NUM_BLOCK ) && ( z != 0 ) ){
      z = z + r[ ir ];
      r[ ir ] = z & __11111111__;
      z = z >> BLOCK_BIT;
      ir++;
   }

   assert( z == 0 ); // If z != 0, while loop stopped because ( ir == NUM_BLOCK ): Overflow exception

   return ( ( ir > lr ) ? ir : lr );

}

int jans::big_int::__diff3set__( ubase_t * r, ubase_t * a, ubase_t * b ){

   // r = a - b

   const int comp = __compare__( a, b );
   assert( comp >= 0 );

   ucarry_t add = 0;
   ucarry_t sub = 0; // serves as carry

   for ( int i = 0; i < comp; i++ ){
      add = a[ i ];
      sub = b[ i ] + sub;
      if ( add >= sub ){
         r[ i ] = ( add - sub );
         sub = 0;
      } else { // base > sub > add >= 0 --> base + add - sub > 0
         r[ i ] = ( ( add + ( ( ucarry_t )( 1UL ) << BLOCK_BIT ) ) - sub );
         sub = 1;
      }
   } // At end of loop sub == 0, as comp >= 0
   for ( int i = comp; i < NUM_BLOCK; i++ ){ r[ i ] = 0; } // Non-set part of the array is cleared here.

   for ( int i = comp; i > 0; i-- ){ if ( r[ i - 1 ] != 0 ){ return i; } }
   return 0;

}

int jans::big_int::__diff1__( ubase_t * r, const ubase_t b ){

   // r = r - b

   ucarry_t add = 0;
   ucarry_t sub = b; // serves as carry

   int ir = 0;
   while ( ( ir < NUM_BLOCK ) && ( sub != 0 ) ){
      add = r[ ir ];
      if ( add >= sub ){
         r[ ir ] = ( add - sub );
         sub = 0;
      } else { // base > sub > add >= 0 --> base + add - sub > 0
         r[ ir ] = ( ( add + ( ( ucarry_t )( 1UL ) << BLOCK_BIT ) ) - sub );
         sub = 1;
      }
      ir++;
   }

   assert( sub == 0 ); // If sub != 0, while loop stopped because ( ir == NUM_BLOCK ): Overflow exception

   for ( int i = NUM_BLOCK; i > 0; i-- ){ if ( r[ i - 1 ] != 0 ){ return i; } }
   return 0;

}

int jans::big_int::__mult3set__( ubase_t * r, ubase_t * a, const int la, ubase_t * b, const int lb ){

   // r = a * b

   __clear__( r );
   const int lr = __mult3add__( r, 0, a, la, b, lb );
   return lr;

}

int jans::big_int::__mult3add__( ubase_t * r, const int lr, ubase_t * a, const int la, ubase_t * b, const int lb ){

   // r += a * b

   assert( ( la + lb - 2 ) < NUM_BLOCK ); // ia + ib <= la + lb - 2 < NUM_BLOCK: Overflow exception

   int lrnew = lr;
   for ( int ib = 0; ib < lb; ib++ ){
      lrnew = jans::big_int::__mult2add__( r, lrnew, a, la, b[ ib ], ib );
   }
   return lrnew;

}

int jans::big_int::__mult2set__( ubase_t * r, ubase_t * a, const int la, const ubase_t b, const int shift ){

   // r[ shift : ] = b * a[ : ]

   __clear__( r );
   const int upper = shift + la;
   assert( ( upper - 1 ) < NUM_BLOCK ); // shift + ia <= shift + la - 1 < NUM_BLOCK: Overflow exception

         ucarry_t z = 0;
   const ucarry_t f = b;

   for ( int ia = 0; ia < la; ia++ ){
      z = z + ( f * a[ ia ] );
      r[ shift + ia ] = z & __11111111__;
      z = z >> BLOCK_BIT;
   }

   if ( z == 0 ){ return upper; }
   assert( upper < NUM_BLOCK ); // If z != 0, we need to place z in upper: Overflow exception
   r[ upper ] = z;
   return ( upper + 1 );

}

int jans::big_int::__mult2add__( ubase_t * r, const int lr, ubase_t * a, const int la, const ubase_t b, const int shift ){

   // r[ shift : ] += b * a[ : ]; This would be lapack "axpy" with a shift

   const int upper = shift + la;
   assert( ( upper - 1 ) < NUM_BLOCK ); // shift + ia <= shift + la - 1 < NUM_BLOCK: Overflow exception

   ucarry_t w = 0;
   ucarry_t x = 0;
   ucarry_t z = 0;

   for ( int ia = 0; ia < la; ia++ ){
      w = r[ shift + ia ];
      x = a[ ia ];
      z = z + w + ( x * b );
      r[ shift + ia ] = z & __11111111__;
      z = z >> BLOCK_BIT;
   }

   if ( z == 0 ){ return ( ( upper > lr ) ? upper : lr ); }

   int ir = upper;
   while ( ( ir < NUM_BLOCK ) && ( z != 0 ) ){
      w = r[ ir ];
      z = z + w;
      r[ ir ] = z & __11111111__;
      z = z >> BLOCK_BIT;
      ir++;
   }

   assert( z == 0 ); // If z != 0, while loop stopped because ( ir == NUM_BLOCK ): Overflow exception

   return ( ( ir > lr ) ? ir : lr );

}

int jans::big_int::__scal1__( ubase_t * r, const int lr, const ubase_t b ){

   // r = r * b

         ucarry_t z = 0;
   const ucarry_t f = b;

   for ( int ir = 0; ir < lr; ir++ ){
      z = z + ( f * r[ ir ] );
      r[ ir ] = z & __11111111__;
      z = z >> BLOCK_BIT;
   }

   if ( z == 0 ){ return lr; }
   assert( lr < NUM_BLOCK ); // If z != 0, we need to place z in lr: Overflow exception
   r[ lr ] = z & __11111111__;
   return ( lr + 1 );

}

ubase_t jans::big_int::__divide__( ubase_t * q, int & lq, const ubase_t div ){

   // Solves for n = q * d + r, with r < d; whereby initially (q, lq) contains (n, ln).

   if ( lq == 0 ){ return 0; }

   ucarry_t rem = 0;
   const int start = lq - 1;
   lq = 0;

   for ( int iq = start; iq >= 0; iq-- ){

      rem = ( rem << BLOCK_BIT ) + q[ iq ];
      q[ iq ] = rem / div;
      if ( q[ iq ] > 0 ){
         rem = rem - ( ( ( ucarry_t ) div ) * q[ iq ] );
         if ( lq == 0 ){ lq = iq + 1; }
      }
   }

   return ( ( ubase_t )( rem ) );

}

void jans::big_int::__divide__( ubase_t * q, int & lq, ubase_t * r, int & lr, ubase_t * d, const int ld ){

/*  Outer case 1: r = r_i b^n + ( b - 1 ) b^( n - 1 ) + ( b - 1 ) b^( n - 2 ) + ...
                  d = d_i b^n +     ( 0 ) b^( n - 1 ) +     ( 0 ) b^( n - 2 ) + ...
         In this case, q_guess = lower( r_i / d_i ) is correct.

    Outer case 2: r = r_i b^n +     ( 0 ) b^( n - 1 ) +     ( 0 ) b^( n - 2 ) + ...
                  d = d_i b^n + ( b - 1 ) b^( n - 1 ) + ( b - 1 ) b^( n - 2 ) + ...
         In this case, q_guess = lower( r_i / d_i ) may be an overestimation.

         In this (extremal outer) case, q_solve = lower( r_i b^n / ( ( d_i + 1 ) b^n - 1 ) ) >= lower( r_i / ( d_i + 1 ) )
*/

   // Solves for n = q * d + r, with r < d; whereby initially (r, lr) contains (n, ln).

   __clear__( q );
   lq = 0;

   const int comp = __compare__( r, d );
   if ( comp <  0 ){ return; } // q = 0 and n = r < d
   if ( comp == 0 ){ // q = 1 and r = 0
      q[ 0 ] = 1;
      lq = 1;
      __clear__( r );
      lr = 0;
      return;
   }

   assert( lr >= ld );
   assert( ld >= 1  );

   const int shift = lr - ld;

   for ( int iq = shift; iq >= 0; iq-- ){

      const int ir = ld - 1 + iq;

      ucarry_t num;
      if ( ir < lr - 1 ){
         num = r[ ir + 1 ];
         num = ( num << BLOCK_BIT ) + r[ ir ];
      } else {
         num = r[ ir ];
      }
      ubase_t q_max = num / ( d[ ld - 1 ] );
      ubase_t q_min = num / ( d[ ld - 1 ] + 1UL );

      ubase_t temp[ NUM_BLOCK ];
      __mult2set__( temp, d, ld, q_max, iq );
      int comp = __compare__( r, temp );
      if ( comp >= 0 ){ q_min = q_max; }

      while ( q_max > q_min + 1 ){
         ubase_t q_test = ( ( ( ucarry_t ) q_max ) + ( ( ucarry_t ) q_min ) ) / 2;  // while condition implies q_min + 1 <= q_test <= q_max - 1
         __mult2set__( temp, d, ld, q_test, iq );
         comp = __compare__( r, temp );
         if ( comp >= 0 ){ q_min = q_test; }
                    else { q_max = q_test; }
      }

      q[ iq ] = q_min; // q_min contains solution
      if ( ( q_min > 0 ) && ( lq == 0 ) ){ lq = iq + 1; }
      if ( comp == 0 ){
         lr = 0;
         __clear__( r );
         return;
      } else {
         if ( q_min > 0 ){
            if ( comp < 0 ){ __mult2set__( temp, d, ld, q_min, iq ); } // Latest temp based on q_max if comp < 0 --> calc temp
            lr = __diff3set__( r, r, temp );
         }
      }

   }

}

int jans::big_int::__gcd__( ubase_t * res, ubase_t * a, const int la, ubase_t * b, const int lb ){

   // Solves for res = gcd( a, b ); a >= b; destroys a & b in the proces
   // Euclidean algorithm

   ubase_t * al = a; int ll = la;
   ubase_t * as = b; int ls = lb;
   int lq = 0;

   while ( ls != 0 ){

      __divide__( res, lq, al, ll, as, ls ); // al(in) = temp * as + al(out)
      ubase_t * a_swap = al; al = as; as = a_swap;
            int l_swap = ll; ll = ls; ls = l_swap;
   }

   __copy__( res, al );
   return ll;

}

int jans::big_int::__ceil_sqrt__( ubase_t * d, ubase_t * num, const int ln ){

   // Solves for d = ceil( sqrt( num ) )
   // Babylonian method: x_new = ( x + N / x ) / 2

   ubase_t q[ NUM_BLOCK ]; __clear__( q ); int lq = 0;
   ubase_t r[ NUM_BLOCK ]; __clear__( r ); int lr = 0;
                           __clear__( d ); int ld = 0;

   // Guess the sqrt and put it in (d, ld)
   {
      int lb = 0;
      for ( int j = 0; j < BLOCK_BIT; j++ ){
         if ( ( num[ j ] >> j ) & 1UL ){ lb = j; }
      }
      lb = ( lb + BLOCK_BIT * ( ln - 1 ) ) / 2;
      int id = lb / BLOCK_BIT;
      int jd = lb - BLOCK_BIT * id;
      d[ id ] = ( 1UL << jd );
      ld = id + 1;
   }

   while ( true ){

      lq = __mult3set__( q, d, ld, d, ld );     // q = d^2
      const int comp1 = __compare__( q, num );
      if ( comp1 == 0 ){ return ld; }
      lq = __sum1__( q, lq, 1 );                // q = d^2 + 1
      lr = __mult2set__( r, d, ld, 2, 0 );      // r = 2 * d

      if ( comp1 > 0 ){                         // d^2 > num
         lq = __diff3set__( q, q, r );          // q = ( d - 1 )^2
         const int comp2 = __compare__( q, num );
         if ( comp2 < 0 ){ return ld; }         // d^2 > num > ( d - 1 )^2
         if ( comp2 == 0 ){                     // d^2 > num = ( d - 1 )^2
            ld = __diff1__( d, 1 );
            return ld;
         }
      }

      if ( comp1 < 0 ){                         // d^2 < num
         __sum3set__( q, q, lq, r, lr );        // q = ( d + 1 )^2
         const int comp2 = __compare__( q, num );
         if ( comp2 >= 0 ){                     // ( d + 1 )^2 >= num > d^2
            ld = __sum1__( d, ld, 1 );
            return ld;
         }
      }

      __copy__( r, num );
      lr = ln;
      __divide__( q, lq, r, lr, d, ld );        // num = q * d + r
      ld = __sum3set__( d, d, ld, q, lq );      // d = d_old + num / d_old
      __shift_down__( d, 1 );                   // d_new = ( d_old + num / d_old ) / 2
      ld = ( ( d[ ld - 1 ] > 0 ) ? ld : ( ld - 1 ) );

   }

}


