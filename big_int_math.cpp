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
#include <stdio.h>
#include <iostream>

int jans::big_int::__sum3set__( ubase_t * r, ubase_t * a, const int la, ubase_t * b, const int lb, const int start ){

   // r[ start : ] = a[ start : ] + b[ start : ]

   const int upper = ( ( la > lb ) ? la : lb );

   ucarry_t x = 0;
   ucarry_t y = 0;
   ucarry_t z = 0;

   for ( int i = start; i < upper; i++ ){
      x = a[ i ];
      y = b[ i ];
      z = z + ( x + y );
      r[ i ] = z & __11111111__;
      z = z >> BLOCK_BIT;
   }

   if ( z != 0 ){ // Normally z should be one in this case
      if ( upper < NUM_BLOCK ){
         r[ upper ] = 1; // z & __11111111__;
         return ( upper + 1 );
      } else {
         assert( false ); // Overflow exception
      }
   }

   return upper;

}

int jans::big_int::__diff3set__( ubase_t * r, ubase_t * a, ubase_t * b ){

   // r[ : ] = a[ : ] - b[ : ]; safe to set one of a or b equal to r, requires a > b

   const int comp = __compare__( a, b );
   if ( comp == 0 ){
      __clear__( r );
      return 0;
   }

   assert( comp > 0 );

   ucarry_t add = 0;
   ucarry_t sub = 0; // serves as carry

   for ( int i = 0; i < comp; i++ ){
      add = a[ i ];
      sub = b[ i ] + sub;
      if ( add >= sub ){
         r[ i ] = ( add - sub );
         sub = 0;
      } else { // base > sub > add >= 0 --> base + add - sub > 0
         r[ i ] = ( ( add + ( ( ucarry_t )( 1 ) << BLOCK_BIT ) ) - sub );
         sub = 1;
      }
   }

   assert( sub == 0 );
   for ( int i = comp; i < NUM_BLOCK; i++ ){ r[ i ] = 0; }

   for ( int i = comp; i > 0; i-- ){ if ( r[ i - 1 ] != 0 ){ return i; } }
   return 0;

}

int jans::big_int::__mult3add__( ubase_t * r, ubase_t * a, const int la, ubase_t * b, const int lb ){

   // r[ : ] += a[ : ] * b[ : ]

   assert( ( la + lb - 2 ) < NUM_BLOCK ); // ia + ib <= la + lb - 2 < NUM_BLOCK: Overflow exception

   int lr = 0;
   for ( int ib = 0; ib < lb; ib++ ){
      lr = jans::big_int::__mult2add__( r, lr, a, la, b[ ib ], ib );
   }
   return lr;

}

int jans::big_int::__mult2set__( ubase_t * r, ubase_t * a, const int la, const ubase_t b, const int shift ){

   // r[ shift : ] = b * a[ : ]; Safe for "scal" operations when ( shift == 0 )

   const int upper = shift + la;
   assert( ( upper - 1 ) < NUM_BLOCK ); // shift + ia <= shift + la - 1 < NUM_BLOCK: Overflow exception

   ucarry_t x = 0;
   ucarry_t y = b;
   ucarry_t z = 0;

   for ( int ia = 0; ia < la; ia++ ){
      x = a[ ia ];
      z = z + ( x * y );
      r[ shift + ia ] = z & __11111111__;
      z = z >> BLOCK_BIT;
   }

   if ( z != 0 ){
      if ( upper < NUM_BLOCK ){
         r[ upper ] = z & __11111111__;
         return ( upper + 1 );
      } else {
         assert( false ); // Overflow exception
      }
   }

   return upper;

}

int jans::big_int::__mult2add__( ubase_t * r, const int lr, ubase_t * a, const int la, const ubase_t b, const int shift ){

   // r[ shift : ] += b * a[ : ]; This would be lapack "axpy" with a shift

   const int upper = shift + la;
   assert( ( upper - 1 ) < NUM_BLOCK ); // shift + ia <= shift + la - 1 < NUM_BLOCK: Overflow exception

   ucarry_t w = 0;
   ucarry_t x = 0;
   ucarry_t y = b;
   ucarry_t z = 0;

   for ( int ia = 0; ia < la; ia++ ){
      w = r[ shift + ia ];
      x = a[ ia ];
      z = z + w + ( x * y );
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

   if ( z != 0 ){ // Then while loop stopped because ( ir == NUM_BLOCK )
      assert( false );
   }

   return ( ( ir > lr ) ? ir : lr );

}

/*ubase_t jans::big_int::__div_helper__( const ubase_t n1, const ubase_t n0, const ubase_t d0 ){

   ucarry_t val = n1;

   val = ( ( val << BLOCK_BIT ) + n0 ) / d0; // val = ( n1 * b + n0 ) / d0 = q1 * b + q0
   val =   ( val >> BLOCK_BIT );             // q1 = q / b

   return ( ubase_t )( val );

}

void jans::big_int::__divide__( ubase_t * q, int & lq, ubase_t * temp, ubase_t * r, int & lr, ubase_t * d, const int ld ){

   // Solves for n = q * d + r, with r < d, whereby initially (r, lr) contains (n, ln).

   __clear__( q ); lq = 0;

   const int comp = __compare__( r, d );
   if ( comp <  0 ){ return; } // q = 0 and r = n < d
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
   if ( shift > 0 ){
      for ( int i = shift - 1 + ld; i >= shift; i-- ){ d[ i ] = d[ i - shift ]; }
      for ( int i = shift - 1;      i >= 0;     i-- ){ d[ i ] = 0; }
   }

   for ( int count = 0; count < shift; count++ ){

      ubase_t q_guess = 0;
      //if ( lr >= 2 ){
         std::cout << "Current d = " << d[ lr - 1 ] << std::endl;
         q_guess = __div_helper__( r[ lr - 1 ], r[ lr - 2 ], d[ lr - 1 ] );
      //} else {

      //}

      bool stop = false;
      while( !stop ){

         __clear__( temp );
         __mult2set__( temp, d, ld + shift - count, q_guess ); // temp = q_guess * shifted( d )
         const int check1 = __compare__( r, temp );
         assert( check1 <= lr );

         if ( check1 <  0 ){ q_guess -= 1; } // r < q_guess * shifted( d )
         if ( check1 == 0 ){ // r == q_guess * shifted( d )   ==>   done
            q[ shift - count ] = q_guess;
            __clear__( r );
            for ( int i = 0; i < ld; i++ ){ d[ i ] = d[ i + shift - count ]; }
            count = shift;
            stop = true;
            lr = 0;
         }
         if ( check1 >  0 ){ // r > q_guess * shifted( d )
            __diff3set__( temp, r, temp ); // temp = r - q_guess * shifted( d )
            const int check2 = __compare__( temp, d );
            if ( check2 >  0 ){ q_guess += 1; } // temp > shifted( d )
            if ( check2 == 0 ){ // temp == shifted( d )
               q[ shift - count ] = q_guess + 1;
               __clear__( r );
               for ( int i = 0; i < ld; i++ ){ d[ i ] = d[ i + shift - count ]; }
               count = shift;
               stop = true;
               lr = 0;
            }
            if ( check2 < 0 ){ // temp < shifted( d )
               q[ shift - count ] = q_guess;
               for ( int i = 0; i < lr; i++ ){ r[ i ] = temp[ i ]; }
               for ( int i = 0; i < ld + shift - count - 1; i++ ){ d[ i ] = d[ i + 1 ]; }
               stop = true;
               lr = lr - 1;
            }
         }

      }
   }

   if ( ld == 1 ){
      q[ 0 ] = r[ 0 ] / d[ 0 ];
      r[ 0 ] = r[ 0 ] - q[ 0 ] * d[ 0 ];
   }
   if ( ld == 2 ){
      ucarry_t num = r[ 1 ]; num = ( ( num << BLOCK_BIT ) + r[ 0 ] );
      ucarry_t div = d[ 1 ]; div = ( ( div << BLOCK_BIT ) + d[ 0 ] );
      q[ 0 ] = num / div;
      num = num - ( div * q[ 0 ] );
      r[ 1 ] = ( num >> BLOCK_BIT );
      r[ 0 ] = ( num & __11111111__ );
   }
   if ( ld >= 2 ){ // Only option left

      ubase_t q_guess = __div_helper__( r[ ld - 1 ], r[ ld - 2 ], d[ ld - 1 ] );
      bool stop = false;

      while( !stop ){
         __clear__( temp );
         __mult2set__( temp, d, ld, q_guess, 0 ); // t[ : ] = q_guess * d[ : ]
         const int check1 = __compare__( r, temp );
         if ( check1 > 0 ){ // num(rest) > q_guess * d
            __diff3set__( temp, r, temp ); // t now contains remainder = num(rest) - q_guess*d
            const int check2 = __compare__( temp, d );
            if ( check2 >= 0 ){ // remainder >= d
               q_guess += 1;
               std::cout << "__divide__ : q_guess++ for final" << std::endl;
               if ( check2 == 0 ){ // remainder == d
                  q[ 0 ] = q_guess;
                  __clear__( r );
                  stop = true;
               }
            }
            if ( check2 < 0 ){ // remainder < d: stop
               q[ 0 ] = q_guess;
               for ( int id = 0; id < ld; id++ ){ r[ id ] = temp[ id ]; }
               stop = true;
            }
         }
         if ( check1 == 0 ){ // num(rest) == q_guess * d: done
            q[ 0 ] = q_guess;
            __clear__( r );
            stop = true;
         }
         if ( check1 < 0 ){ // num(rest) < q_guess * d
            q_guess -= 1;
            std::cout << "__divide__ : q_guess-- for final" << std::endl;
         }
      }
   }

   for ( int i = 0; i < NUM_BLOCK; i++ ){ if ( q[ i ] != 0 ){ lq = ( i + 1 ); } }
   for ( int i = 0; i < NUM_BLOCK; i++ ){ if ( r[ i ] != 0 ){ lr = ( i + 1 ); } }

}*/


