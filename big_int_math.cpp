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

int jans::big_int::__plus_one__( ubase_t * r, const int lr ){

   // r++

   ucarry_t z = 1;

   int ir = 0;
   while ( ( ir < NUM_BLOCK ) && ( z != 0 ) ){
      z = z + r[ ir ];
      r[ ir ] = z & __11111111__;
      z = z >> BLOCK_BIT;
      ir++;
   }

   if ( z != 0 ){ assert( false ); } // Then while loop stopped because ( ir == NUM_BLOCK ): Overflow exception

   return ( ( ir > lr ) ? ir : lr );

}

int jans::big_int::__diff3set__( ubase_t * r, ubase_t * a, ubase_t * b ){

   // r = a - b

   // if ( ( r != a ) && ( r != b ) ){ __clear__( r ); } ---> Not necessary: for i < comp, it will be overwritten; for i >= comp, it will be set to zero.
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
   }
   for ( int i = comp; i < NUM_BLOCK; i++ ){ r[ i ] = 0; }

   assert( sub == 0 );
   for ( int i = comp; i > 0; i-- ){ if ( r[ i - 1 ] != 0 ){ return i; } }
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

int jans::big_int::__mult2set__( ubase_t * r, ubase_t * a, const int la, const ubase_t b ){

   // r = a * b

   __clear__( r );
   const int lr = __mult2add__( r, 0, a, la, b, 0 );
   return lr;

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

   if ( z != 0 ){ assert( false ); } // Then while loop stopped because ( ir == NUM_BLOCK ): Overflow exception

   return ( ( ir > lr ) ? ir : lr );

}

int jans::big_int::__scal1__( ubase_t * r, const int lr, const ubase_t b ){

   // r = r * b

   ucarry_t x = 0;
   ucarry_t z = 0;

   for ( int ir = 0; ir < lr; ir++ ){
      x = r[ ir ];
      z = z + ( x * b );
      r[ ir ] = z & __11111111__;
      z = z >> BLOCK_BIT;
   }

   if ( z != 0 ){
      if ( lr < NUM_BLOCK ){
         r[ lr ] = z & __11111111__;
         return ( lr + 1 );
      } else {
         assert( false ); // Overflow exception
      }
   }

   return lr;

}

void jans::big_int::__divide__( ubase_t * q, int & lq, ubase_t * r, int & lr, ubase_t * d, const int ld ){

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

   //int jrl = 0; { int j = 0; while( ( j < BLOCK_BIT ) && ( jrl == 0 ) ){ 
   int jrl = 0; for ( int j = 0; j < BLOCK_BIT; j++ ){ if ( ( r[ lr - 1 ] >> j ) & 1U ){ jrl = j; } }
   int jdl = 0; for ( int j = 0; j < BLOCK_BIT; j++ ){ if ( ( d[ ld - 1 ] >> j ) & 1U ){ jdl = j; } }
   const int shift_bit = ( lr - ld ) * BLOCK_BIT + ( jrl - jdl );

   std::cout << "lr  = " << lr << std::endl;
   std::cout << "ld  = " << ld << std::endl;
   std::cout << "jrl = " << jrl << std::endl;
   std::cout << "jdl = " << jdl << std::endl;
   std::cout << "shift_bit = " << shift_bit << std::endl;

   assert( shift_bit >= 0 );

   if ( shift_bit > 0 ){ __shift_up__( d, shift_bit ); }

   for ( int jq = shift_bit; jq >= 0; jq-- ){

      const int larger = __compare__( r, d );
      if ( larger >= 0 ){
         lr = __diff3set__( r, r, d );
         q[ jq / BLOCK_BIT ] |= ( 1U << ( jq % BLOCK_BIT ) );
         if ( lq == 0 ){ lq = ( jq / BLOCK_BIT ) + 1; }
      }
      if ( jq > 0 ){ __shift_down__( d, 1 ); }

   }

}


