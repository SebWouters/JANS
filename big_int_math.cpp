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

int jans::big_int::__add3kernel__( ubase_t * r, ubase_t * a, const int la, ubase_t * b, const int lb, const int start ){

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
         r[ upper ] = z & __11111111__;
         return ( upper + 1 );
      } else {
         assert( false ); // Overflow exception
      }
   }

   return upper;

}

int jans::big_int::__multiply3kernel__( ubase_t * r, ubase_t * t, ubase_t * a, const int la, ubase_t * b, const int lb ){

   // r[ : ] = a[ : ] * b[ : ]

   assert( ( la + lb - 2 ) < NUM_BLOCK ); // ia + ib <= la + lb - 2 < NUM_BLOCK: Overflow exception

   int lr = 0;
   int lt = 0;
   for ( int ib = 0; ib < lb; ib++ ){
      lt = jans::big_int::__multiply2kernel__( t, a, la, b[ ib ], ib );
      lr = jans::big_int::__add3kernel__( r, t, lt, r, lr, ib );
   }
   return lr;

}

int jans::big_int::__multiply2kernel__( ubase_t * r, ubase_t * a, const int la, const ubase_t b, const int shift ){

   // r[ shift : ] = a[ : ] * b

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

/*int jans::big_int::__daxpy2kernel__( ubase_t * r, const int lr, ubase_t * a, const int la, const ubase_t b, const int shift ){

   const int upper = shift + la
   assert( ( upper - 1 ) < NUM_BLOCK ); // shift + ia <= shift + la - 1 < NUM_BLOCK: Overflow exception

   ucarry_t z = 0;
   ucarry_t y = b;

   for ( int ia = 0; ia < la; ia++ ){
      ucarry_t w = r[ shift + ia ];
      ucarry_t x = a[ ia ];
      z = z + w + ( x * y );
      r[ shift + ia ] = z & __11111111__;
      z = z >> BLOCK_BIT;
   }

   if ( z != 0 ){
      if ( upper < NUM_BLOCK ){
         if ( lr <= upper ){ // r[ upper ] == 0
            r[ upper ] = z & __11111111__;
            return ( upper + 1 );
         } else { // r[ upper ] != 0
            // dump carry
         }
      } else {
         assert( false ); // Overflow exception
      }
   }

   return upper;

}*/


