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
#include <math.h>

#include "sieve.h"

void jans::sieve::__startup1__(){

   jans::big_int work;
   ubase_t rem;

   jans::big_int::prod( work, target, 2 );      // work   = 2 * N
   jans::big_int::ceil_sqrt( mpqs_q, work );    // mpqs_q = ceil( sqrt( 2N ) )
   rem = jans::big_int::div( work, mpqs_q, M ); // work   = floor( ceil( sqrt( 2N ) ) / M )
   jans::big_int::ceil_sqrt( mpqs_q, work );    // mpqs_q = ceil( sqrt( floor( ceil( sqrt( 2N ) ) / M ) ) )

   rem = jans::big_int::div( work, mpqs_q, 4 );
   jans::big_int::minus( mpqs_q, 1 + rem ); // mpqs_q % 4 = 3

}

void jans::sieve::__startup2__( const ubase_t bound ){

   /*
    *   Map: number = 2 * cnt + 1 <= bound = 2 * floor( bound / 2 )     = 2 * floor( ( bound - 1 ) / 2 ) + 1 (bound even)
    *                                      = 2 * floor( bound / 2 ) + 1 = 2 * floor( ( bound - 1 ) / 2 ) + 1 (bound odd)
    */

   const ubase_t cntmax = ( bound - 1 ) / 2;
   unsigned char * helper = new unsigned char[ cntmax + 1 ];
   for ( ubase_t cnt = 0; cnt <= cntmax; cnt++ ){ helper[ cnt ] = 1; }

   num_primes = 1; // 2 is a prime, even numbers not looped

   for ( ubase_t cnt = 1; cnt <= cntmax; cnt++ ){
      if ( helper[ cnt ] == 1 ){
         const ubase_t number = 2 * cnt + 1;
         const bool    ok     = ( __legendre_symbol__( target, number ) == 1 );
         const ubase_t start  = ( ( ok ) ? 3 : 1 );
         const ubase_t stop   = bound / number;
         if ( ok ){ num_primes++; }
         for ( ubase_t factor = start; factor <= stop; factor += 2 ){
            const ubase_t map = factor * cnt + ( factor / 2 );
            helper[ map ] = 0;
         }
      }
   }

   primes = new ubase_t[ num_primes ];
   roots  = new ubase_t[ num_primes ];
   logval = new  double[ num_primes ];

   primes[ 0 ] = 2;
    roots[ 0 ] = 1;
   logval[ 0 ] = log( 2.0 );

   int check = 1;

   for ( ubase_t cnt = 1; cnt <= cntmax; cnt++ ){
      if ( helper[ cnt ] == 1 ){
         const ubase_t number = 2 * cnt + 1;
         const ubase_t a      = __root_quadratic_residue__( target, number );
              primes[ check ] = number;
               roots[ check ] = a;
              logval[ check ] = log( ( double ) number );
              check++;
      }
   }
   assert( check == num_primes );

   delete [] helper;

}

