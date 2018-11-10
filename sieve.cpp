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
#include <stdio.h>
#include <iostream>
#include <math.h>

#include "sieve.h"

jans::sieve::sieve( const ubase_t bound, big_int & num ){

   this->bound = bound;
   n.copy( num );
   __startup__();

}

jans::sieve::~sieve(){

   delete [] primes;
   delete [] roots;
   delete [] logval;

}

void jans::sieve::__startup__(){

   const int num_blk = bound / ( 2 * BLOCK_BIT ) + 1;
   assert( num_blk < ( 1UL << BLOCK_BIT ) / ( 2 * BLOCK_BIT ) ); // Check that ( 2 * BLOCK_BIT * num_blk ) does not overflow ubase_t
   assert( 2 * num_blk * BLOCK_BIT - 1 >= bound );

   ubase_t * helper = new ubase_t[ num_blk ];
   for ( int i = 0; i < num_blk; i++ ){
      helper[ i ] = __11111111__;
   }

   num_primes = 1; // 2 is a prime, even numbers not looped

   for ( ubase_t ip = 0; ip < num_blk; ip++ ){
      for ( ubase_t jp = 0; jp < BLOCK_BIT; jp++ ){

         const ubase_t number = 2 * ( BLOCK_BIT * ip + jp ) + 1;
         bool process = ( ( helper[ ip ] >> jp ) & 1U );
         process = ( ( process ) && ( number > 2 ) && ( number <= bound ) );

         if ( process ){
            const bool ok = ( __legendre_symbol__( n, number ) == 1 );
            if ( ok ){ num_primes++; }
            const ubase_t start = ( ( ok ) ? 3 : 1 );
            const ubase_t stop  = bound / number;

            for ( ubase_t factor = start; factor <= stop; factor += 2 ){ // only consider odd factors
               const ubase_t it = factor * ip + ( ( 2 * jp + 1 ) * factor - 1 ) / ( 2 * BLOCK_BIT );
               const ubase_t jt =             ( ( ( 2 * jp + 1 ) * factor - 1 ) % ( 2 * BLOCK_BIT ) ) / 2;
               helper[ it ] &= ~( 1U << jt );
            }
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

   for ( ubase_t ip = 0; ip < num_blk; ip++ ){
      for ( ubase_t jp = 0; jp < BLOCK_BIT; jp++ ){

         const ubase_t number = 2 * ( BLOCK_BIT * ip + jp ) + 1;
         bool process = ( ( helper[ ip ] >> jp ) & 1U );
         process = ( ( process ) && ( number > 2 ) && ( number <= bound ) );

         if ( process ){
            // Obtain root
            const ubase_t root = 0;

            primes[ check ] = number;
             roots[ check ] = root;
            logval[ check ] = log( ( double ) number );
            check++;

         }
      }
   }
   assert( check == num_primes );

   delete [] helper;

}

int jans::sieve::__legendre_symbol__( big_int & num, const ubase_t d_prime ){

   // Algorithm 2.3.5, Crandall & Pomerance

   assert( ( d_prime % 2 ) != 0 );

   int t = 1;
   big_int work;
   ubase_t a = jans::big_int::div( work, num, d_prime );
   ubase_t m = d_prime;

   while ( a != 0 ){
      while ( ( a % 2 ) == 0 ){
         a = a / 2;
         ubase_t temp = m % 8;
         if ( ( temp == 3 ) || ( temp == 5 ) ){ t = -t; }
      }
      ubase_t swap = a;
      a = m;
      m = swap;
      if ( ( ( a % 4 ) == 3 ) && ( m % 4 == 3 ) ){ t = -t; }
      a = a % m;
   }
   if ( m == 1 ){ return t; }
   return 0;

}


