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
#include <stdlib.h>

#include "sieve.h"

jans::sieve::sieve( const ubase_t bound, big_int & num ){

   this->bound = bound;
   target.copy( num );
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
            const bool ok = ( __legendre_symbol__( target, number ) == 1 );
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

            const ubase_t a = __root_quadratic_residue__( target, number );
            primes[ check ] = number;
             roots[ check ] = a;
            logval[ check ] = log( ( double ) number );
            check++;

         }
      }
   }
   assert( check == num_primes );

   delete [] helper;

}

int jans::sieve::__legendre_symbol__( big_int & num, const ubase_t p ){

   big_int work;
   const ubase_t rem = jans::big_int::div( work, num, p );
   return __legendre_symbol__( rem, p );

}

int jans::sieve::__legendre_symbol__( const ubase_t num, const ubase_t p ){

   // Algorithm 2.3.5, Crandall & Pomerance

   assert( ( p % 2 ) != 0 );

   int t = 1;
   ubase_t a = num % p;
   ubase_t m = p;
   ubase_t s = 0;

   while ( a != 0 ){
      while ( ( a % 2 ) == 0 ){
         a = a / 2;
         s = m % 8;
         if ( ( s == 3 ) || ( s == 5 ) ){ t = -t; }
      }
      s = a;
      a = m;
      m = s;
      if ( ( ( a % 4 ) == 3 ) && ( m % 4 == 3 ) ){ t = -t; }
      a = a % m;
   }
   if ( m == 1 ){ return t; }
   return 0;

}

ubase_t jans::sieve::__power__( const ubase_t num, const ubase_t pow, const ubase_t mod ){

   ucarry_t res = 1;
   ucarry_t temp = num % mod; // temp = num^{2^j} % mod
   for ( int j = 0; j < BLOCK_BIT; j++ ){
      if ( ( pow >> j ) & 1U ){ res = ( res * temp ) % mod; }
      temp = ( temp * temp ) % mod;
   }
   return res;

}

ubase_t jans::sieve::__root_quadratic_residue__( big_int & num, const ubase_t p ){

   big_int work;
   const ubase_t rem = jans::big_int::div( work, num, p ); // rem = num % p
   return __root_quadratic_residue__( rem, p );

}

ubase_t jans::sieve::__root_quadratic_residue__( const ubase_t num, const ubase_t p ){

   // Tonelli–Shanks
   // https://ipfs.io/ipfs/QmXoypizjW3WknFiJnKLwHCnL72vedxjQkDDP1mXWo6uco/wiki/Tonelli%E2%80%93Shanks_algorithm.html

   const ubase_t rem = num % p;

   if ( p % 4 == 3 ){
      return __power__( rem, ( p + 1 ) / 4, p );
   }

   ubase_t S = 0;
   ubase_t Q = p - 1;
   while ( Q % 2 == 0 ){
      S++;
      Q = Q / 2;
   } // p - 1 = 2^S * Q(odd)

   ubase_t z = 0;
   int leg_sym = 0;
   while ( leg_sym != -1 ){
      z = ( rand() % ( p - 2 ) ) + 2; // random number in [ 2 .. p - 1 ]
      leg_sym = __legendre_symbol__( z, p );
   }

   ubase_t c = __power__(   z, Q, p );
   ubase_t t = __power__( rem, Q, p );
   ubase_t M = S;
   ubase_t R = __power__( rem, ( Q + 1 ) / 2, p );

   while ( t != 1 ){

      ubase_t i = 0;
      ucarry_t b = t;
      while ( ( b != 1 ) && ( i < M ) ){
         b = ( b * b ) % p;
         i++;
      }
      assert( i < M );

      ubase_t j = 0;
      b = c;
      while ( j < M - i - 1 ){
         b = ( b * b ) % p; // c^{2^j}
         j++;
      }

      R = ( R * b ) % p;
      t = ( ( ( t * b ) % p ) * b ) % p;
      c = ( b * b ) % p;
      M = i;

   }

   ucarry_t test = R;
   test = ( test * test ) % p;
   assert( test == rem );
   return R;

}







