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

jans::sieve::sieve( jans::big_int & num, const ubase_t B, const ubase_t M, const int extra ){

   //check_bounds_M
   assert( M >= B );
   this->M = M; // Sieve for x in [-M, M]
   target.copy( num );

   __startup1__(); // Sets initial mpqs_p
   __startup2__( B ); // Sets num_primes & creates primes, roots, and logvals

   extra_sz = extra;
   powspace = num_primes + 1; // Positive and negative Q(x)
   linspace = powspace + extra_sz;
   lincount = 0;
   xvalues  = new jans::big_int[ linspace ];
   powers   = new ubase_t[ linspace * powspace ];
   for ( int cnt = 0; cnt < ( linspace * powspace ); cnt++ ){ powers[ cnt ] = 0; }

}

jans::sieve::~sieve(){

   delete [] primes;
   delete [] roots;
   delete [] logval;

   delete [] xvalues;
   delete [] powers;

}

void jans::sieve::run( jans::big_int & p, jans::big_int & q, const double grace ){

   jans::big_int a;
   jans::big_int b;

   __next_mpqs_p__( a, b ); // 0 <= b < a /2 

   // compute inverse of a modulo p_i for factor base
   // store a x + b in list instead of x, store additionally p in list
   // check bounds for ( ax + b ) --> with a * x = sqrt( 2 N )/ M * M --> upper limit for a...

   //return;

   __sieving_grace__( M, grace );
   assert( lincount == linspace );

   unsigned char * helper = new unsigned char[ extra_sz * linspace ];

   __solve_gaussian__( helper, powers, powspace, linspace );
   __factor__( helper, p, q );

   delete [] helper;

}

void jans::sieve::__next_mpqs_p__( jans::big_int & a, jans::big_int & b ){

   bool okprime = false;

   ubase_t num_attempts = 0;

   while ( okprime == false ){

      num_attempts += 1;
      jans::big_int::minus( mpqs_p, 4 ); // To keep p % 4 == 3  and  p * p <= sqrt(2N)/M
      okprime = true;

      // Check 1: factor base does not divide mpqs_p --> if not probably prime then
      for ( int ip = 0; ip < num_primes; ip++ ){
         const ubase_t rem = jans::big_int::div( a, mpqs_p, primes[ ip ] );
         if ( rem == 0 ){
            okprime = false;
            ip = num_primes;
         }
      }

      // Check 2: (n/p) == 1
      if ( okprime == true ){
         const int symbol = __legendre_symbol__( target, mpqs_p );
         if ( symbol != 1 ){
            okprime = false;
         }
      }

      // Check 3: calculate b under assumption mpqs_p prime and check b * b == n ( mod mpqs_p * mpqs_p )
      if ( okprime == true ){

         std::cout << "Candidate p = " << mpqs_p.write( 10 ) << std::endl;

         jans::big_int work1;
         jans::big_int work2;
         jans::big_int work3;
         jans::big_int h0;
         jans::big_int h1;
         jans::big_int h2;

         // h0 = n ^ ( ( p - 3 ) / 4 ) mod p
         jans::big_int::diff( work1, mpqs_p, 3 ); // work1 = p - 3
         ubase_t rem = jans::big_int::div( work2, work1, 4 ); // work2 = ( p - 3 ) / 4
         assert( rem == 0 );
         jans::big_int::power( h0, target, work2, mpqs_p );

         // h1 = n ^ ( ( p + 1 ) / 4 ) mod p = n * h0 mod p
         jans::big_int::prodmod( work2, h1, target, h0, mpqs_p );

         // h2 = (2 * h1)^{-1}( ( n - h1^2 ) / p ) mod p
         jans::big_int::sum( work2, mpqs_p, 1 );                    // work2 = p + 1
         rem = jans::big_int::div( work3, work2, 2 );               // work3 = ( p + 1 ) / 2 : work3 = 2^{-1} mod p
         assert( rem == 0 );
         jans::big_int::prodmod( work2, work1, work3, h0, mpqs_p ); // work1 = (2 * h1)^{-1} mod p
         jans::big_int::prod( work2, h1, h1 );
         jans::big_int::diff( work3, target, work2 );
         jans::big_int::div( work2, h2, work3, mpqs_p );            // work2 = ( n - h1^2 ) / p
         assert( jans::big_int::equal( h2, 0 ) );
         jans::big_int::prodmod( work3, h2, work1, work2, mpqs_p ); // h2 = [ (2 * h1)^{-1} * ( n - h1^2 ) / p ] mod p

         // a = p * p
         jans::big_int::prod( a, mpqs_p, mpqs_p );

         // b = h1 + p*h2 <= ( p - 1 ) + p * ( p - 1 ) = p * p - 1 < a
         jans::big_int::prod( work3, h2, mpqs_p );
         jans::big_int::sum( b, h1, work3 );
         jans::diff::( work2, a, b );
         if ( jans::big_int::smaller( work2, b ) ){
            b.copy( work2 );
         }

         // check b * b % a == target % a
         jans::big_int::prod( work1, b, b );
         jans::big_int::div( work3, work2, work1,  a ); // work2 = ( b * b ) % a
         jans::big_int::div( work3, work1, target, a ); // work1 = target % a
         okprime = jans::big_int::equal( work1, work2 );
         assert( okprime );

      }
   }

}

void jans::sieve::__factor__( unsigned char * helper, jans::big_int & p, jans::big_int & q ){

   int sol = 0;

   jans::big_int x;
   jans::big_int y;
   jans::big_int work1;
   jans::big_int work2;

   while ( sol < extra_sz ){

      y.copy( 1 );

      for ( int ip = 0; ip < num_primes; ip++ ){
         ubase_t pow = 0;
         for ( int vec = 0; vec < linspace; vec++ ){
            if ( helper[ vec + linspace * sol ] == 1 ){ pow += powers[ 1 + ip + powspace * vec ]; }
         }
         assert( pow % 2 == 0 );
         pow = pow / 2;
         for ( ubase_t it = 0; it < pow; it++ ){
            jans::big_int::prod( work1, y, primes[ ip ] );
            jans::big_int::div( work2, y, work1, target );
         }
      }

      x.copy( 1 );

      for ( int vec = 0; vec < linspace; vec++ ){
         if ( helper[ vec + linspace * sol ] == 1 ){
            jans::big_int::prod( work1, x, xvalues[ vec ] );
            jans::big_int::div( work2, x, work1, target );
         }
      }

      jans::big_int::diff( work2, target, y ); // work2 = -y mod N = N - y
      if ( ( jans::big_int::equal( x, y ) == false ) && ( jans::big_int::equal( x, work2 ) == false ) ){

         if ( jans::big_int::smaller( x, y ) ){
            jans::big_int::diff( work1, y, x );
            jans::big_int::gcd( work2, work1, target );
         } else {
            jans::big_int::diff( work1, x, y );
            jans::big_int::gcd( work2, work1, target );
         }

         jans::big_int::div( work1, x, target, work2 );
         assert( jans::big_int::equal( x, 0 ) ); // Remainder zero

         if ( ! ( jans::big_int::equal( work1, 1 ) || ( jans::big_int::equal( work2, 1 ) ) ) ){

            p.copy( work2 );
            q.copy( work1 );
            return;

         }
      }

      sol++;
   }

   std::cout << "Increase number of additional smooth ( x * x - N ) to be found!" << std::endl;

}

void jans::sieve::__sieving_grace__( const ubase_t blk_size, const double grace ){

   jans::big_int work1; jans::big_int::prod( work1, target, 2 ); // 2N
   jans::big_int work2;
   jans::big_int lower_up; jans::big_int::ceil_sqrt( lower_up, target ); // ceil( sqrt( N ) )
   jans::big_int upper_up;
   jans::big_int limit_up; jans::big_int::ceil_sqrt( limit_up, work1 ); // ceil( sqrt( 2N ) )
   jans::big_int upper_dn; upper_dn.copy( lower_up ); jans::big_int::minus( upper_dn, 1 ); // ceil( sqrt( N ) ) - 1
   jans::big_int lower_dn;

   double  * sumlog    = new  double[ blk_size   ];
   ubase_t * shift1_up = new ubase_t[ num_primes ];
   ubase_t * shift2_up = new ubase_t[ num_primes ];
   ubase_t * shift1_dn = new ubase_t[ num_primes ];
   ubase_t * shift2_dn = new ubase_t[ num_primes ];
   ubase_t * helper    = new ubase_t[ num_primes ];

   for ( int ip = 0; ip < num_primes; ip++ ){
      ubase_t pri = primes[ ip ];
      ubase_t rem = jans::big_int::div( work1, lower_up, pri ); // ceil( sqrt( N ) ) % p
      ubase_t pos =     pri + roots[ ip ] - rem; while( pos >= pri ){ pos -= pri; }
      ubase_t neg = 2 * pri - roots[ ip ] - rem; while( neg >= pri ){ neg -= pri; }
      shift1_up[ ip ] = pos;
      shift2_up[ ip ] = neg;
      shift1_dn[ ip ] = pri - pos - 1;
      shift2_dn[ ip ] = pri - neg - 1;
   }

   jans::big_int::sum( upper_up, lower_up, blk_size );
   if ( jans::big_int::smaller( upper_up, limit_up ) == false ){ upper_up.copy( limit_up ); }

   if ( jans::big_int::smaller( upper_dn, blk_size ) ){ lower_dn.copy( 0 ); }
   else { jans::big_int::diff( lower_dn, upper_dn, blk_size ); }

   while ( ( lincount < linspace ) && ( ( jans::big_int::equal( lower_up, upper_up ) == false ) || ( jans::big_int::equal( lower_dn, upper_dn ) == false ) ) ){

      if ( jans::big_int::equal( lower_up, upper_up ) == false ){

         jans::big_int::diff( work1, upper_up, lower_up );
         const ubase_t loopsize = work1.get_blk( 0 );

         __fill_sumlog__( sumlog, shift1_up, shift2_up, loopsize );

         int count1 = 0;
         int count2 = 0;
         work1.copy( lower_up );
         for ( ubase_t cnt = 0; cnt < loopsize; cnt++ ){
            jans::big_int::xx_min_num( work2, work1, target );
            const double reference = jans::big_int::logarithm( work2 ) - grace;
            if ( sumlog[ cnt ] > reference ){
               count1++;
               const bool smooth = __extract__( work2, helper );
               if ( smooth ){
                  count2++;
                  xvalues[ lincount ].copy( work1 );
                  std::cout << "Smooth ( x * x - N ) no. " << lincount << " for value x = " << xvalues[ lincount ].write( 10 ) << std::endl;
                  powers[ lincount * powspace + 0 ] = 0; // Because positive Q(x)
                  for ( int ip = 0; ip < num_primes; ip++ ){
                     powers[ lincount * powspace + 1 + ip ] = helper[ ip ];
                  }
                  lincount++;
                  if ( lincount == linspace ){ cnt = loopsize; }
               }
            }
            jans::big_int::plus( work1, 1 );
         }
         std::cout << "In the present block (up) " << count1 << " / " << loopsize << " were selected based on sums of logarithms." << std::endl;
         std::cout << "In the present block (up) " << count2 << " / " << count1   << " of the candidates were smooth." << std::endl;

         lower_up.copy( upper_up );
         jans::big_int::sum( upper_up, lower_up, blk_size );
         if ( jans::big_int::smaller( upper_up, limit_up ) == false ){ upper_up.copy( limit_up ); }

      }

      if ( ( lincount < linspace ) && ( jans::big_int::equal( lower_dn, upper_dn ) == false ) ){

         jans::big_int::diff( work1, upper_dn, lower_dn );
         const ubase_t loopsize = work1.get_blk( 0 );

         __fill_sumlog__( sumlog, shift1_dn, shift2_dn, loopsize );

         int count1 = 0;
         int count2 = 0;
         work1.copy( upper_dn );
         for ( ubase_t cnt = 0; cnt < loopsize; cnt++ ){
            jans::big_int::min_xx_plus_num( work2, work1, target );
            const double reference = jans::big_int::logarithm( work2 ) - grace;
            if ( sumlog[ cnt ] > reference ){
               count1++;
               const bool smooth = __extract__( work2, helper );
               if ( smooth ){
                  count2++;
                  xvalues[ lincount ].copy( work1 );
                  std::cout << "Smooth ( x * x - N ) no. " << lincount << " for value x = " << xvalues[ lincount ].write( 10 ) << std::endl;
                  powers[ lincount * powspace + 0 ] = 1; // Because negative Q(x)
                  for ( int ip = 0; ip < num_primes; ip++ ){
                     powers[ lincount * powspace + 1 + ip ] = helper[ ip ];
                  }
                  lincount++;
                  if ( lincount == linspace ){ cnt = loopsize; }
               }
            }
            jans::big_int::minus( work1, 1 );
         }
         std::cout << "In the present block (dn) " << count1 << " / " << loopsize << " were selected based on sums of logarithms." << std::endl;
         std::cout << "In the present block (dn) " << count2 << " / " << count1   << " of the candidates were smooth." << std::endl;

         upper_dn.copy( lower_dn );
         if ( jans::big_int::smaller( upper_dn, blk_size ) ){ lower_dn.copy( 0 ); }
         else { jans::big_int::diff( lower_dn, upper_dn, blk_size ); }

      }

   }

   delete [] sumlog;
   delete [] shift1_up;
   delete [] shift2_up;
   delete [] shift1_dn;
   delete [] shift2_dn;
   delete [] helper;

}

void jans::sieve::__fill_sumlog__( double * sumlog, ubase_t * shift1, ubase_t * shift2, const ubase_t loopsize ) const{

   for ( ubase_t cnt = 0; cnt < loopsize; cnt++ ){ sumlog[ cnt ] = 0.0; }

   for ( int ip = 0; ip < num_primes; ip++ ){
      const ubase_t prime = primes[ ip ];
      const double  log_p = logval[ ip ];
            ubase_t index = shift1[ ip ];
      while ( index < loopsize ){
         sumlog[ index ] += log_p;
         index += prime;
      }
      shift1[ ip ] = index - loopsize;
   }

   for ( int ip = 1; ip < num_primes; ip++ ){ // Not for prime 2
      const ubase_t prime = primes[ ip ];
      const double  log_p = logval[ ip ];
            ubase_t index = shift2[ ip ];
      while ( index < loopsize ){
         sumlog[ index ] += log_p;
         index += prime;
      }
      shift2[ ip ] = index - loopsize;
   }

}

