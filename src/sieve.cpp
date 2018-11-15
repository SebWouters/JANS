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

jans::sieve::sieve( jans::big_int & num, const ubase_t factorbound, const ubase_t sievespace, const int congruences ){

   //check_bounds_M
   assert( sievespace >= factorbound );
   this->M = sievespace; // Sieve for x in [-M, M]
   target.copy( num );

   __startup1__(); // Sets initial mpqs_q
   __startup2__( factorbound ); // Sets num_primes & creates primes, roots, and logvals

   extra    = congruences;
   powspace = num_primes + 1; // Positive and negative Q(x)
   linspace = powspace + extra;

   // Don't create overflows
   assert( ( ( ( ucarry_t )( linspace ) ) * ( ( ucarry_t )( powspace ) ) ) <= ( INT_MAX ) );

   lincount = 0;
   xvalues  = new jans::big_int[ linspace ];
   pvalues  = new jans::big_int[ linspace ];
   powers   = new ubase_t[ linspace * powspace ];
   for ( int cnt = 0; cnt < ( linspace * powspace ); cnt++ ){ powers[ cnt ] = 0; }

}

jans::sieve::~sieve(){

   delete [] primes;
   delete [] roots;
   delete [] logval;

   delete [] xvalues;
   delete [] pvalues;
   delete [] powers;

}

void jans::sieve::run( jans::big_int & sol_p, jans::big_int & sol_q, const double threshold ){

   jans::big_int a;
   jans::big_int b;

   const ubase_t size = 2 * M + 1;
   double  * sumlog = new double[ size ];
   ubase_t * shift1 = new ubase_t[ num_primes ];
   ubase_t * shift2 = new ubase_t[ num_primes ];

   while ( lincount < linspace ){
      __next_mpqs_q__( a, b ); // 0 <= b < a /2
      __calculate_shifts__( shift1, shift2, a, b );
      __sieve_sumlog__( size, sumlog, shift1, shift2 );
      __check_sumlog__( size, sumlog, shift1, threshold, a, b );
   }

   delete [] shift1;
   delete [] shift2;
   delete [] sumlog;

   unsigned char * helper = new unsigned char[ extra * linspace ];
   __solve_gaussian__( helper, powers, powspace, linspace );
   __factor__( helper, sol_p, sol_q );
   delete [] helper;

}

void jans::sieve::__next_mpqs_q__( jans::big_int & a, jans::big_int & b ){

   bool okprime = false;

   while ( okprime == false ){

      jans::big_int::minus( mpqs_q, 4 ); // To keep p % 4 == 3  and  p * p <= sqrt(2N)/M
      okprime = true;

      // Check 1: factor base does not divide mpqs_q --> if not probably prime then
      for ( int ip = 0; ip < num_primes; ip++ ){
         const ubase_t rem = jans::big_int::div( a, mpqs_q, primes[ ip ] );
         if ( rem == 0 ){
            okprime = false;
            ip = num_primes;
         }
      }

      // Check 2: Miller Rabin
      if ( okprime == true ){
         okprime = jans::big_int::miller_rabin( mpqs_q, 20 );
      }

      // Check 3: (n/p) == 1
      if ( okprime == true ){
         const int symbol = __legendre_symbol__( target, mpqs_q );
         if ( symbol != 1 ){
            okprime = false;
         }
      }

      // Check 4: calculate b under assumption mpqs_q prime and check b * b == n ( mod mpqs_q * mpqs_q )
      if ( okprime == true ){

         jans::big_int work1;
         jans::big_int work2;
         jans::big_int work3;
         jans::big_int h0;
         jans::big_int h1;
         jans::big_int h2;

         // h0 = n ^ ( ( p - 3 ) / 4 ) mod q
         jans::big_int::diff( work1, mpqs_q, 3 ); // work1 = q - 3
         ubase_t rem = jans::big_int::div( work2, work1, 4 ); // work2 = ( q - 3 ) / 4
         jans::big_int::power( h0, target, work2, mpqs_q );

         // h1 = n ^ ( ( p + 1 ) / 4 ) mod q = n * h0 mod q
         jans::big_int::prodmod( work2, h1, target, h0, mpqs_q );

         // h2 = (2 * h1)^{-1}( ( n - h1^2 ) / q ) mod q
         jans::big_int::sum( work2, mpqs_q, 1 );                    // work2 = q + 1
         rem = jans::big_int::div( work3, work2, 2 );               // work3 = ( q + 1 ) / 2 : work3 = 2^{-1} mod q
         jans::big_int::prodmod( work2, work1, work3, h0, mpqs_q ); // work1 = (2 * h1)^{-1} mod q (h0 is inverse of h1 if q prime)
         jans::big_int::prod( work2, h1, h1 );
         jans::big_int::diff( work3, target, work2 );
         jans::big_int::div( work2, h2, work3, mpqs_q );            // work2 = ( n - h1^2 ) / q
         if ( jans::big_int::equal( h2, 0 ) == false ){
            okprime = false;
         } else {
            jans::big_int::prodmod( work3, h2, work1, work2, mpqs_q ); // h2 = [ (2 * h1)^{-1} * ( n - h1^2 ) / q ] mod q

            // a = q * q
            jans::big_int::prod( a, mpqs_q, mpqs_q );

            // b = h1 + q*h2 <= ( q - 1 ) + q * ( q - 1 ) = q * q - 1 < a
            jans::big_int::prod( work3, h2, mpqs_q );
            jans::big_int::sum( b, h1, work3 );
            jans::big_int::diff( work2, a, b );
            if ( jans::big_int::smaller( work2, b ) ){
               b.copy( work2 ); // Hence 0 <= b < a/2
            }

            // check b * b % a == target % a
            jans::big_int::prod( work1, b, b );
            jans::big_int::diff( work3, target, work1 );
            jans::big_int::div( work2, work1, work3, a ); // work2 contains abs_c = (n - b*b)/a, work1 remainder (should be zero)
            okprime = jans::big_int::equal( work1, 0 );
         }
      }
   }

}

void jans::sieve::__calculate_shifts__( ubase_t * shift1, ubase_t * shift2, jans::big_int & a, jans::big_int & b ){

   jans::big_int work;
   for ( int ip = 0; ip < num_primes; ip++ ){
      ubase_t pri = primes[ ip ];
      ubase_t inv = __inv_x_mod_p__( a, pri );
      ubase_t rmb = jans::big_int::div( work, b, pri );
      ubase_t rmm = M % pri;
      shift1[ ip ] = ( ( ( ucarry_t )( inv ) ) * ( ( ucarry_t )(     pri + roots[ ip ] - rmb ) ) + ( ( ucarry_t )( rmm ) ) ) % pri;
      shift2[ ip ] = ( ( ( ucarry_t )( inv ) ) * ( ( ucarry_t )( 2 * pri - roots[ ip ] - rmb ) ) + ( ( ucarry_t )( rmm ) ) ) % pri;
   }

}

void jans::sieve::__factor__( unsigned char * helper, jans::big_int & p, jans::big_int & q ){

   int sol = 0;

   jans::big_int x;
   jans::big_int y;
   jans::big_int work1;
   jans::big_int work2;
   jans::big_int work3;

   while ( sol < extra ){

      y.copy( 1 );

      for ( int ip = 0; ip < num_primes; ip++ ){
         ubase_t pow = 0;
         for ( int vec = 0; vec < linspace; vec++ ){
            if ( helper[ vec + linspace * sol ] == 1 ){ pow += powers[ 1 + ip + powspace * vec ]; }
         }
         assert( pow % 2 == 0 );
         pow = pow / 2;
         work2.copy( pow );
         work1.copy( primes[ ip ] );
         jans::big_int::power( work3, work1, work2, target ); // work3 = prime ^ ( pow ) % N
         jans::big_int::prod( work1, y, work3 );
         jans::big_int::div( work2, y, work1, target ); // y *= work3 % N
      }
      for ( int vec = 0; vec < linspace; vec++ ){
         if ( helper[ vec + linspace * sol ] == 1 ){
            jans::big_int::prod( work1, y, pvalues[ vec ] );
            jans::big_int::div( work2, y, work1, target ); // y *= p(a,vec) % N
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

   std::cerr << "   Error: Increase -Z, --congruences" << std::endl;

}

void jans::sieve::__check_sumlog__( const ubase_t size, double * sumlog, ubase_t * helper, const double threshold, jans::big_int & a, jans::big_int & b ){

   // 0 <= b < a/2
   // ( a * x + 2 * b ) * x is non-negative ( check with x in [-M, -1] and [0, M] resp. )
   // ( a * x + b ) is negative for x in [-M, -1] and non-negative for x in [0, M] resp. )

   jans::big_int work1;
   jans::big_int work2;
   jans::big_int abs_c;
   jans::big_int::prod( work1, b, b );
   jans::big_int::diff( work2, target, work1 );
   jans::big_int::div( abs_c, work1, work2, a );
   assert( jans::big_int::equal( work1, 0 ) );

   int cnt_sumlog = 0;
   int cnt_smooth = 0;

   ubase_t cnt = 0;
   while ( ( cnt < size ) && ( lincount < linspace ) ){

      // work1 = a * x * x + 2 * b * x >= 0
      const ubase_t abs_x = ( ( cnt < M ) ? ( M - cnt ) : ( cnt - M ) );
      jans::big_int::prod( work2, a,     abs_x );
      jans::big_int::prod( work1, work2, abs_x );
      jans::big_int::prod( work2, b, 2 * abs_x );
      if ( cnt < M ){ jans::big_int::diff( work1, work1, work2 ); }
               else { jans::big_int::sum(  work1, work1, work2 ); }

      // work2 = abs( a * x * x + 2 * b * x + c ) and negative contains sign indication
      const bool negative = jans::big_int::smaller( work1, abs_c );
      if ( negative ){ jans::big_int::diff( work2, abs_c, work1 ); }
                else { jans::big_int::diff( work2, work1, abs_c ); }

      const double reference = log( jans::big_int::i2f( work2 ) ) - threshold;
      if ( sumlog[ cnt ] > reference ){
         cnt_sumlog++;
         const bool smooth = __extract__( work2, helper );
         if ( smooth ){
            cnt_smooth++;
            jans::big_int::prod( work1, a, abs_x );
            if ( cnt < M ){ jans::big_int::diff( work1, work1, b ); }
                     else { jans::big_int::sum(  work1, work1, b ); }
            xvalues[ lincount ].copy( work1 );
            pvalues[ lincount ].copy( mpqs_q );
            powers[ lincount * powspace + 0 ] = ( ( negative ) ? 1 : 0 );
            for ( int ip = 0; ip < num_primes; ip++ ){
               powers[ lincount * powspace + 1 + ip ] = helper[ ip ];
            }
            lincount++;
         }
      }
      cnt++;
   }

   std::cout << "For q = " << mpqs_q.write( 10 ) << ", sieving retains " << cnt_sumlog << " / " << cnt
                                       << " and trial division retains " << cnt_smooth << " / " << cnt_sumlog << "." << std::endl;
   std::cout << "Obtained / required B-smooth numbers = " << lincount << " / " << linspace << "." << std::endl;

}

void jans::sieve::__sieve_sumlog__( const ubase_t size, double * sumlog, ubase_t * shift1, ubase_t * shift2 ) const{

   for ( ubase_t cnt = 0; cnt < size; cnt++ ){ sumlog[ cnt ] = 0.0; }

   for ( int ip = 0; ip < num_primes; ip++ ){
      const ubase_t prime = primes[ ip ];
      const double  log_p = logval[ ip ];
            ubase_t index = shift1[ ip ];
      while ( index < size ){
         sumlog[ index ] += log_p;
         index += prime;
      }
   }

   for ( int ip = 1; ip < num_primes; ip++ ){ // Not for prime 2
      const ubase_t prime = primes[ ip ];
      const double  log_p = logval[ ip ];
            ubase_t index = shift2[ ip ];
      while ( index < size ){
         sumlog[ index ] += log_p;
         index += prime;
      }
   }

}

ucarry_t jans::sieve::optimal_factorbound( jans::big_int & number ){

   const double   ln_val   = log( jans::big_int::i2f( number ) );
   const ucarry_t optbound = ceil( exp( 0.5 * sqrt( ln_val * log( ln_val ) ) ) );
   return optbound;

}

