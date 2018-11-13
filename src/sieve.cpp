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

jans::sieve::sieve( const ubase_t bound, jans::big_int & num, const int extra ){

   this->bound = bound;
   target.copy( num );

   __startup__(); // Creates primes, roots, and logvals

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

void jans::sieve::run( jans::big_int & p, jans::big_int & q, const ubase_t blk_size, const double grace ){

   __sieving_grace__( blk_size, grace );
   assert( lincount == linspace );

   unsigned char * helper = new unsigned char[ extra_sz * linspace ];

   __solve_gaussian__( helper );
   __factor__( helper, p, q );

   delete [] helper;

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

void jans::sieve::__solve_gaussian__( unsigned char * helper ) const{

   unsigned char * lin_contributes = new unsigned char[ linspace ];
   unsigned char * pow_contributes = new unsigned char[ powspace ];
   for ( ubase_t lin = 0; lin < linspace; lin++ ){ lin_contributes[ lin ] = 1; }
   for ( ubase_t pow = 0; pow < powspace; pow++ ){ pow_contributes[ pow ] = 1; }

   ubase_t vecspace = linspace;
   for ( ubase_t pow = 0; pow < powspace; pow++ ){
      ubase_t num_odds = 0;
      ubase_t last = linspace;
      for ( ubase_t lin = 0; lin < linspace; lin++ ){
         if ( ( powers[ pow + powspace * lin ] % 2 ) == 1 ){
            num_odds++;
            last = lin;
            if ( num_odds == 2 ){ lin = linspace; }
         }
      }
      if ( num_odds == 1 ){
         vecspace--;
         lin_contributes[ last ] = 0;
      }
   }

   std::cout << "Removed " << linspace - vecspace << " of the " << linspace << " vectors with a unique odd prime power." << std::endl;

   ubase_t redspace = powspace;
   for ( ubase_t pow = 0; pow < powspace; pow++ ){
      ubase_t num_odds = 0;
      ubase_t lin = 0;
      for ( ubase_t vec = 0; vec < vecspace; vec++ ){
         while ( lin_contributes[ lin ] == 0 ){ lin++; }
         if ( ( powers[ pow + powspace * lin ] % 2 ) == 1 ){
            num_odds++;
            vec = vecspace;
         }
         lin++;
      }
      if ( num_odds == 0 ){
         redspace--;
         pow_contributes[ pow ] = 0;
      }
   }

   std::cout << "Removed " << powspace - redspace << " of the " << powspace << " primes with only even powers." << std::endl;

   const ubase_t d_vec = vecspace;
   const ubase_t d_lin = linspace;
   const ubase_t d_red = redspace;
   const ubase_t d_pow = powspace;
   const ubase_t d_row = redspace + linspace;
   const ubase_t d_sol = linspace - powspace;

   unsigned char * matrix = new unsigned char[ d_row * d_vec ];

   /*
    * matrix = [  d_red x d_vec  ]
    *          [  -------------  ]
    *          [  d_lin x d_vec  ]
    */

   for ( ubase_t vec = 0; vec < d_vec; vec++ ){
      for ( ubase_t row = 0; row < d_row; row++ ){
         matrix[ row + d_row * vec ] = 0;
      }
   }

   {
      ubase_t lin = 0;
      for ( ubase_t vec = 0; vec < d_vec; vec++ ){
         while ( lin_contributes[ lin ] == 0 ){ lin++; }
         ubase_t pow = 0;
         for ( ubase_t red = 0; red < d_red; red++ ){
            while ( pow_contributes[ pow ] == 0 ){ pow++; }
            matrix[ red + d_row * vec ] = ( powers[ pow + d_pow * lin ] % 2 );
            pow++;
         }
         matrix[ d_red + lin + d_row * vec ] = 1;
         lin++;
      }
   }

   delete [] lin_contributes;
   delete [] pow_contributes;

   ubase_t start_vec = 0;
   for ( ubase_t red = 0; red < d_red; red++ ){
      bool found   = false;
      ubase_t iter = start_vec;
      while ( ( found == false ) && ( iter < d_vec ) ){
         if ( matrix[ red + d_row * iter ] == 1 ){ found = true; }
         else { iter++; }
      }
      if ( found == true ){
         if ( iter != start_vec ){
            for ( ubase_t row = red; row < d_row; row++ ){
               matrix[ row + d_row * start_vec ] = ( matrix[ row + d_row * iter ] ) ^ ( matrix[ row + d_row * start_vec ] );
            }
         }
         for ( ubase_t vec = start_vec + 1; vec < d_vec; vec++ ){
            if ( matrix[ red + d_row * vec ] == 1 ){
               for ( ubase_t row = red; row < d_row; row++ ){
                  matrix[ row + d_row * vec ] = ( matrix[ row + d_row * vec ] ) ^ ( matrix[ row + d_row * start_vec ] );
               }
            }
         }
         start_vec++;
      }
   }

   for ( ubase_t sol = 0; sol < d_sol; sol++ ){
      for ( ubase_t row = 0; row < d_lin; row++ ){
         helper[ row + d_lin * sol ] = matrix[ d_red + row + d_row * ( d_vec - d_sol + sol ) ];
      }
   }

   delete [] matrix;

}

void jans::sieve::__sieving_grace__( const ubase_t blk_size, const double grace ){

   jans::big_int work1; jans::big_int::prod( work1, target, 2 ); // 2N
   jans::big_int work2;
   jans::big_int lower_up; jans::big_int::ceil_sqrt( lower_up, target ); // ceil( sqrt( N ) )
   jans::big_int upper_up;
   jans::big_int limit_up; jans::big_int::ceil_sqrt( limit_up, work1 ); // ceil( sqrt( 2N ) )
   jans::big_int upper_dn; upper_dn.copy( lower_up ); jans::big_int::minus1( upper_dn ); // ceil( sqrt( N ) ) - 1
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
            jans::big_int::plus1( work1 );
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
            jans::big_int::minus1( work1 );
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

