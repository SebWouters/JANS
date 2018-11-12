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

   __startup__();

   extra_sz = extra;
   linspace = num_primes + extra_sz;
   lincount = 0;
   xvalues  = new jans::big_int[ linspace ];
   powers   = new ubase_t[ linspace * num_primes ];
   for ( int cnt = 0; cnt < ( linspace * num_primes ); cnt++ ){ powers[ cnt ] = 0; }

}

jans::sieve::~sieve(){

   delete [] primes;
   delete [] roots;
   delete [] logval;

   delete [] xvalues;
   delete [] powers;

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

void jans::sieve::run( jans::big_int & p, jans::big_int & q, const int blk_size, const double grace ){

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
            if ( helper[ vec + linspace * sol ] == 1 ){ pow += powers[ ip + num_primes * vec ]; }
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

void jans::sieve::__solve_gaussian__( unsigned char * helper ){

   const int d_row = linspace;
   const int d_col = num_primes;

   unsigned char * matrix = new unsigned char[ d_row * ( d_col + d_row ) ];

   for ( int row = 0; row < d_row; row++ ){
      for ( int col = 0; col < ( d_col + d_row ); col++ ){
         matrix[ row + d_row * col ] = 0;
      }
      matrix[ row + d_row * ( d_col + row ) ] = 1;
   }

   for ( int row = 0; row < d_row; row++ ){
      for ( int col = 0; col < d_col; col++ ){
         matrix[ row + d_row * col ] = powers[ col + d_col * row ] % 2;
      }
   }

   int start_row = 0;
   for ( int p = 0; p < d_col; p++ ){
      bool found = false;
      int iter = start_row;
      while ( ( found == false ) && ( iter < d_row ) ){
         if ( matrix[ iter + d_row * p ] == 1 ){ found = true; }
         else { iter++; }
      }
      if ( found == true ){
         if ( iter != start_row ){
            for ( int col = 0; col < ( d_col + d_row ); col++ ){
               matrix[ start_row + d_row * col ] = ( matrix[ iter + d_row * col ] + matrix[ start_row + d_row * col ] ) % 2;
            }
         }
         for ( int row = start_row + 1; row < d_row; row++ ){
            if ( matrix[ row + d_row * p ] == 1 ){
               for ( int col = 0; col < ( d_col + d_row ); col++ ){
                  matrix[ row + d_row * col ] = ( matrix[ row + d_row * col ] + matrix[ start_row + d_row * col ] ) % 2;
               }
            }
         }
         start_row++;
      }
   }

   for ( int sol = 0; sol < ( d_row - d_col ); sol++ ){
      for ( int col = 0; col < d_row; col++ ){
         helper[ col + d_row * sol ] = matrix[ start_row + sol + d_row * ( d_col + col ) ];
      }
   }

   delete [] matrix;

}

void jans::sieve::__sieving_grace__( const int blk_size, const double grace ){

   jans::big_int work1; jans::big_int::prod( work1, target, 2 ); // 2N
   jans::big_int work2;
   jans::big_int lower; jans::big_int::ceil_sqrt( lower, target ); // ceil( sqrt( N ) )
   jans::big_int upper;
   jans::big_int limit; jans::big_int::ceil_sqrt( limit, work1 ); // ceil( sqrt( 2N ) )

   double  * sumlog   = new  double[ blk_size   ];
   ubase_t * offset_p = new ubase_t[ num_primes ];
   ubase_t * offset_n = new ubase_t[ num_primes ];
   ubase_t * helper   = new ubase_t[ num_primes ];

   for ( int ip = 0; ip < num_primes; ip++ ){

      ubase_t pri = primes[ ip ];
      ubase_t rem = jans::big_int::div( work1, lower, pri ); // ceil( sqrt( N ) ) % p
      ubase_t pos =     pri + roots[ ip ] - rem; while( pos >= pri ){ pos -= pri; } offset_p[ ip ] = pos;
      ubase_t neg = 2 * pri - roots[ ip ] - rem; while( neg >= pri ){ neg -= pri; } offset_n[ ip ] = neg;

   }

   jans::big_int::sum( upper, lower, blk_size );
   if ( jans::big_int::smaller( upper, limit ) == false ){ upper.copy( limit ); }

   while ( ( lincount < linspace ) && ( jans::big_int::smaller( lower, limit ) ) ){

      jans::big_int::diff( work1, upper, lower );
      const ubase_t loopsize = work1.get_blk( 0 );
      for ( ubase_t cnt = 0; cnt < loopsize; cnt++ ){ sumlog[ cnt ] = 0.0; }

      for ( int ip = 0; ip < num_primes; ip++ ){
         const ubase_t prime   =   primes[ ip ];
         const double logvalue =   logval[ ip ];
         ubase_t idx           = offset_p[ ip ];
         while ( idx < loopsize ){
            sumlog[ idx ] += logvalue;
            idx += prime;
         }
         offset_p[ ip ] = idx - loopsize;
      }

      for ( int ip = 1; ip < num_primes; ip++ ){ // Not for prime 2
         const ubase_t prime   =   primes[ ip ];
         const double logvalue =   logval[ ip ];
         ubase_t idx           = offset_n[ ip ];
         while ( idx < loopsize ){
            sumlog[ idx ] += logvalue;
            idx += prime;
         }
         offset_n[ ip ] = idx - loopsize;
      }

      int count1 = 0;
      int count2 = 0;
      for ( ubase_t cnt = 0; cnt < loopsize; cnt++ ){
         jans::big_int::sum( work1, lower, cnt );
         jans::big_int::xx_min_num( work2, work1, target );
         const double reference = jans::big_int::logarithm( work2 ) - grace;
         if ( sumlog[ cnt ] > reference ){
            count1++;
            const bool smooth = __extract__( work2, helper );
            if ( smooth ){
               count2++;
               xvalues[ lincount ].copy( work1 );
               std::cout << "Smooth ( x * x - N ) no. " << lincount << " for value x = " << xvalues[ lincount ].write( 10 ) << std::endl;
               for ( int ip = 0; ip < num_primes; ip++ ){
                  powers[ lincount * num_primes + ip ] = helper[ ip ];
               }
               lincount++;
               if ( lincount == linspace ){ cnt = loopsize; }
            }
         }
      }
      std::cout << "In the present block " << count1 << " / " << loopsize << " were selected based on sums of logarithms." << std::endl;
      std::cout << "In the present block " << count2 << " / " << count1   << " of the candidates were smooth." << std::endl;

      lower.copy( upper );
      jans::big_int::sum( upper, lower, blk_size );
      if ( jans::big_int::smaller( upper, limit ) == false ){ upper.copy( limit ); }

   }

   delete [] sumlog;
   delete [] offset_p;
   delete [] offset_n;
   delete [] helper;

}

int jans::sieve::__legendre_symbol__( jans::big_int & num, const ubase_t p ){

   jans::big_int work;
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
   return ( ( ubase_t )( res ) );

}

ubase_t jans::sieve::__root_quadratic_residue__( jans::big_int & num, const ubase_t p ){

   jans::big_int work;
   const ubase_t rem = jans::big_int::div( work, num, p ); // rem = num % p
   return __root_quadratic_residue__( rem, p );

}

ubase_t jans::sieve::__root_quadratic_residue__( const ubase_t num, const ubase_t p ){

   // Tonelli–Shanks

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

   ubase_t M = S;
   ubase_t c = __power__(   z, Q, p );             // c^{2^{M-1}} = z^{(p-1)/2} = -1 mod p
   ubase_t t = __power__( rem, Q, p );             // t^{2^{M-1}} = n^{(p-1)/2} = +1 mod p
   ubase_t R = __power__( rem, ( Q + 1 ) / 2, p ); // R^2 = n.t mod p

   while ( t != 1 ){

      /*
          Find min. 0 < i < M so that t^{2^i} mod p = 1
          Guaranteed to exist: - while loop condition
                               - t^{2^{M-1}} = +1 mod p
      */
      ubase_t i = 0;
      ucarry_t b = t;
      while ( ( b != 1 ) && ( i < M ) ){
         i++;
         b = ( b * b ) % p; // t^{2^i} mod p
      }

      // Calculate c^{2^{M-i-1}} mod p
      ubase_t j = 0;
      b = c;
      while ( j < M - i - 1 ){
         j++;
         b = ( b * b ) % p; // c^{2^j} mod p
      }

      M = i;                             // M decreases each iteration!
      c = ( b * b ) % p;                 // c'^{2^{M'-1}} =   (b^2)^{2^{i-1}} = c^{2^{M-1}}             = -1   mod p
      t = ( ( ( t * b ) % p ) * b ) % p; // t'^{2^{M'-1}} = (t.b^2)^{2^{i-1}} = t^{2^{i-1}}.c^{2^{M-1}} = +1   mod p
      R = ( R * b ) % p;                 // R'^2          = (R.b)^2           = n.t.b^2                 = n.t' mod p

   }

   return R;

}

bool jans::sieve::__extract__( big_int & x, ubase_t * helper ) const{

   for ( int ip = 0; ip < num_primes; ip++ ){
      helper[ ip ] = jans::big_int::extract_pow_p( x, primes[ ip ] );
      if ( jans::big_int::equal( x, 1 ) ){
         for ( int it = ip + 1; it < num_primes; it++ ){ helper[ it ] = 0; }
         return true;
      }
   }

   return false;

}

