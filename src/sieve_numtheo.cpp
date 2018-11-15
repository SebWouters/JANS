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

#include "sieve.h"

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

int jans::sieve::__legendre_symbol__( jans::big_int & num, jans::big_int & p ){

   // Algorithm 2.3.5, Crandall & Pomerance

   assert( jans::big_int::even( p ) == false );
   jans::big_int w;
   jans::big_int q;
   jans::big_int a; jans::big_int::div( w, a, num, p ); // a = num % p
   jans::big_int m; m.copy( p ); // m = p

   int t = 1;
   ubase_t rem;

   while ( jans::big_int::equal( a, 0 ) == false ){
      while ( jans::big_int::even( a ) ){
         rem = jans::big_int::div( w, a, 2 );
         a.copy( w );
         assert( rem == 0 );
         rem = jans::big_int::div( w, m, 8 );
         if ( ( rem == 3 ) || ( rem == 5 ) ){ t = -t; }
      }
      if ( ( jans::big_int::div( w, a, 4 ) == 3 ) && ( jans::big_int::div( q, m, 4 ) == 3 ) ){ t = -t; }
      w.copy( m ); // w = m(old) = a(new)
      m.copy( a );
      a.copy( w );
      jans::big_int::div( q, a, w, m ); // a = a % m
   }
   if ( jans::big_int::equal( m, 1 ) ){ return t; }
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

ubase_t jans::sieve::__inv_x_mod_p__( jans::big_int & x, const ubase_t p ){

   jans::big_int quot;
   const ubase_t rem = jans::big_int::div( quot, x, p ); // rem = x % p
   assert( jans::big_int::equal( quot, 0 ) == false );

   int     a = 0;   // a = u_ini = 0
   ubase_t g = p;   // g = w_ini = p
   int     u = 1;   // u = a_ini - q * u_ini = 1 - q * 0 = 1
   ubase_t w = rem; // w = g_ini - q * w_ini = x - q * p = x % p

   while ( w > 0 ){
      ubase_t q = g / w;
      ubase_t r = g - q * w;
      ubase_t s = a - q * u;
      a = u;
      g = w;
      u = s;
      w = r;
   }
   assert( g == 1 );

   while ( a < 0 ){ a = a + p; }
   ubase_t inverse_x = a;
   assert( 1 == ( ( ( ucarry_t )( rem ) ) * inverse_x ) % p );
   return inverse_x;

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
      z = ( ( jans::big_int::random_ubase_t() ) % ( p - 2 ) ) + 2; // random number in [ 2 .. p - 1 ]
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
      while ( j + i + 1 < M ){
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

