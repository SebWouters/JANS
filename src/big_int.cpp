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

#include "big_int.h"

int jans::big_int::NUM_BLOCK = 0;

bool jans::big_int::nb_set = false;

jans::big_int::big_int(){

   assert( nb_set );

   lead = 0;
   data = new ubase_t[ NUM_BLOCK ];
   __clear__( data );

}

jans::big_int::~big_int(){

   delete [] data;

}

void jans::big_int::copy( big_int & tocopy ){

   lead = tocopy.lead;
   __copy__( data, tocopy.data );

}

void jans::big_int::copy( const ubase_t value ){

   __clear__( data );
   data[ 0 ] = value;
   lead = ( ( value == 0 ) ? 0 : 1 );

}

void jans::big_int::sanity_check(){

   assert( RAND_MAX >= ( 1UL << CHAR_BIT ) );
   assert( sizeof( ucarry_t ) >= 2 * sizeof( ubase_t ) );

}

void jans::big_int::set_num_block( const int factor ){

   assert( nb_set == false );
   assert( factor > 0 );
   NUM_BLOCK = ( ( factor * BASE_UNIT ) / BLOCK_BIT );
   nb_set = true;

}

bool jans::big_int::equal( big_int & n1, big_int & n2 ){

   return ( ( n1.lead == n2.lead ) && ( __compare__( n1.data, n2.data ) == 0 ) );

}

bool jans::big_int::equal( big_int & n1, const ubase_t n2 ){

   if ( n2 == 0 ){ return ( n1.lead == 0 ); }
   return ( ( n1.lead == 1 ) && ( n1.data[ 0 ] == n2 ) );

}

bool jans::big_int::smaller( big_int & n1, big_int & n2 ){

   return ( __compare__( n1.data, n2.data ) < 0 );

}

bool jans::big_int::smaller( big_int & n1, const ubase_t n2 ){

   if ( n1.lead > 1 ){ return false; }
   return ( n1.data[ 0 ] < n2 );

}

bool jans::big_int::even( big_int & n1 ){

   return ( ( n1.data[ 0 ] % 2 ) == 0 );

}

void jans::big_int::plus( big_int & res, const ubase_t val ){

   res.lead = __sum1__( res.data, res.lead, val );

}

void jans::big_int::sum( big_int & res, big_int & a, big_int & b ){

   res.lead = __sum3set__( res.data, a.data, a.lead, b.data, b.lead );

}

void jans::big_int::sum( big_int & res, big_int & a, const ubase_t b ){

   res.lead = a.lead;
   __copy__( res.data, a.data );
   res.lead = __sum1__( res.data, res.lead, b );

}

void jans::big_int::minus( big_int & res, const ubase_t val ){

   res.lead = __diff1__( res.data, val );

}

void jans::big_int::diff( big_int & res, big_int & a, big_int & b ){

   res.lead = __diff3set__( res.data, a.data, b.data );

}

void jans::big_int::diff( big_int & res, big_int & a, const ubase_t b ){

   __copy__( res.data, a.data );
   res.lead = __diff1__( res.data, b );

}

void jans::big_int::prod( big_int & res, big_int & a, big_int & b ){

   res.lead = __mult3set__( res.data, a.data, a.lead, b.data, b.lead );

}

void jans::big_int::prod( big_int & res, big_int & a, const ubase_t b ){

   res.lead = __mult2set__( res.data, a.data, a.lead, b, 0 );

}

void jans::big_int::div( big_int & q, big_int & r, big_int & n, big_int & d ){

   r.lead = n.lead;
   __copy__( r.data, n.data );
   __divide__( q.data, q.lead, r.data, r.lead, d.data, d.lead );

}

ubase_t jans::big_int::div( big_int & q, big_int & n, const ubase_t d ){

   q.lead = n.lead;
   __copy__( q.data, n.data );
   ubase_t rem = __divide__( q.data, q.lead, d );
   return rem;

}

void jans::big_int::gcd( big_int & res, big_int & a, big_int & b ){

   const int comp = __compare__( a.data, b.data );
   if ( comp == 0 ){ // a == b
      __copy__( res.data, a.data );
      res.lead = a.lead;
      return;
   }

   ubase_t a_cpy[ NUM_BLOCK ]; __copy__( a_cpy, a.data );
   ubase_t b_cpy[ NUM_BLOCK ]; __copy__( b_cpy, b.data );

   if ( comp > 0 ){ res.lead = __gcd__( res.data, a_cpy, a.lead, b_cpy, b.lead ); } // a > b
             else { res.lead = __gcd__( res.data, b_cpy, b.lead, a_cpy, a.lead ); } // b > a

}

void jans::big_int::prodmod( big_int & q, big_int & r, big_int & a, big_int & b, big_int & m ){

   r.lead = __mult3set__( r.data, a.data, a.lead, b.data, b.lead );
   __divide__( q.data, q.lead, r.data, r.lead, m.data, m.lead );

}

void jans::big_int::power( big_int & res, big_int & base, big_int & expo, big_int & mod ){

   big_int work1;
   big_int work2;
   big_int work3;
   big_int junk;

   div( junk, work1, expo, mod ); // work1 = expo % mod
   div( junk, work2, base, mod ); // work2 = base % mod

   res.copy( 1 );

   for ( int ie = 0; ie < work1.lead; ie++ ){
      for ( int je = 0; je < BLOCK_BIT; je++ ){
         const bool to_multiply = ( ( work1.data[ ie ] >> je ) & 1U );
         if ( to_multiply ){
            prod( work3, work2, res );
            div( junk, res, work3, mod ); // res = ( work2 * res ) % mod
         }
         prod( work3, work2, work2 );
         div( junk, work2, work3, mod ); // work2 = ( base ) ^ ( 2 ^ ( BLOCK_BIT * ie + je + 1 ) )
      }
   }

}

void jans::big_int::ceil_sqrt( big_int & res, big_int & n ){

   long double number = i2f( n );
   number = sqrt( number );
   f2i( res, number );
   res.lead = __ceil_sqrt__( res.data, res.lead, n.data, n.lead );

}

long double jans::big_int::i2f( big_int & x ){

   long double result = 0.0;
   long double base   = 1.0;

   for ( int i = 0; i < x.lead; i++ ){
      result = result + ( base * x.data[ i ] );
      base   = base * ( 1UL << BLOCK_BIT );
   }

   return result;

}

void jans::big_int::f2i( big_int & x, const long double number ){

   long double remainder = number;
   long double index     = log2( remainder );

   __clear__( x.data );
   x.lead = 0;

   while ( index >= 0.0 ){

      const int idx = ( int )( index );
      const int ix  = idx / BLOCK_BIT;
      const int jx  = idx % BLOCK_BIT;

      x.data[ ix ] |= ( 1U << jx );
      if ( ix + 1 > x.lead ){ x.lead = ix + 1; }

      long double set = 1.0;
      for ( int cnt = 0; cnt < ix; cnt++ ){ set = set * ( 1UL << BLOCK_BIT ); }
      set = set * ( 1U << jx );

      remainder = remainder - set;
      index = log2( remainder );

   }

}

ubase_t jans::big_int::extract_pow_p( big_int & x, const ubase_t p ){

   assert( p > 1 );

   ubase_t pow = 0;
   ubase_t rem = 0;

   ubase_t work[ NUM_BLOCK ];
   __copy__( work, x.data );
   int lw = x.lead;

   while ( ( rem == 0 ) && ( lw > 0 ) ){
      rem = __divide__( work, lw, p );
      if ( rem == 0 ){
         pow++;
         x.lead = lw;
         __copy__( x.data, work );
      }
   }

   return pow;

}

ubase_t jans::big_int::random_ubase_t(){

   ubase_t number = 0;
   for ( int j = 0; j < ( BLOCK_BIT / CHAR_BIT ); j++ ){
      number |= ( ( rand() & ( ( 1U << CHAR_BIT ) - 1 ) ) << ( j * CHAR_BIT ) );
   }
   return number;

}

void jans::big_int::randomize( big_int & n, const int ln ){

   n.lead = ln;
   __clear__( n.data );

   for ( int in = 0; in < ln; in++ ){
      n.data[ in ] = random_ubase_t();
   }

}

bool jans::big_int::miller_rabin( big_int & n, const ubase_t attempts ){

   if ( equal( n, 2 ) ){ return true; } // prime
   if ( equal( n, 3 ) ){ return true; } // prime
   if ( even( n ) ){ return false; } // composite

   big_int u;
   big_int work;
   big_int temp;
   big_int junk;
   big_int check;

   // n - 1 = 2^r * u
   u.copy( n );
   jans::big_int::minus( u, 1 );
   const ubase_t r = extract_pow_p( u, 2 );

   // check = n - 1
   check.copy( n );
   jans::big_int::minus( check, 1 );

   for ( ubase_t cnt = 0; cnt < attempts; cnt++ ){

      // random temp in [ 2 ... n - 2 ]
      do {
         randomize( temp, check.lead );
         __divide__( junk.data, junk.lead, temp.data, temp.lead, check.data, check.lead );
      } while ( ( equal( temp, 0 ) ) || ( equal( temp, 1 ) ) );

      // work = temp ^ u % n
      power( work, temp, u, n );

      bool ctu = ( ( equal( work, 1 ) == false ) && ( equal( work, check ) == false ) );

      for ( ubase_t i = 0; ( ( i < r - 1 ) && ( ctu ) ); i++ ){
         prod( temp, work, work );
         div( junk, work, temp, n );
         ctu = ( equal( work, check ) == false );
      }

      if ( ctu ){ return false; } // composite

   }

   return true;

}

