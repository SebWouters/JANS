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

   assert( sizeof( ucarry_t ) >= 2 * sizeof( ubase_t ) );

}

void jans::big_int::set_num_block( const int factor ){

   assert( nb_set == false );
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

ubase_t jans::big_int::get_blk( const int i ){

   return data[ i ];

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

ubase_t jans::big_int::gcd( big_int & a, const ubase_t b ){

   assert( ( a.lead > 1 ) || ( a.data[ 0 ] >= b ) );

   ubase_t cpy[ NUM_BLOCK ];
   __copy__( cpy, a.data );
   int lq = a.lead;

   ubase_t new_a = b;
   ubase_t new_b = __divide__( cpy, lq, b ); // a mod b

   while( new_b != 0 ){
      ubase_t swap = new_b;
      new_b = new_a % new_b;
      new_a = swap;
   }

   return new_a;

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
         const bool to_multiply = ( ( work1.get_blk( ie ) >> je ) & 1U );
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

   res.lead = __ceil_sqrt__( res.data, n.data, n.lead );

}

double jans::big_int::logarithm( big_int & x ){

   long double result = 0.0;
   long double base   = 1.0;

   for ( int i = 0; i < x.lead; i++ ){
      result = result + ( base * x.data[ i ] );
      base   = base * ( 1UL << BLOCK_BIT );
   }

   result = log( result );
   return ( ( double )( result ) );

}

void jans::big_int::xx_min_num( big_int & res, big_int & x, big_int & num ){

   res.lead = __mult3set__( res.data, x.data, x.lead, x.data, x.lead );
   res.lead = __diff3set__( res.data, res.data, num.data );

}

void jans::big_int::min_xx_plus_num( big_int & res, big_int & x, big_int & num ){

   res.lead = __mult3set__( res.data, x.data, x.lead, x.data, x.lead );
   res.lead = __diff3set__( res.data, num.data, res.data );

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

