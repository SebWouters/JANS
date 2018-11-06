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

#include "big_int.h"

/*

TODO:

 * sieve : gcd (Euclidean algorithm): gcd(a,b) --> if b divides a, then b, else gcd(b,a mod b)
 * sieve : sqrt( N ) : perhaps the Babylonian method: x_new = ( x + N / x ) / 2?
 * sieve : multiplication mod N; ( a * b ) mod N = ( a mod N * b mod N ) mod N, or is there something fancier?
 * sieve : With N = 2^32, N/ln(N) yields 193635250 primes. I guess base_t division will suffice for subset of primes?
 * sieve : Perhaps make an int jans::big_int::extract( big_int & num, int & ln, base_t factor ); x = extract( n_in, f ); --> n_out * f^x = n_in and ( n_out mod f ) == 0
 *
 * modular addition / subtraction
 * modular multiplication (Kochanski, Montgomery)

*/

jans::big_int::big_int(){

   sign = true;
   lead = 0;
   __clear__( data );

}

jans::big_int::~big_int(){}

void jans::big_int::sanity_check(){

   std::cout << "num_bits(ubase_t)  = " << sizeof(ubase_t)  * CHAR_BIT << std::endl;
   std::cout << "num_bits(ucarry_t) = " << sizeof(ucarry_t) * CHAR_BIT << std::endl;
   assert( sizeof(ucarry_t) > sizeof(ubase_t) );

}

bool jans::big_int::equal( big_int & n1, big_int & n2 ){

   return ( ( n1.sign == n2.sign ) && ( n1.lead == n2.lead ) && ( __compare__( n1.data, n2.data ) == 0 ) );

}

void jans::big_int::diff( big_int & res, big_int & a, big_int & b ){

   const int comp = __compare__( a.data, b.data );

   if ( ( a.sign == true  ) && ( b.sign == true  ) ){
      if ( comp == 0 ){ res.sign = true;  res.lead = 0; __clear__( res.data ); }
      if ( comp >  0 ){ res.sign = true;  res.lead = __diff3set__( res.data, a.data, b.data ); } // a - b =   ( |a| - |b| )
      if ( comp <  0 ){ res.sign = false; res.lead = __diff3set__( res.data, b.data, a.data ); } // a - b = - ( |b| - |a| )
   }
   if ( ( a.sign == false ) && ( b.sign == false ) ){
      if ( comp == 0 ){ res.sign = true;  res.lead = 0; __clear__( res.data ); }
      if ( comp >  0 ){ res.sign = false; res.lead = __diff3set__( res.data, a.data, b.data ); } // a - b = - ( |a| - |b| )
      if ( comp <  0 ){ res.sign = true;  res.lead = __diff3set__( res.data, b.data, a.data ); } // a - b =   ( |b| - |a| )
   }
   if ( ( a.sign == true  ) && ( b.sign == false ) ){ res.sign = true;  res.lead = __sum3set__( res.data, a.data, a.lead, b.data, b.lead ); }
   if ( ( a.sign == false ) && ( b.sign == true  ) ){ res.sign = false; res.lead = __sum3set__( res.data, a.data, a.lead, b.data, b.lead ); }

}

void jans::big_int::sum( big_int & res, big_int & a, big_int & b ){

   const int comp = __compare__( a.data, b.data );

   if ( ( a.sign == true  ) && ( b.sign == false ) ){
      if ( comp == 0 ){ res.sign = true;  res.lead = 0; __clear__( res.data ); }
      if ( comp >  0 ){ res.sign = true;  res.lead = __diff3set__( res.data, a.data, b.data ); } // a + b =   ( |a| - |b| )
      if ( comp <  0 ){ res.sign = false; res.lead = __diff3set__( res.data, b.data, a.data ); } // a + b = - ( |b| - |a| )
   }
   if ( ( a.sign == false ) && ( b.sign == true  ) ){
      if ( comp == 0 ){ res.sign = true;  res.lead = 0; __clear__( res.data ); }
      if ( comp >  0 ){ res.sign = false; res.lead = __diff3set__( res.data, a.data, b.data ); } // a + b = - ( |a| - |b| )
      if ( comp <  0 ){ res.sign = true;  res.lead = __diff3set__( res.data, b.data, a.data ); } // a + b =   ( |b| - |a| )
   }
   if ( ( a.sign == true  ) && ( b.sign == true  ) ){ res.sign = true;  res.lead = __sum3set__( res.data, a.data, a.lead, b.data, b.lead ); }
   if ( ( a.sign == false ) && ( b.sign == false ) ){ res.sign = false; res.lead = __sum3set__( res.data, a.data, a.lead, b.data, b.lead ); }

}

void jans::big_int::prod( big_int & res, big_int & a, big_int & b ){

   res.sign = !( a.sign != b.sign );
   res.lead = __mult3set__( res.data, a.data, a.lead, b.data, b.lead );

}

void jans::big_int::div( big_int & q, big_int & r, big_int & n, big_int & d ){

   // 0 <= r < d

   assert( d.sign == true );

   int lq = 0;
   int lr = n.lead;
   __copy__( r.data, n.data );
   __divide__( q.data, lq, r.data, lr, d.data, d.lead );
   q.sign = true;
   r.sign = true;

   if ( n.sign == false ){

      // - |n| = - d * ( return( q ) + 1 ) + ( d - return( r ) )
      if ( lr > 0 ){
         lr = __diff3set__( r.data, d.data, r.data );
         lq = __plus_one__( q.data, lq );
      }
      q.sign = ( ( lq > 0 ) ? false : true );

   }

   q.lead = lq;
   r.lead = lr;

}

void jans::big_int::shift_up( const int k ){

   __shift_up__( data, k );
   lead = 0;
   for ( int i = 0; i < NUM_BLOCK; i++ ){ if ( data[ i ] != 0 ){ lead = i + 1; } }

}

void jans::big_int::shift_down( const int k ){

   __shift_down__( data, k );
   lead = 0;
   for ( int i = 0; i < NUM_BLOCK; i++ ){ if ( data[ i ] != 0 ){ lead = i + 1; } }

}





