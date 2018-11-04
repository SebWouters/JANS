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

 * multiplication itf & test (e.g. prime factors RSA-100)

 * bin to str(dec)
 * modulo
 * division
 * gcd (Euclidean algorithm)
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

void jans::big_int::__clear__( ubase_t * a ){

   for ( int i = 0; i < NUM_BLOCK; i++ ){ a[ i ] = 0; }

}

bool jans::big_int::equal( big_int & n1, big_int & n2 ){

   return ( ( n1.sign == n2.sign ) && ( n1.lead == n2.lead ) && ( __compare__( n1.data, n2.data ) == 0 ) );

}

int jans::big_int::__compare__( ubase_t * a, ubase_t * b ){

   for ( int i = NUM_BLOCK - 1; i >= 0; i-- ){
      if ( a[ i ] > b[ i ] ){ return (  i + 1 ); }
      if ( a[ i ] < b[ i ] ){ return ( -i - 1 ); }
   }
   return 0;

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
   if ( ( a.sign == true  ) && ( b.sign == false ) ){ res.sign = true;  res.lead = __sum3set__( res.data, a.data, a.lead, b.data, b.lead, 0 ); }
   if ( ( a.sign == false ) && ( b.sign == true  ) ){ res.sign = false; res.lead = __sum3set__( res.data, a.data, a.lead, b.data, b.lead, 0 ); }

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
      if ( comp >  0 ){ res.sign = false; res.lead = __diff3set__( res.data, a.data, b.data ); } // a - b = - ( |a| - |b| )
      if ( comp <  0 ){ res.sign = true;  res.lead = __diff3set__( res.data, b.data, a.data ); } // a - b =   ( |b| - |a| )
   }
   if ( ( a.sign == true  ) && ( b.sign == true  ) ){ res.sign = true;  res.lead = __sum3set__( res.data, a.data, a.lead, b.data, b.lead, 0 ); }
   if ( ( a.sign == false ) && ( b.sign == false ) ){ res.sign = false; res.lead = __sum3set__( res.data, a.data, a.lead, b.data, b.lead, 0 ); }

}

