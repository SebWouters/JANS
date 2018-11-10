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

void jans::big_int::prod( big_int & res, big_int & a, const ubase_t b ){

   res.sign = a.sign;
   res.lead = __mult2set__( res.data, a.data, a.lead, b, 0 );

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
         lq = __sum1__( q.data, lq, 1 );
      }
      q.sign = ( ( lq > 0 ) ? false : true );

   }

   q.lead = lq;
   r.lead = lr;

}

ubase_t jans::big_int::div( big_int & q, big_int & n, const ubase_t d ){

   // 0 <= r < d

   int lq = n.lead;
   __copy__( q.data, n.data );
   ubase_t rem = __divide__( q.data, lq, d );
   q.sign = true;

   if ( n.sign == false ){

      // - |n| = - d * ( return( q ) + 1 ) + ( d - return( r ) )
      if ( rem > 0 ){
         rem = d - rem;
         lq = __sum1__( q.data, lq, 1 );
      }
      q.sign = ( ( lq > 0 ) ? false : true );

   }

   q.lead = lq;
   return rem;

}

void jans::big_int::gcd( big_int & res, big_int & a, big_int & b ){

   assert( a.sign );
   assert( b.sign );
   res.sign = true;

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

   assert( a.sign );
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





