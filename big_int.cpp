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
 * addition itf & test
 * bin to str(dec)
 * modulo
 * division
 * subtraction
 * gcd (Euclidean algorithm)
 * modular addition / subtraction
 * modular multiplication (Kochanski, Montgomery)

*/

const char jans::big_int::__conversion__[ 16 ] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };

jans::big_int::big_int(){

   sign = true;
   __clear__( data );

}

jans::big_int::~big_int(){}

void jans::big_int::sanity_check(){

   assert( sizeof(unsigned int) > sizeof(unsigned char) );

}

void jans::big_int::__clear__( unsigned char * a ){

   for ( int i = 0; i < NUM_BLOCK; i++ ){ a[ i ] &= __00000000__; }

}

void jans::big_int::__add3kernel__( unsigned char * res, unsigned char * a, unsigned char * b, const int is ){

   unsigned int z = 0;
   for ( int i = is; i < NUM_BLOCK; i++ ){
      unsigned int x = a[ i ];
      unsigned int y = b[ i ];
      z = x + y + z;
      res[ i ] = z & __11111111__;
      z = z >> BLOCK_BIT;
   }
   assert( z == 0 ); // Else overflow exception

}

void jans::big_int::__multiply3kernel__( unsigned char * res, unsigned char * temp, unsigned char * a, unsigned char * b ){

   // Block version of Russian peasant
   for ( int ib = 0; ib < NUM_BLOCK; ib++ ){
      if ( b[ ib ] != __00000000__ ){
         jans::big_int::__multiply3kernel__( temp, a, b[ ib ], ib );
         jans::big_int::__add3kernel__( res, temp, res, ib );
      }
   }

}

void jans::big_int::__multiply3kernel__( unsigned char * res, unsigned char * a, const unsigned char b, const int ib ){

   unsigned int z = 0;
   unsigned int y = b;
   for ( int ia = 0; ia < NUM_BLOCK - ib; ia++ ){
      unsigned int x = a[ ia ];
      z = z + ( x * y );
      res[ ia + ib ] = z & __11111111__;
      z = z >> BLOCK_BIT;
   }
   for ( int ia = NUM_BLOCK - ib; ia < NUM_BLOCK; ia++ ){
      assert( a[ ia ] != __00000000__ ); // Else overflow exception
   }
   assert( z == 0 ); // Else overflow exception

}

bool jans::big_int::compare( big_int & n1, big_int & n2 ){

   return ( ( n1.sign == n2.sign ) & __compare__( n1.data, n2.data ) );

}

bool jans::big_int::__compare__( unsigned char * a, unsigned char * b ){

   for ( int i = 0; i < NUM_BLOCK; i++ ){
      if ( a[ i ] != b[ i ] ){ return false; }
   }
   return true;

}

void jans::big_int::set( const std::string number, const unsigned char base ){

   assert( ( base == 2 ) || ( base == 10 ) || ( base == 16 ) ); // It can be done more efficient for base 2 and 16, but it's for free.
   __clear__( data );

   unsigned char shift[ NUM_BLOCK ]; __clear__( shift );
   unsigned char temp [ NUM_BLOCK ]; __clear__( temp  );
   shift[ 0 ] = 1;

   for ( int i = number.size() - 1; i >= 0; i-- ){
      const unsigned char digit = __convert_c2i__( number.at( i ) );
      assert( digit != 17 );
      if ( digit != __00000000__ ){
         __multiply3kernel__( temp, shift, digit, 0 );
         __add3kernel__( data, temp, data, 0 );
      }
      __multiply3kernel__( shift, shift, base, 0 );
   }

}

unsigned char jans::big_int::__convert_c2i__( const char c ){

   for ( unsigned char i = 0; i < 16; i++ ){ if ( c == __conversion__[ i ] ){ return i; } }
   return 17; // Error code

}

char jans::big_int::__convert_i2c__( const unsigned char c ){

   if ( c >= 16 ){ return 'z'; } // Error code
   return __conversion__[ c ];

}

std::string jans::big_int::str( const unsigned char base ){

   assert( ( base == 2 ) || ( base == 10 ) || ( base == 16 ) );

   if ( ( base == 2 ) || ( base == 16 ) ){

      const int log2base  = ( ( base == 2 ) ? 1 : 4 );
      const int text_size = n_bits() / log2base;
      char text[ text_size ];

      for ( int i = 0; i < NUM_BLOCK; i++ ){
         for ( int j = 0; j < ( BLOCK_BIT / log2base ); j++ ){
            text[ text_size - 1 - ( ( BLOCK_BIT / log2base ) * i + j ) ]
               = __conversion__[ ( ( data[ i ] & ( ( base - 1 ) << ( log2base * j ) ) ) >> ( log2base * j ) ) ];
         }
      }

      std::string reduction( text, 0, text_size ); // Clear clutter at end
      const int start  = reduction.find_first_not_of( '0' );
      const int length = reduction.size() - start;
      reduction = reduction.substr( start, length );
      return reduction;

   }

   return "error";

}





