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

#include "big_int.h"

const char jans::big_int::__conversion__[ 16 ] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f' };

ubase_t jans::big_int::__convert_c2i__( const char c ){

   for ( ubase_t i = 0; i < 16; i++ ){ if ( c == __conversion__[ i ] ){ return i; } }
   return 17; // Error code

}

char jans::big_int::__convert_i2c__( const ubase_t c ){

   if ( c >= 16 ){ return 'z'; } // Error code
   return __conversion__[ c ];

}

void jans::big_int::read( const std::string number, const ubase_t base ){

   assert( ( base == 2 ) || ( base == 10 ) || ( base == 16 ) );
   __clear__( data );
   lead = 0;

   ubase_t shift[ NUM_BLOCK ];
   __clear__( shift );
   shift[ 0 ] = 1;
   int ls = 1;

   for ( int i = number.size() - 1; i >= 0; i-- ){
      const ubase_t digit = __convert_c2i__( number.at( i ) );
      assert( digit != 17 );
      if ( digit != 0 ){
         lead = __mult2add__( data, lead, shift, ls, digit, 0 );
      }
      ls = __scal1__( shift, ls, base );
   }

}

std::string jans::big_int::write( const ubase_t base ) const{

   assert( ( base == 2 ) || ( base == 10 ) || ( base == 16 ) );

   if ( lead == 0 ){ return "0"; }

   if ( ( base == 2 ) || ( base == 16 ) ){

      const int log2base  = ( ( base == 2 ) ? 1 : 4 );
      const int text_size = ( lead * BLOCK_BIT ) / log2base;
      char text[ text_size ];

      for ( int i = 0; i < lead; i++ ){
         for ( int j = 0; j < ( BLOCK_BIT / log2base ); j++ ){
            text[ text_size - 1 - ( ( BLOCK_BIT / log2base ) * i + j ) ]
               = __conversion__[ ( ( data[ i ] & ( ( base - 1 ) << ( log2base * j ) ) ) >> ( log2base * j ) ) ];
         }
      }

      std::string reduction;
      reduction.append( text, text_size ); // Clear clutter at end
      const int start  = reduction.find_first_not_of( '0' );
      const int length = reduction.size() - start;
      reduction = reduction.substr( start, length );
      return reduction;

   }

   const double log2base = 3.32192809489;
   const int text_size = ( ( int )( ( lead * BLOCK_BIT ) / log2base ) ) + 4;
   char text[ text_size ];
   for ( int c = 0; c < text_size; c++ ){ text[ c ] = '0'; }

   ubase_t num[ NUM_BLOCK ];
   __copy__( num, data );
   int ln = lead;
   int c = 0;
   while ( ( ln > 0 ) || ( num[ 0 ] != 0 ) ){
      text[ text_size - 1 - c ] = __convert_i2c__( __divide__( num, ln, base ) );
      c++;
   }

   std::string reduction;
   reduction.append( text, text_size ); // Clear clutter at end
   const int start  = reduction.find_first_not_of( '0' );
   const int length = reduction.size() - start;
   reduction = reduction.substr( start, length );
   return reduction;

}

