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

#ifndef jans_big_int
#define jans_big_int

#include <string>
#include <limits.h>

#define BLOCK_BIT ( sizeof(unsigned char) * CHAR_BIT )
#define NUM_BLOCK 128

#define __00000000__ (  ((unsigned char)(0)) )
#define __11111111__ ( ~((unsigned char)(0)) )

#define SET_BIT(A,k) ( A[ ( ( k ) / ( sizeof(unsigned char) * CHAR_BIT ) ) ] |=  ( 1 << ( ( k ) % ( sizeof(unsigned char) * CHAR_BIT ) ) ) )
#define DEL_BIT(A,k) ( A[ ( ( k ) / ( sizeof(unsigned char) * CHAR_BIT ) ) ] &= ~( 1 << ( ( k ) % ( sizeof(unsigned char) * CHAR_BIT ) ) ) )
#define GET_BIT(A,k) ( A[ ( ( k ) / ( sizeof(unsigned char) * CHAR_BIT ) ) ] &   ( 1 << ( ( k ) % ( sizeof(unsigned char) * CHAR_BIT ) ) ) )

namespace jans{

   class big_int{

      public:

         big_int();

         virtual ~big_int();

         void set( const std::string number, const unsigned char base );

         static bool compare( big_int & n1, big_int & n2 );

         static int n_bits(){ return ( NUM_BLOCK * BLOCK_BIT ); }

         static int n_blocks(){ return NUM_BLOCK; }

         static void sanity_check();

         std::string str( const unsigned char base );

         std::string str_hex();

      private:

         unsigned char data[ NUM_BLOCK ];

         bool sign;

         static void __clear__( unsigned char * a );

         static void __add3kernel__( unsigned char * res, unsigned char * a, unsigned char * b, const int is );

         static void __multiply3kernel__( unsigned char * res, unsigned char * temp, unsigned char * a, unsigned char * b );

         static void __multiply3kernel__( unsigned char * res, unsigned char * a, const unsigned char b, const int ib );

         static const char __conversion__[ 16 ];

         static unsigned char __convert_c2i__( const char c );

         static char __convert_i2c__( const unsigned char c );

         static bool __compare__( unsigned char * a, unsigned char * b );

   };

}

#endif

