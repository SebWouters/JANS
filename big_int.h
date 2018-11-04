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

#ifndef JANS_BIG_INT
#define JANS_BIG_INT

#include <string>
#include <limits.h>

#define ubase_t  unsigned int
#define ucarry_t unsigned long long

#define BLOCK_BIT ( sizeof(unsigned int) * CHAR_BIT )
#define NUM_BLOCK ( 1024 / BLOCK_BIT )

#define __00000000__ (  ((unsigned int)(0)) )
#define __11111111__ ( ~((unsigned int)(0)) )

#define SET_BIT(A,k) ( A[ ( ( k ) / ( sizeof(unsigned int) * CHAR_BIT ) ) ] |=  ( 1 << ( ( k ) % ( sizeof(unsigned int) * CHAR_BIT ) ) ) )
#define DEL_BIT(A,k) ( A[ ( ( k ) / ( sizeof(unsigned int) * CHAR_BIT ) ) ] &= ~( 1 << ( ( k ) % ( sizeof(unsigned int) * CHAR_BIT ) ) ) )
#define GET_BIT(A,k) ( A[ ( ( k ) / ( sizeof(unsigned int) * CHAR_BIT ) ) ] &   ( 1 << ( ( k ) % ( sizeof(unsigned int) * CHAR_BIT ) ) ) )

namespace jans{

   class big_int{

      public:

         big_int();

         virtual ~big_int();

         static bool compare( big_int & n1, big_int & n2 );

         static int n_bits(){ return ( NUM_BLOCK * BLOCK_BIT ); }

         static int n_blocks(){ return NUM_BLOCK; }

         static void sanity_check();

         void read( const std::string number, const ubase_t base );

         std::string write( const ubase_t base );

      private:

         ubase_t data[ NUM_BLOCK ];

         bool sign;

         int lead; // Upper bound for loops over the blocks: lead = 1 + max{i}( data[ i ] != 0 )

         static void __clear__( ubase_t * a );

         static bool __compare__( ubase_t * a, ubase_t * b );

         /********
          *  IO  *
          ********/

         static const char __conversion__[ 16 ];

         static ubase_t __convert_c2i__( const char c );

         static char __convert_i2c__( const ubase_t c );

         /*************************
          *  Basic math routines  *
          *************************/

         // r[ start : ] = a[ start : ] + b[ start : ]
         static int __add3kernel__( ubase_t * r, ubase_t * a, const int la, ubase_t * b, const int lb, const int start = 0 );

         // r[ : ] = a[ : ] * b[ : ]
         static int __multiply3kernel__( ubase_t * r, ubase_t * temp, ubase_t * a, const int la, ubase_t * b, const int lb );

         // r[ shift : ] = a[ : ] * b
         static int __multiply2kernel__( ubase_t * r, ubase_t * a, const int la, const ubase_t b, const int shift );

   };

}

#endif

