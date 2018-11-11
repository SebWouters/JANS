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

#define __11111111__ ( ~((unsigned int)(0)) )

//#define SET_BIT(A,k) (   A[ ( ( k ) / ( sizeof(unsigned int) * CHAR_BIT ) ) ] |=  ( 1U << ( ( k ) % ( sizeof(unsigned int) * CHAR_BIT ) ) ) )
//#define DEL_BIT(A,k) (   A[ ( ( k ) / ( sizeof(unsigned int) * CHAR_BIT ) ) ] &= ~( 1U << ( ( k ) % ( sizeof(unsigned int) * CHAR_BIT ) ) ) )
//#define GET_BIT(A,k) ( ( A[ ( ( k ) / ( sizeof(unsigned int) * CHAR_BIT ) ) ] >>          ( ( k ) % ( sizeof(unsigned int) * CHAR_BIT ) ) ) & 1U )

namespace jans{

   class big_int{

      public:

         big_int();

         virtual ~big_int();

         void copy( big_int & tocopy );

         void set( const ubase_t value );

         ubase_t get_blk( const int i ); // data[ i ]

         void read( const std::string number, const ubase_t base );

         std::string write( const ubase_t base ) const;

         static bool equal( big_int & n1, big_int & n2 ); // ( n1 == n2 )

         static bool equal( big_int & n1, const ubase_t n2 );

         static bool smaller( big_int & n1, big_int & n2 ); // ( n1 < n2 )

         static void sanity_check();

         // Generic math operations

         static void sum( big_int & res, big_int & a, big_int & b );

         static void sum( big_int & res, big_int & a, const ubase_t b );

         static void diff( big_int & res, big_int & a, big_int & b );

         static void diff( big_int & res, big_int & a, const ubase_t b );

         static void prod( big_int & res, big_int & a, big_int & b );

         static void prod( big_int & res, big_int & a, const ubase_t b );

         static void div( big_int & q, big_int & r, big_int & n, big_int & d );

         static ubase_t div( big_int & q, big_int & n, const ubase_t d ); // Returns remainder

         static double logarithm( big_int & x );

         // Number theory operations

         static void gcd( big_int & res, big_int & a, big_int & b );

         static ubase_t gcd( big_int & a, const ubase_t b ); // Returns gcd( a, b )

         static void ceil_sqrt( big_int & res, big_int & n );

         static void xx_min_num( big_int & res, big_int & x, big_int & num );

         static ubase_t extract_pow_2( big_int & x );

         static ubase_t extract_pow_p( big_int & x, const ubase_t p );

      private:

         ubase_t data[ NUM_BLOCK ];

         int lead; // Upper bound for loops over the blocks: lead = 1 + max{i}( data[ i ] != 0 )

         /******************************
          *  Private static functions  *
          ******************************/

         static void __clear__( ubase_t * a );

         static void __copy__( ubase_t * r, const ubase_t * a );

         static int __compare__( const ubase_t * a, const ubase_t * b );

         static void __shift_up__( ubase_t * a, const int k );

         static void __shift_down__( ubase_t * a, const int k );

         /********
          *  IO  *
          ********/

         static const char __conversion__[ 16 ];

         static ubase_t __convert_c2i__( const char c );

         static char __convert_i2c__( const ubase_t c );

         /*************************
          *  Basic math routines  *
          *************************/

         // r = a + b; if (( r != a ) && ( r != b )){ __clear__( r ); }
         static int __sum3set__( ubase_t * r, ubase_t * a, const int la, ubase_t * b, const int lb );

         // r = r + b
         static int __sum1__( ubase_t * r, const int lr, const ubase_t b );

         // r = a - b
         static int __diff3set__( ubase_t * r, ubase_t * a, ubase_t * b );

         // r = r - b
         static int __diff1__( ubase_t * r, const ubase_t b );

         // r = a * b; calls __clear__( r )
         static int __mult3set__( ubase_t * r, ubase_t * a, const int la, ubase_t * b, const int lb );

         // r += a * b
         static int __mult3add__( ubase_t * r, const int lr, ubase_t * a, const int la, ubase_t * b, const int lb );

         // r[ shift : ] = a[ : ] * b; calls __clear__( r )
         static int __mult2set__( ubase_t * r, ubase_t * a, const int la, const ubase_t b, const int shift );

         // r[ shift : ] += b * a[ : ]; this would be lapack "axpy" with a shift
         static int __mult2add__( ubase_t * r, const int lr, ubase_t * a, const int la, const ubase_t b, const int shift );

         // r = r * b
         static int __scal1__( ubase_t * r, const int lr, const ubase_t b );

         // Solves for n = q * d + r, with r < d; whereby initially (q, lq) contains (n, ln).
         static ubase_t __divide__( ubase_t * q, int & lq, const ubase_t div );

         // Solves for n = q * d + r, with r < d; whereby initially (r, lr) contains (n, ln).
         static void __divide__( ubase_t * q, int & lq, ubase_t * r, int & lr, ubase_t * d, const int ld );

         // Solves for res = gcd( a, b ); a >= b; destroys a & b in the proces
         static int __gcd__( ubase_t * res, ubase_t * a, const int la, ubase_t * b, const int lb );

         // Solves for d = ceil( sqrt( num ) )
         static int __ceil_sqrt__( ubase_t * d, ubase_t * num, const int ln );

   };

}

#endif

