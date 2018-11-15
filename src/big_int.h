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

#define BLOCK_BIT ( sizeof( unsigned int ) * CHAR_BIT )
#define BASE_UNIT 256

#define __11111111__ ( ~( ( unsigned int )( 0 ) ) )

namespace jans{

   class big_int{

      public:

         big_int();

         virtual ~big_int();

         void copy( big_int & tocopy );

         void copy( const ubase_t value );

         static void sanity_check();

         static void set_num_block( const int factor );

         // IO

         void read( const std::string number, const ubase_t base );

         std::string write( const ubase_t base ) const;

         static ubase_t convert_c2i( const char c );

         static char convert_i2c( const ubase_t c );

         static long double i2f( big_int & x );

         static void f2i( big_int & x, const long double number );

         // Basic operations

         static bool equal( big_int & n1, big_int & n2 ); // ( n1 == n2 )

         static bool equal( big_int & n1, const ubase_t n2 );

         static bool smaller( big_int & n1, big_int & n2 ); // ( n1 < n2 )

         static bool smaller( big_int & n1, const ubase_t n2 );

         static bool even( big_int & n1 );

         // Generic math operations

         static void plus( big_int & res, const ubase_t val ); // res += val

         static void minus( big_int & res, const ubase_t val ); // res -= val

         static void sum( big_int & res, big_int & a, big_int & b );

         static void sum( big_int & res, big_int & a, const ubase_t b );

         static void diff( big_int & res, big_int & a, big_int & b );

         static void diff( big_int & res, big_int & a, const ubase_t b );

         static void prod( big_int & res, big_int & a, big_int & b );

         static void prod( big_int & res, big_int & a, const ubase_t b );

         static void div( big_int & q, big_int & r, big_int & n, big_int & d );

         static ubase_t div( big_int & q, big_int & n, const ubase_t d ); // Returns remainder

         // Number theory operations

         static void gcd( big_int & res, big_int & a, big_int & b );

         static void power( big_int & res, big_int & base, big_int & expo, big_int & mod );

         static void prodmod( big_int & q, big_int & r, big_int & a, big_int & b, big_int & m ); // { q, r } = { ( a * b ) / m, ( a * b ) % m }

         static void ceil_sqrt( big_int & res, big_int & n );

         static ubase_t extract_pow_p( big_int & x, const ubase_t p );

         static ubase_t random_ubase_t();

         static void randomize( big_int & n );

         static bool miller_rabin( big_int & n, const ubase_t attempts );

      private:

         ubase_t * data;

         int lead; // Upper bound for loops over the blocks: lead = 1 + max{i}( data[ i ] != 0 )

         static int NUM_BLOCK;

         static bool nb_set;

         // Basic functionality

         static void __clear__( ubase_t * a );

         static void __copy__( ubase_t * r, const ubase_t * a );

         static int __compare__( const ubase_t * a, const ubase_t * b );

         static void __shift_up__( ubase_t * a, const int k );

         static void __shift_down__( ubase_t * a, const int k );

         static const char __conversion__[ 16 ];

         // Internal math routines

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
         static int __ceil_sqrt__( ubase_t * d, const int ld, ubase_t * num, const int ln );

   };

}

#endif

