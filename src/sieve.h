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

#ifndef JANS_SIEVE
#define JANS_SIEVE

#include "big_int.h"

namespace jans{

   class sieve{

      public:

         sieve( jans::big_int & num, const ubase_t B, const ubase_t M, const int extra );

         virtual ~sieve();

         void run( jans::big_int & p, jans::big_int & q, const double grace );

      private:

         ubase_t M;

         jans::big_int target; // N

         jans::big_int mpqs_p; // number near ( 2N )^0.25 / sqrt( M )

         // Factor base information

         int num_primes;

         ubase_t * primes;

         ubase_t * roots;

         double * logval;

         // Intermediate sieve results

         int extra_sz;

         int powspace; // num_primes + 1 ( to account for sign Q(x) )

         int linspace; // powspace + extra_sz

         int lincount;

         jans::big_int * xvalues;

         jans::big_int * pvalues;

         ubase_t * powers;

         // Helper funcionality

         void __startup1__();

         void __startup2__( const ubase_t bound );

         static int __legendre_symbol__( jans::big_int & num, jans::big_int & p );

         static int __legendre_symbol__( jans::big_int & num, const ubase_t p );

         static int __legendre_symbol__( const ubase_t num, const ubase_t p );

         static ubase_t __power__( const ubase_t num, const ubase_t pow, const ubase_t mod );

         static ubase_t __inv_x_mod_p__( jans::big_int & x, const ubase_t p );

         static ubase_t __root_quadratic_residue__( jans::big_int & num, const ubase_t p );

         static ubase_t __root_quadratic_residue__( const ubase_t num, const ubase_t p );

         bool __extract__( big_int & x, ubase_t * powers ) const;

         void __check_sumlog__( const ubase_t size, double * sumlog, ubase_t * helper, const double grace, jans::big_int & a, jans::big_int & b );

         static void __solve_gaussian__( unsigned char * out, ubase_t * vectors, const ubase_t d_pow, const ubase_t d_lin );

         void __factor__( unsigned char * helper, jans::big_int & p, jans::big_int & q );

         void __next_mpqs_p__( jans::big_int & a, jans::big_int & b );

         void __calculate_shifts__( ubase_t * shift1, ubase_t * shift2, jans::big_int & a, jans::big_int & b );

         void __sieve_sumlog__( const ubase_t size, double * sumlog, ubase_t * shift1, ubase_t * shift2 ) const;

   };

}

#endif

