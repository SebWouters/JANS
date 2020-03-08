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
#include "solver/solver.h"
#include <vector>

namespace jans{

   class sieve{

      public:

         sieve( jans::big_int & num, const ubase_t factorbound, const ubase_t sievespace, const int extra );

         virtual ~sieve();

         static ucarry_t optimal_factorbound( jans::big_int & number );

         void run( jans::big_int & sol_p, jans::big_int & sol_q, const double threshold );

      private:

         ubase_t M;

         jans::big_int target; // N

         // Factor base information

         int num_primes;

         ubase_t * primes;

         ubase_t * roots;

         double * logval;

         // Intermediate sieve results

         int extra;

         int powspace; // num_primes + 1 ( to account for sign Q(x) )

         int required; // powspace + extra_sz

         //int lincount;

         //jans::big_int * xvalues;

         //jans::big_int * pvalues;

         //ubase_t ** powers;

         std::vector<smooth_number> factorization;

         // Helper funcionality

         static int __legendre_symbol__( jans::big_int & num, jans::big_int & p );

         static int __legendre_symbol__( jans::big_int & num, const ubase_t p );

         static int __legendre_symbol__( const ubase_t num, const ubase_t p );

         static ubase_t __power__( const ubase_t num, const ubase_t pow, const ubase_t mod );

         static ubase_t __inv_x_mod_p__( jans::big_int & x, const ubase_t p );

         static ubase_t __root_quadratic_residue__( jans::big_int & num, const ubase_t p );

         static ubase_t __root_quadratic_residue__( const ubase_t num, const ubase_t p );

         bool __extract__( big_int & x, ubase_t * powers ) const;

         // The core routines, in order

         void __startup1__( jans::big_int & mpqs_q );

         void __startup2__( const ubase_t bound );

         bool __check_mpqs_q__( jans::big_int & a, jans::big_int & b, jans::big_int & mpqs_q );

         void __calculate_shifts__( ubase_t * shift1, ubase_t * shift2, jans::big_int & a, jans::big_int & b );

         void __sieve_sumlog__( const ubase_t size, double * sumlog, ubase_t * shift1, ubase_t * shift2 ) const;

         void __check_sumlog__( const ubase_t size, double * sumlog, ubase_t * helper, const double threshold, jans::big_int & a, jans::big_int & b, jans::big_int & mpsqs_q );

         void __factor__(const std::vector<std::vector<uint32_t>>& nullspace, jans::big_int & p, jans::big_int & q);

   };

}

#endif

