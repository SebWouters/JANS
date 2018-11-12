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

         sieve( const ubase_t bound, jans::big_int & num, const int extra );

         virtual ~sieve();

         void run( jans::big_int & p, jans::big_int & q, const int blk_size, const double grace );

      private:

         ubase_t bound;

         jans::big_int target; // N

         // Factor base information

         int num_primes;

         ubase_t * primes;

         ubase_t * roots;

         double * logval;

         // Intermediate sieve results

         int extra_sz;

         int linspace; // num_primes + extra_sz

         int lincount;

         jans::big_int * xvalues;

         ubase_t * powers;

         // Helper funcionality

         void __startup__();

         static int __legendre_symbol__( jans::big_int & num, const ubase_t p );

         static int __legendre_symbol__( const ubase_t num, const ubase_t p );

         static ubase_t __power__( const ubase_t num, const ubase_t pow, const ubase_t mod );

         static ubase_t __root_quadratic_residue__( jans::big_int & num, const ubase_t p );

         static ubase_t __root_quadratic_residue__( const ubase_t num, const ubase_t p );

         bool __extract__( big_int & x, ubase_t * powers ) const;

         void __sieving_grace__( const int blk_size, const double grace );

         void __solve_gaussian__( unsigned char * helper );

         void __factor__( unsigned char * helper, jans::big_int & p, jans::big_int & q );

   };

}

#endif

