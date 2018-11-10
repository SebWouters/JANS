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

         sieve( const ubase_t bound, big_int & num );

         virtual ~sieve();

      private:

         ubase_t bound;

         big_int n;

         int num_primes;

         ubase_t * primes;

         ubase_t * roots;

         double * logval;

         void __startup__();

         static int __legendre_symbol__( big_int & num, const ubase_t d_prime );

   };

}

#endif

