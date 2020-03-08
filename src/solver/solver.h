/*
   JANS: just another number sieve
   Copyright (C) 2018-2020 Sebastian Wouters

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

#pragma once

#include "../big_int.h"
#include <vector>

namespace jans {

    typedef struct
    {
        uint32_t index; // of the prime
        uint32_t power; // of the prime: idea: only retain if power==1?
    } prime_factor;

    typedef struct
    {
        jans::big_int xval; // Best to pull out?
        jans::big_int pval;
        std::vector<prime_factor> factors;
        bool negative;
    } smooth_number;

    namespace solver {

         std::vector<smooth_number> clean(const std::vector<smooth_number>& list, const uint32_t num_primes);

         std::vector<uint32_t> power_contributions(const std::vector<smooth_number>& list, const uint32_t num_primes);

         std::vector<std::vector<uint32_t>> gaussian(const std::vector<smooth_number>& list, const uint32_t num_primes);

    }
}


