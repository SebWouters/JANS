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

#include <vector>

namespace jans {

    /*
        Solver for a vector space GF(2) x GF(2) x GF(2) x ... x GF(2):
          - full_vector is never represented
          - basis_size  is the length of a full_vector
          - basis_index is an index vis-a-vis a full_vector (0 <= basis_index < basis_size)
     */

    // A sparse_vector contains all basis_index for which full_vector[basis_index] == 1
    using sparse_vector = std::vector<uint32_t>;

    namespace gf2solver
    {
        std::vector<uint32_t>   space_contributions(const std::vector<sparse_vector>& space, const uint32_t basis_size);
        std::vector<uint32_t>   basis_contributions(const std::vector<sparse_vector>& space, const uint32_t basis_size);
        std::vector<std::vector<uint32_t>> gaussian(const std::vector<sparse_vector>& space, const uint32_t basis_size);

        std::vector<uint8_t> random(const uint32_t size);
        std::vector<uint8_t> zeroes(const uint32_t size);
        uint8_t dotproduct(const std::vector<uint8_t>& first, const std::vector<uint8_t>& second);
        bool update_x(std::vector<uint8_t>& x, const std::vector<uint8_t>& v, const std::vector<uint8_t>& Av, const std::vector<uint8_t>& b);

        std::vector<uint8_t>  first_half(const std::vector<sparse_vector>& space, const std::vector<uint8_t>& input, const uint32_t basis_size);
        std::vector<uint8_t> second_half(const std::vector<sparse_vector>& space, const std::vector<uint8_t>& half,  const uint32_t basis_size);
        std::vector<uint8_t>      matvec(const std::vector<sparse_vector>& space, const std::vector<uint8_t>& input, const uint32_t basis_size);
    }
}


