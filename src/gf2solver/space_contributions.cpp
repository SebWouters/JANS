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

#include <stdio.h>
#include <iostream>
#include <algorithm>

#include "gf2solver.h"

/*
    Returns the indices of sparse_vectors in space which can contribute to the nullspace of space
*/
std::vector<uint32_t> jans::gf2solver::space_contributions(const std::vector<sparse_vector>& space, const uint32_t basis_size)
{
    std::vector<uint8_t> contributes;
    contributes.reserve(space.size());
    for (uint32_t space_idx = 0U; space_idx < space.size(); ++space_idx){ contributes.push_back(1U); }

    for (uint32_t basis_idx = 0U; basis_idx < basis_size; ++basis_idx)
    {
        uint8_t num_odd = 0U;
        size_t space_idx_last = space.size();
        for (const sparse_vector& vec : space)
        {
            if (std::find(vec.begin(), vec.end(), basis_idx) != vec.end())
            {
                ++num_odd;
                // std::vector obeys &v[n] = &v[0] + n
                space_idx_last = &vec - &space[0];
            }
            if (num_odd >= 2U)
            {
                break;
            }
        }
        if (num_odd == 1U)
        {
            contributes[space_idx_last] = 0U;
        }
    }

    std::vector<uint32_t> relevant;
    relevant.reserve(space.size());
    for (uint32_t space_idx = 0; space_idx < space.size(); ++space_idx)
    {
        if (contributes[space_idx])
        {
            relevant.push_back(space_idx);
        }
    }

    std::cout << "jans::gf2solver::space_contributions: " << space.size() - relevant.size() << " of the " << space.size() << " space vectors cannot contribute to the nullspace." << std::endl;

    return relevant;
}

