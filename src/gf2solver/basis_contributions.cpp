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
    Returns the basis_indices (0 <= basis_index < basis_size) which can contribute to the nullspace of space
*/
std::vector<uint32_t> jans::gf2solver::basis_contributions(const std::vector<sparse_vector>& space, const uint32_t basis_size)
{
    std::vector<uint32_t> contributions;
    contributions.reserve(basis_size);

    for (uint32_t basis_idx = 0U; basis_idx < basis_size; ++basis_idx)
    {
        for (const sparse_vector& vec : space){
            if (std::find(vec.begin(), vec.end(), basis_idx) != vec.end())
            {
                contributions.push_back(basis_idx);
                break;
            }
        }
    }

    std::cout << "jans::solver::basis_contributions: " << basis_size - contributions.size() << " of the " << basis_size << " basis indices cannot contribute to the nullspace." << std::endl;

    return contributions;
}

