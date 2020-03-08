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
    Returns std::vector<nullvector>, with nullvector = { index : XOR_index( space[index] ) == null_vector }
*/
std::vector<std::vector<uint32_t>> jans::gf2solver::gaussian(const std::vector<sparse_vector>& space, const uint32_t basis_size)
{
    const std::vector<uint32_t> contributions = basis_contributions(space, basis_size);

    const uint32_t d_vec = space.size();
    const uint32_t d_bas = contributions.size();
    const uint32_t d_row = d_bas + d_vec;

    const uint32_t N_row = (d_row >> 3U) + (((d_row & 7U) > 0U) ? 1U : 0U);
    const uint64_t X_row = static_cast<uint64_t>(N_row);

    uint8_t * matrix = new uint8_t[X_row * d_vec];
    /*
     * matrix = [  d_bas x d_vec  ]
     *          [  -------------  ]
     *          [  d_vec x d_vec  ]
     */

    for (uint32_t vec = 0U; vec < d_vec; ++vec)
    {
        for (uint32_t row = 0U; row < N_row; ++row)
        {
            matrix[row + X_row * vec] = 0U;
        }
        uint32_t row = 0U;
        const sparse_vector& space_vec = space[vec];
        for (const uint32_t basis_idx : contributions)
        {
            if (std::find(space_vec.begin(), space_vec.end(), basis_idx) != space_vec.end())
            {
                matrix[(row >> 3U) + X_row * vec] |= static_cast<uint8_t>(1U) << (row & 7U);
            }
            ++row;
        }
        matrix[((d_bas + vec) >> 3U) + X_row * vec] |= static_cast<uint8_t>(1U) << ((d_bas + vec) & 7U);
    }

    uint32_t start = 0;
    for (uint32_t bas = 0U; bas < d_bas; ++bas)
    {
        bool found = false;
        uint32_t iter = start;
        const uint32_t ix = bas >> 3U;
        const uint8_t  iy = bas  & 7U;
        while ((!found) && (iter < d_vec)){
            if ((matrix[ix + X_row * iter] >> iy) & static_cast<uint8_t>(1U))
            {
                found = true;
            } else {
                ++iter;
            }
        }
        if (found){
            if (iter != start){
                for (uint32_t row = ix; row < N_row; ++row)
                {
                    // Everything above specific bit has been cleared: bit-wise swapping (also for preceding bits) is OK.
                    uint8_t swap = matrix[row + X_row * iter];
                    matrix[row + X_row * iter ] = matrix[row + X_row * start];
                    matrix[row + X_row * start] = swap;
                }
            }
            #pragma omp parallel for schedule(static)
            for (uint32_t vec = iter + 1; vec < d_vec; ++vec){
                if ((matrix[ix + X_row * vec] >> iy) & static_cast<uint8_t>(1U))
                {
                    for (uint32_t row = ix; row < N_row; ++row){
                        // Everything above specific bit has been cleared: bit-wise XOR (also for preceding bits) is OK.
                        matrix[row + X_row * vec] ^= matrix[row + X_row * start];
                    }
                }
            }
            ++start;
        }
    }

    // Now matrix [0:d_bas, start:d_vec] == 0

    std::cout << "jans::gf2solver::gaussian: Found " << d_vec - start << " vectors in the nullspace." << std::endl;

    std::vector<std::vector<uint32_t>> nullspace;
    for (uint32_t sol = start; sol < d_vec; ++sol)
    {
        std::vector<uint32_t> nullvector;
        for (uint32_t row = d_bas; row < d_row; ++row){
            if ((matrix[(row >> 3U) + X_row * sol] >> (row & 7U)) & static_cast<uint8_t>(1U))
            {
                // item (of list) = row - d_bas
                nullvector.push_back(row - d_bas);
            }
        }
        nullspace.push_back(nullvector);
    }

    delete [] matrix;

    return nullspace;
}

