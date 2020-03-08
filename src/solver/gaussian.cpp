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
#include <assert.h>

#include "solver.h"

/*
    Retrieve the nullspace via Gaussian elimination
    IMPORTANT: Run jans::solver::clean before passing to this function
*/
std::vector<std::vector<uint32_t>> jans::solver::gaussian(const std::vector<smooth_number>& list, const uint32_t num_primes)
{
    const std::vector<uint32_t> contributions = power_contributions(list, num_primes);

    const uint32_t d_lin = list.size();
    const uint32_t d_pow = contributions.size() + 1U;
    const uint32_t d_row = d_pow + d_lin;

    const uint32_t N_row = (d_row >> 3U) + (((d_row & 7U) > 0U) ? 1U : 0U);
    const uint64_t X_row = static_cast<uint64_t>(N_row);
    uint8_t * matrix = new uint8_t[X_row * d_lin];

    /*
     * matrix = [  d_pow x d_lin  ]
     *          [  -------------  ]
     *          [  d_lin x d_lin  ]
     */

    for (uint32_t lin = 0U; lin < d_lin; ++lin)
    {
        for (uint32_t row = 0U; row < N_row; ++row){ matrix[row + X_row * lin] = 0U; }

        matrix[0U + X_row * lin] |= list[lin].negative ? 1U : 0U;

        uint32_t red = 0U;
        for (const uint32_t ip : contributions)
        {
            for (const prime_factor& pf : list[lin].factors)
            {
                if ((pf.index == ip) && (pf.power & 1U))
                {
                    matrix[((1U + red) >> 3U) + X_row * lin] |= static_cast<uint8_t>(1U) << ((1U + red) & 7U);
                }
            }
            ++red;
        }
        matrix[((d_pow + lin) >> 3U) + X_row * lin] |= static_cast<uint8_t>(1U) << ((d_pow + lin) & 7U);
    }

    uint32_t start = 0;
    for (uint32_t pow = 0U; pow < d_pow; ++pow)
    {
        bool found   = false;
        uint32_t iter = start;
        const uint32_t ix = pow >> 3U;
        const uint8_t  iy = pow  & 7U;
        while ((!found) && (iter < d_lin)){
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
            for (uint32_t vec = iter + 1; vec < d_lin; ++vec){
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

    // Now matrix [0:d_pow, start:d_lin] == 0

    std::cout << "jans::solver::gaussian: Found " << d_lin - start << " vectors in the nullspace." << std::endl;
    std::vector<std::vector<uint32_t>> nullspace;

    for (uint32_t sol = start; sol < d_lin; ++sol)
    {
        std::vector<uint32_t> nullvector;
        for (uint32_t row = d_pow; row < d_row; ++row){
            if ((matrix[(row >> 3U) + X_row * sol] >> (row & 7U)) & static_cast<uint8_t>(1U))
            {
                // item (of list) = row - d_pow
                nullvector.push_back(row - d_pow);
            }
        }
        nullspace.push_back(nullvector);
    }

    delete [] matrix;

    return nullspace;
}

