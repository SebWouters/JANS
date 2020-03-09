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
#include <random>

#include "gf2solver.h"

/*
    Return a random vector: TODO: should be bit-random in range [0, 2**NBIT[
*/
std::vector<uint8_t> jans::gf2solver::random(const uint32_t size)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);

    std::vector<uint8_t> random_vector;
    random_vector.reserve(size);
    for (uint32_t iter = 0U; iter < size; ++iter)
    {
        random_vector.push_back(dis(gen));
        //if (dis(gen))
        //{
        //    random_vector.push_back(iter);
        //}
    }
    return random_vector;
}

/*
    Return an empty vector
*/
std::vector<uint8_t> jans::gf2solver::zeroes(const uint32_t size)
{
    std::vector<uint8_t> result;
    result.reserve(size);
    for (uint32_t iter = 0U; iter < size; ++iter)
    {
        result.push_back(0U);
    }
    return result;
}

/*
    Dot product TODO: should become bitwise matrix
*/
uint8_t jans::gf2solver::dotproduct(const std::vector<uint8_t>& first, const std::vector<uint8_t>& second)
{
    assert(first.size() == second.size());
    uint8_t result = 0U;
    for (uint32_t idx = 0U; idx < first.size(); ++idx)
    {
        result ^= (first[idx] ^ second[idx]);
    }
}

/*
    Return half[basis_idx, bit] = sum_{space_idx} space[space_idx][basis_idx] * input[space_idx, bit]
*/
std::vector<uint8_t> jans::gf2solver::first_half(const std::vector<sparse_vector>& space, const std::vector<uint8_t>& input, const uint32_t basis_size)
{
    assert(input.size() == space.size());
    std::vector<uint8_t> half = zeroes(basis_size);

    #pragma omp parallel
    {
        std::vector<uint8_t> temp = zeroes(basis_size);

        #pragma omp for schedule(static)
        for (uint32_t space_idx = 0U; space_idx < space.size(); ++space_idx)
        {
            const sparse_vector& vec = space[space_idx];
            for (const uint32_t basis_idx : vec)
            {
                temp[basis_idx] ^= input[space_idx];
            }
        }

        #pragma omp critical
        {
            for (uint32_t basis_idx = 0; basis_idx < basis_size; ++basis_idx)
            {
                half[basis_idx] ^= temp[basis_idx];
            }
        }
    }

    return half;
}


/*
    Return output[space_idx, bit] = sum_{basis_idx} space[space_idx][basis_idx] * half[basis_idx, bit]
*/
std::vector<uint8_t> jans::gfsolver::second_half(const std::vector<sparse_vector>& space, const std::vector<uint8_t>& half, const uint32_t basis_size)
{
    assert(half.size() == basis_size);
    std::vector<uint8_t> output = zeroes(space.size());

    #pragma omp parallel for schedule(static)
    for (uint32_t space_idx = 0; space_idx < space.size(); ++space_idx)
    {
        const sparse_vector& vec = space[space_idx];
        for (const uint32_t basis_idx : vec)
        {
            output[space_idx] ^= half[basis_idx];
        }
    }
    return output;
}

/*
    Return output[space_idx, bit] = space^T * space * input[space_idx, bit]
*/
std::vector<uint8_t> jans::gf2solver::matvec(const std::vector<sparse_vector>& space, const std::vector<uint8_t>& input, const uint32_t basis_size)
{
    std::vector<uint8_t> half   =  first_half(space, input, basis_size);
    std::vector<uint8_t> output = second_half(space, half,  basis_size);
    return output;
}

/*
    Update x: TODO: bitwise
*/
bool jans::gf2solver::update_x(std::vector<uint8_t>& x, const std::vector<uint8_t>& v, const std::vector<uint8_t>& Av, const std::vector<uint8_t>& b)
{
    const uint8_t denominator = dotproduct(v, Av);
    if (denominator == 0U)
    {
        return true; // Found solution!
    }
    const uint8_t numerator = dotproduct(v, b);
    assert(x.size() == v.size());
    if (numerator)
    {
        for (uint32_t idx = 0; idx < v.size(); ++idx)
        {
            x[idx] ^= v[idx];
        }
    }
    return false;
}

/*
    Compute new v and shift down: TODO: bitwise
*/
void jans::gf2solver::update_v(std::vector<uint8_t>& Av1, std::vector<uint8_t>& v1, std::vector<uint8_t>& Av2, std::vector<uint8_t>& v2)
{
    const uint8_t c1 = dotproduct(Av1, Av1); // denominator 1 as algo not yet stopped :-)
    const uint8_t c2 = dotproduct(Av2, Av1); // denominator 1 as algo not yet stopped :-)

    std::vector<uint8_t> v0 = std::copy(Av1);
    if (c1)
    {
        for (uint32_t idx = 0U; idx < v0.size(); ++idx)
        {
            v0[idx] ^= v1[idx];
        }
    }
    if (c2)
    {
        for (uint32_t idx = 0U; idx < v0.size(); ++idx)
        {
            v0[idx] ^= v2[idx];
        }
    }
    v2.swap(v1);
    Av2.swap(Av1);
    v1.swap(v0);
}

/*
    Returns std::vector<nullvector>, with nullvector = { index : XOR_index( space[index] ) == null_vector }
    TODO: bitwise and block_lanczos...
*/
std::vector<std::vector<uint32_t>> jans::gf2solver::lanczos(const std::vector<sparse_vector>& space, const uint32_t basis_size)
{
    const std::vector<uint32_t> contributions = basis_contributions(space, basis_size);

    const std::vector<uint8_t> y_rand = random(space.size());
    const std::vector<uint8_t> b_vec  = matvec(space, y_rand, basis_size);

    std::vector<uint8_t>   x = zeroes(space.size());
    std::vector<uint8_t>  v1 = std::copy(b_vec);
    std::vector<uint8_t> Av1 = matvec(space, vi, basis_size);
    std::vector<uint8_t>  v2 = zeroes(space.size());
    std::vector<uint8_t> Av2 = zeroes(space.size());
    while (!update_x(x, v1, Av1, b))
    {
        update_v(Av1, v1, Av2, v2);
        Av1 = matvec(space, v1, basis_size);
    }

    // Two null vectors in the end:
    v1;
    for (uint32_t idx = 0; idx < v2.size(); ++idx)
    {
        v2[idx] = x[idx] ^ y_rand[idx];
    }

    Here be test whether in nullspace :-)

    std::cout << "jans::gf2solver::lanczos: Found 2 vectors in the nullspace." << std::endl;

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

