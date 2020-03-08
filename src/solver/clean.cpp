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

#include "solver.h"

/*
    Remove items from list which have a UNIQUE odd power of a prime, i.e. there are no other items in list with an odd power of that prime
    These items cannot be multiplied with other items of list to yield an even power of the prime
*/
std::vector<jans::smooth_number> jans::solver::clean(const std::vector<smooth_number>& list, const uint32_t num_primes)
{
    std::vector<uint8_t> contributes;
    contributes.reserve(list.size());
    for (size_t item = 0U; item < list.size(); ++item){ contributes.push_back(1U); }

    for (uint32_t ip = 0U; ip < num_primes; ++ip)
    {
        uint8_t num_odds = 0U;
        size_t index_last = list.size();
        for (const smooth_number& sn : list)
        {
            for (const prime_factor& pf : sn.factors)
            {
                if ((pf.index == ip) && ((pf.power & 1U) == 1U))
                {
                    ++num_odds;
                    index_last = &sn - &list[0]; // std::vector obeys &v[n] = &v[0] + n
                    if (num_odds >= 2U)
                    {
                        break;
                    }
                }
            }
            if (num_odds >= 2U)
            {
                break;
            }
        }
        if (num_odds == 1U)
        {
            contributes[index_last] = 0U;
        }
    }

    std::vector<smooth_number> relevant;
    for (size_t item = 0; item < list.size(); ++item)
    {
        if (contributes[item] == 1U)
        {
            relevant.push_back(list[item]);
        }
    }

    std::cout << "jans::solver::clear: Removed " << list.size() - relevant.size() << " of the " << list.size() << " vectors with a unique odd prime power." << std::endl;

    return relevant;
}

