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
    Return an std::vector with the prime indices which can contribute
*/
std::vector<uint32_t> jans::solver::power_contributions(const std::vector<smooth_number>& list, const uint32_t num_primes)
{
    std::vector<uint32_t> contributions;
    contributions.reserve(num_primes);

    for (uint32_t ip = 0U; ip < num_primes; ++ip)
    {
        uint8_t num_odds = 0U;
        for (const smooth_number& sn : list){
            for (const prime_factor& pf : sn.factors)
            {
                if ((pf.index == ip) && ((pf.power & 1U) == 1U))
                {
                    ++num_odds;
                    break;
                }
            }
            if (num_odds >= 1U)
            {
                break;
            }
        }
        if (num_odds >= 1U)
        {
            contributions.push_back(ip);
        }
    }

    std::cout << "jans::solver::power_contributions: Removed " << num_primes - contributions.size() << " of the " << num_primes << " primes with only even powers." << std::endl;

    return contributions;
}

