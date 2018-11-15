JANS: just another number sieve
===============================

Copyright (C) 2018 Sebastian Wouters <sebastianwouters@gmail.com>

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

Information
-----------

JANS is a Multiple Polynomial Quadratic Sieve (MPQS)
to factor a large number N. The sieving function is:

    f(x) = ax^2 + 2bx + c with x in [-M, M].

      --> a = q^2 <= sqrt(2N)/M
      --> q a prime
      --> q [mod 4] = 3

The MPQS attempts to factor f(x) over primes p <= F.
With sufficient factorizations, the program builds Z
congruences

    Y^2 = X^2  [mod N].

With probability 1 - (1/2)^Z a factor gcd( X-Y, N )
of N is retrieved which is non-trivial.

Arguments
---------

    -N, --number=integer
           Number to factorize.

    -F, --factorbound=integer
           Upper bound for factor base primes p <= F.

    -M, --sievespace=integer
           x in [-M, M] whereby M >= F.

    -Z, --congruences=integer
           Number of congruences to construct (default 11).

    -T, --threshold=float
           Threshold for attempting trial division (default 8.0).

    -B, --bits=integer
           Large integer bit precision. Should be a multiple of 256 (default 1024).

    -v, --version
           Print the version.

    -h, --help
           Display this help.

Examples
--------

    $ jans -N 61421677127643670816789 -F 1350 -M 10000
    $ jans -N 11111111111111111111111111111 -F 3652 -M 10000
    $ jans -N 4205940337640636327774357502033476724941 -F 25499 -M 40000

