/*
   JANS: just another number sieve
   Copyright (C) 2018 Sebastian Wouters

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

#include <math.h>
#include <getopt.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <iostream>

#include "big_int.h"
#include "sieve.h"

void print_help(){

std::cout << "\n"
"JANS: just another number sieve\n"
"Copyright (C) 2018 Sebastian Wouters\n"
"\n"
"Usage: jans [OPTIONS]\n"
"\n"
"   INFO\n"
"       JANS is a Multiple Polynomial Quadratic Sieve (MPQS)\n"
"       to factor a large number N. The sieving function is:\n"
"\n"
"           f(x) = ax^2 + 2bx + c with x in [-M, M].\n"
"\n"
"             --> a = q^2 <= sqrt(2N)/M\n"
"             --> q a prime\n"
"             --> q [mod 4] = 3\n"
"\n"
"       The MPQS attempts to factor f(x) over primes p <= F.\n"
"       With sufficient factorizations, the program builds Z\n"
"       congruences\n"
"\n"
"           Y^2 = X^2  [mod N].\n"
"\n"
"       With probability 1 - (1/2)^Z a factor gcd( X-Y, N ) \n"
"       of N is retrieved which is non-trivial.\n"
"\n"
"   ARGUMENTS\n"
"       -N, --number=integer\n"
"              Number to factorize.\n"
"\n"
"       -F, --factorbound=integer\n"
"              Upper bound for factor base primes p <= F.\n"
"\n"
"       -M, --sievespace=integer\n"
"              x in [-M, M] whereby M >= F.\n"
"\n"
"       -Z, --congruences=integer\n"
"              Number of congruences to construct (default 11).\n"
"\n"
"       -T, --threshold=float\n"
"              Threshold for attempting trial division (default 8.0).\n"
"\n"
"       -B, --bits=integer\n"
"              Large integer bit precision. Should be a multiple of " << BASE_UNIT << " (default 1024).\n"
"\n"
"       -v, --version\n"
"              Print the version.\n"
"\n"
"       -h, --help\n"
"              Display this help.\n"
"\n"
"   EXAMPLES\n"
"       $ jans -N 61421677127643670816789 -F 1350 -M 10000\n"                                           // Linspace = 138
"       $ jans -N 11111111111111111111111111111 -F 3652 -M 10000\n"                                     // Linspace = 262
"       $ jans -N 4205940337640636327774357502033476724941 -F 25499 -M 40000\n"                         // Linspace = 1379
"       $ jans -N 392318858461667547569595655490009919272404068553904357377 -F 295292 -M 300000\n"      // Linspace = 12837
"       $ jans -N 3571082522473766674484304975778527401895200115726120795842576355509746402614775567 -F 1240000 -M 1240000 -T 12.0\n"
//"       $ jans -N 3571082522473766674484304975778527401895200115726120795842576355509746402614775567 -F 6434337 -M 6500000 -T 15.0\n" // Linspace = 220164
"       $ jans -N 1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139 -F 1240000 -M 1240000\n"
" " << std::endl;

}

int main( int argc, char ** argv ){

   ubase_t factorbound = 0;
   ubase_t sievespace  = 0;
   ubase_t congruences = 11;
   double  threshold   = 8.0;
   ubase_t bits        = 1024;

   std::string temp_str;
   long long temp_int = -1;

   struct option long_options[] =
   {
      {"number",      required_argument, 0, 'N'},
      {"factorbound", required_argument, 0, 'F'},
      {"sievespace",  required_argument, 0, 'M'},
      {"congruences", required_argument, 0, 'Z'},
      {"threshold",   required_argument, 0, 'T'},
      {"bits",        required_argument, 0, 'B'},
      {"version",     no_argument,       0, 'v'},
      {"help",        no_argument,       0, 'h'},
      {0, 0, 0, 0}
   };

   int option_index = 0;
   int c;
   while (( c = getopt_long( argc, argv, "hvN:F:M:Z:T:B:", long_options, &option_index )) != -1 ){
      switch( c ){
         case 'h':
         case '?':
            print_help();
            return 0;
            break;
         case 'v':
            std::cout << "jans version 0.3.0" << std::endl;
            return 0;
            break;
         case 'N':
            temp_str = optarg;
            break;
         case 'F':
            temp_int = atol( optarg );
            if ( temp_int < 1 ){
               std::cerr << "   Error: -F, --factorbound should be a non-zero positive integer" << std::endl;
               return 7;
            }
            factorbound = temp_int;
            break;
         case 'M':
            temp_int = atol( optarg );
            if ( temp_int < 1 ){
               std::cerr << "   Error: -M, --sievespace should be a non-zero positive integer" << std::endl;
               return 7;
            }
            sievespace = temp_int;
            break;
         case 'Z':
            temp_int = atol( optarg );
            if ( temp_int < 1 ){
               std::cerr << "   Error: -Z, --congruences should be a non-zero positive integer" << std::endl;
               return 7;
            }
            congruences = temp_int;
            break;
         case 'T':
            threshold = atof( optarg );
            if ( threshold < 0.0 ){
               std::cerr << "   Error: -T, --threshold should be a non-negative float" << std::endl;
               return 7;
            }
            break;
         case 'B':
            temp_int = atol( optarg );
            if ( ( temp_int < 1 ) || ( ( temp_int % BASE_UNIT ) != 0 ) || ( temp_int < BASE_UNIT ) ){
               std::cerr << "   Error: -B, --bits should be a non-zero positive multiple of " << BASE_UNIT << std::endl;
               return 7;
            }
            bits = temp_int;
            break;
      }
   }

   jans::big_int::set_num_block( bits / BASE_UNIT );
   jans::big_int::sanity_check();
   jans::big_int number;

   const int text_size = temp_str.length();
   if ( text_size > 0 ){
      for ( int cnt = 0; cnt < text_size; cnt++ ){
         const ubase_t digit = jans::big_int::convert_c2i( temp_str[ cnt ] );
         if ( digit > 9 ){
            std::cerr << "   Error: -N, --number should only contain digits 0 to 9" << std::endl;
            return 11;
         }
      }
      std::cout << "Parsed integer to factor = " << temp_str << std::endl;
      number.read( temp_str, 10 );
   }

   if ( jans::big_int::equal( number, 0 ) ){
      std::cerr << "   Error: -N, --number should be specified" << std::endl;
      return 11;
   }

   if ( factorbound == 0 ){
      const ucarry_t optbound = jans::sieve::optimal_factorbound( number );
      std::cerr << "   Error: -F, --factorbound should be specified" << std::endl;
      std::cerr << "   A well-educated suggestion is -F " << optbound << std::endl;
      return 11;
   }

   if ( sievespace == 0 ){
      std::cerr << "   Error: -M, --sievespace should be specified" << std::endl;
      return 11;
   }

   if ( sievespace < factorbound ){
      std::cerr << "   Error: -M, --sievespace should not be smaller than -F, --factorfound" << std::endl;
      return 11;
   }

   std::cout << "Parsed command: " << std::endl;
   std::cout << "./jans -N " << number.write( 10 )
                   << " -F " << factorbound
                   << " -M " << sievespace
                   << " -Z " << congruences
                   << " -T " << threshold
                   << " -B " << bits << std::endl;

   jans::big_int sol_p;
   jans::big_int sol_q;

   jans::sieve mysieve( number, factorbound, sievespace, congruences );
   mysieve.run( sol_p, sol_q, threshold );

   std::cout << "Factored N = P x Q with" << std::endl;
   std::cout << "      N = " << number.write( 10 ) << std::endl;
   std::cout << "      P = " <<  sol_p.write( 10 ) << std::endl;
   std::cout << "      Q = " <<  sol_q.write( 10 ) << std::endl;

   return 0;

}

