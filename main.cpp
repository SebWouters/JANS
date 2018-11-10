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

#include "big_int.h"
#include <stdio.h>
#include <iostream>

bool test1( const bool print ){

   std::string rsa_100_n_bin = "101100100011010101100110101111010001111100100000011010101100110111001001011011010001110010101111100100000101111110001110111111011110101011100001010100001110011010111101110010011011101101001111011111110111110110011001001000100111010001010101011101110000001011011101110001110001111010010100001110111101111100010111100101100011111011";
   std::string rsa_100_p_bin = "110011111101111010100000111010110100110101010001111011000011000000000100110101001011001111101101100110011011110011100011000111100110101010000000111110010010011110111";
   std::string rsa_100_q_bin = "110110110111100010100000111111001100011101110101101100001100110111100000000110110000000100010000000010110000010100101111101110101111010100000011111001101111100011101";

   std::string rsa_100_n_dec = "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139";
   std::string rsa_100_p_dec = "37975227936943673922808872755445627854565536638199";
   std::string rsa_100_q_dec = "40094690950920881030683735292761468389214899724061";

   std::string rsa_100_n_hex = "2c8d59af47c81ab3725b472be417e3bf7ab85439af726ed3dfdf66489d155dc0b771c7a50ef7c5e58fb";
   std::string rsa_100_p_hex = "19fbd41d69aa3d86009a967db3379c63cd501f24f7";
   std::string rsa_100_q_hex = "1b6f141f98eeb619bc0360220160a5f75ea07cdf1d";

   jans::big_int n1; n1.read( rsa_100_n_bin, 2  );
   jans::big_int n2; n2.read( rsa_100_n_dec, 10 );
   jans::big_int n3; n3.read( rsa_100_n_hex, 16 );

   const bool equal12 = jans::big_int::equal( n1, n2 );
   const bool equal13 = jans::big_int::equal( n1, n3 );

   std::string n_bin = n3.write( 2  ); const bool equal_bin = ( n_bin.compare( rsa_100_n_bin ) == 0 );
   std::string n_dec = n1.write( 10 ); const bool equal_dec = ( n_dec.compare( rsa_100_n_dec ) == 0 );
   std::string n_hex = n2.write( 16 ); const bool equal_hex = ( n_hex.compare( rsa_100_n_hex ) == 0 );

   if ( print ){

      std::cout << "Test 1: reading and writing of numbers in bin, dec & hex" << std::endl;
      std::cout << "rsa_100_n_dec = " << rsa_100_n_dec << std::endl;
      std::cout << "n_dec         = " << n_dec << std::endl;
      std::cout << "Equal(n1,n2) = " << equal12 << std::endl;
      std::cout << "Equal(n1,n3) = " << equal13 << std::endl;
      std::cout << "Equal(bin)   = " << equal_bin << std::endl;
      std::cout << "Equal(dec)   = " << equal_dec << std::endl;
      std::cout << "Equal(hex)   = " << equal_hex << std::endl;

   }

}

bool test2( const bool print ){

   std::string rsa_100_n_dec = "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139";
   std::string rsa_100_q_dec = "40094690950920881030683735292761468389214899724061";
   std::string rsa_100_p_dec = "37975227936943673922808872755445627854565536638199";
   std::string str_sum       = "78069918887864554953492608048207096243780436362260";
   std::string str_diff      =  "2119463013977207107874862537315840534649363085862";
   std::string str_quot_n_s  = "19503094784939146331495321820061269240767748307052";
   std::string str_rem_n_s   = "65185044735636719643690537913604470077414307348619";

   jans::big_int n; n.read( rsa_100_n_dec, 10 );
   jans::big_int q; q.read( rsa_100_q_dec, 10 );
   jans::big_int p; p.read( rsa_100_p_dec, 10 );
   jans::big_int s; s.read( str_sum,       10 );
   jans::big_int d; d.read( str_diff,      10 );
   jans::big_int x; x.read( str_quot_n_s,  10 );
   jans::big_int y; y.read( str_rem_n_s,   10 );
   jans::big_int z; // Initializes to zero

   jans::big_int sum;
   jans::big_int diff;
   jans::big_int prod;
   jans::big_int quot;
   jans::big_int rem;

   jans::big_int::sum(   sum,      p, q ); const bool eq_sum   = jans::big_int::equal(  sum, s );
   jans::big_int::diff( diff,      q, p ); const bool eq_diff  = jans::big_int::equal( diff, d );
   jans::big_int::prod( prod,      p, q ); const bool eq_prod  = jans::big_int::equal( prod, n );
   jans::big_int::div(  quot, rem, n, p ); const bool eq_quot1 = jans::big_int::equal( quot, q );
                                           const bool eq_rem1  = jans::big_int::equal(  rem, z );
   jans::big_int::div(  quot, rem, n, s ); const bool eq_quot2 = jans::big_int::equal( quot, x );
                                           const bool eq_rem2  = jans::big_int::equal(  rem, y );

   if ( print ){

      std::cout << "Test2: Basic math operations: addition, subtraction, multiplication, division, modulo" << std::endl;
      std::cout << "Equal(sum)   = " << eq_sum   << std::endl;
      std::cout << "Equal(diff)  = " << eq_diff  << std::endl;
      std::cout << "Equal(prod)  = " << eq_prod  << std::endl;
      std::cout << "Equal(quot1) = " << eq_quot1 << std::endl;
      std::cout << "Equal(rem1)  = " << eq_rem1  << std::endl;
      std::cout << "Equal(quot2) = " << eq_quot2 << std::endl;
      std::cout << "Equal(rem2)  = " << eq_rem2  << std::endl;

   }

   return ( eq_sum && eq_diff && eq_prod && eq_quot1 && eq_rem1 && eq_quot2 && eq_rem2 );

}

int main()
{

   jans::big_int::sanity_check();
   test1( true );
   test2( true );
   //test3();

   return 0;

}
