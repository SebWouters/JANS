/*
   JANS: Just Another Number Sieve
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

void test1(){

   std::string rsa_100_n_bin = "101100100011010101100110101111010001111100100000011010101100110111001001011011010001110010101111100100000101111110001110111111011110101011100001010100001110011010111101110010011011101101001111011111110111110110011001001000100111010001010101011101110000001011011101110001110001111010010100001110111101111100010111100101100011111011";
   std::string rsa_100_p_bin = "110011111101111010100000111010110100110101010001111011000011000000000100110101001011001111101101100110011011110011100011000111100110101010000000111110010010011110111";
   std::string rsa_100_q_bin = "110110110111100010100000111111001100011101110101101100001100110111100000000110110000000100010000000010110000010100101111101110101111010100000011111001101111100011101";

   std::string rsa_100_n_dec = "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139";
   std::string rsa_100_p_dec = "37975227936943673922808872755445627854565536638199";
   std::string rsa_100_q_dec = "40094690950920881030683735292761468389214899724061";

   std::string rsa_100_n_hex = "2c8d59af47c81ab3725b472be417e3bf7ab85439af726ed3dfdf66489d155dc0b771c7a50ef7c5e58fb";
   std::string rsa_100_p_hex = "19fbd41d69aa3d86009a967db3379c63cd501f24f7";
   std::string rsa_100_q_hex = "1b6f141f98eeb619bc0360220160a5f75ea07cdf1d";

   jans::big_int n1; n1.set( rsa_100_n_bin, 2  );
   jans::big_int n2; n2.set( rsa_100_n_dec, 10 );
   jans::big_int n3; n3.set( rsa_100_n_hex, 16 );

   const bool equal12 = jans::big_int::compare( n1, n2 );
   const bool equal13 = jans::big_int::compare( n1, n3 );

   std::string n_bin = n3.str( 2  ); const bool equal_bin = ( n_bin.compare( rsa_100_n_bin ) == 0 );
   std::string n_hex = n2.str( 16 ); const bool equal_hex = ( n_hex.compare( rsa_100_n_hex ) == 0 );
   std::cout << "Equal(n1,n2) = " << equal12 << std::endl;
   std::cout << "Equal(n1,n3) = " << equal13 << std::endl;
   std::cout << "Equal(bin)   = " << equal_bin << std::endl;
   std::cout << "Equal(hex)   = " << equal_hex << std::endl;

}

int main()
{

   jans::big_int::sanity_check();
   test1();

   return 0;

}
