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

#include <stdio.h>
#include <iostream>

#include "sieve.h"

void jans::sieve::__solve_gaussian__( unsigned char * out, ubase_t * vectors, const ubase_t d_pow, const ubase_t d_lin ){

   unsigned char * lin_contributes = new unsigned char[ d_lin ];
   unsigned char * pow_contributes = new unsigned char[ d_pow ];
   for ( ubase_t lin = 0; lin < d_lin; lin++ ){ lin_contributes[ lin ] = 1; }
   for ( ubase_t pow = 0; pow < d_pow; pow++ ){ pow_contributes[ pow ] = 1; }

   ubase_t vecspace = d_lin;
   for ( ubase_t pow = 0; pow < d_pow; pow++ ){
      ubase_t num_odds = 0;
      ubase_t last = d_lin;
      for ( ubase_t lin = 0; lin < d_lin; lin++ ){
         if ( ( vectors[ pow + d_pow * lin ] % 2 ) == 1 ){
            num_odds++;
            last = lin;
            if ( num_odds == 2 ){ lin = d_lin; }
         }
      }
      if ( num_odds == 1 ){
         vecspace--;
         lin_contributes[ last ] = 0;
      }
   }

   std::cout << "Removed " << d_lin - vecspace << " of the " << d_lin << " vectors with a unique odd prime power." << std::endl;

   ubase_t redspace = d_pow;
   for ( ubase_t pow = 0; pow < d_pow; pow++ ){
      ubase_t num_odds = 0;
      ubase_t lin = 0;
      for ( ubase_t vec = 0; vec < vecspace; vec++ ){
         while ( lin_contributes[ lin ] == 0 ){ lin++; }
         if ( ( vectors[ pow + d_pow * lin ] % 2 ) == 1 ){
            num_odds++;
            vec = vecspace;
         }
         lin++;
      }
      if ( num_odds == 0 ){
         redspace--;
         pow_contributes[ pow ] = 0;
      }
   }

   std::cout << "Removed " << d_pow - redspace << " of the " << d_pow << " primes with only even powers." << std::endl;

   const ubase_t d_vec = vecspace;
   const ubase_t d_red = redspace;
   const ubase_t d_row = d_red + d_lin;
   const ubase_t d_sol = d_lin - d_pow;

   unsigned char * matrix = new unsigned char[ d_row * d_vec ];

   /*
    * matrix = [  d_red x d_vec  ]
    *          [  -------------  ]
    *          [  d_lin x d_vec  ]
    */

   for ( ubase_t vec = 0; vec < d_vec; vec++ ){
      for ( ubase_t row = 0; row < d_row; row++ ){
         matrix[ row + d_row * vec ] = 0;
      }
   }

   {
      ubase_t lin = 0;
      for ( ubase_t vec = 0; vec < d_vec; vec++ ){
         while ( lin_contributes[ lin ] == 0 ){ lin++; }
         ubase_t pow = 0;
         for ( ubase_t red = 0; red < d_red; red++ ){
            while ( pow_contributes[ pow ] == 0 ){ pow++; }
            matrix[ red + d_row * vec ] = ( vectors[ pow + d_pow * lin ] % 2 );
            pow++;
         }
         matrix[ d_red + lin + d_row * vec ] = 1;
         lin++;
      }
   }

   delete [] lin_contributes;
   delete [] pow_contributes;

   ubase_t start = 0;
   for ( ubase_t red = 0; red < d_red; red++ ){
      bool found   = false;
      ubase_t iter = start;
      while ( ( found == false ) && ( iter < d_vec ) ){
         if ( matrix[ red + d_row * iter ] == 1 ){ found = true; }
         else { iter++; }
      }
      if ( found == true ){
         if ( iter != start ){
            for ( ubase_t row = red; row < d_row; row++ ){
               unsigned char swap            = matrix[ row + d_row * iter  ];
               matrix[ row + d_row * iter  ] = matrix[ row + d_row * start ];
               matrix[ row + d_row * start ] = swap;
            }
         }
         for ( ubase_t vec = iter + 1; vec < d_vec; vec++ ){
            if ( matrix[ red + d_row * vec ] == 1 ){
               for ( ubase_t row = red; row < d_row; row++ ){
                  matrix[ row + d_row * vec ] = ( matrix[ row + d_row * vec ] ) ^ ( matrix[ row + d_row * start ] );
               }
            }
         }
         start++;
      }
   }

   for ( ubase_t sol = 0; sol < d_sol; sol++ ){
      for ( ubase_t row = 0; row < d_lin; row++ ){
         out[ row + d_lin * sol ] = matrix[ d_red + row + d_row * ( d_vec - d_sol + sol ) ];
      }
   }

   delete [] matrix;

}

