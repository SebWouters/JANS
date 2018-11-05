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

void jans::big_int::__clear__( ubase_t * a ){

   for ( int i = 0; i < NUM_BLOCK; i++ ){ a[ i ] = 0; }

}

void jans::big_int::__copy__( ubase_t * r, ubase_t * a ){

   for ( int i = 0; i < NUM_BLOCK; i++ ){ r[ i ] = a[ i ]; }

}

int jans::big_int::__compare__( ubase_t * a, ubase_t * b ){

   for ( int i = NUM_BLOCK - 1; i >= 0; i-- ){
      if ( a[ i ] > b[ i ] ){ return (  i + 1 ); }
      if ( a[ i ] < b[ i ] ){ return ( -i - 1 ); }
   }
   return 0;

}

void jans::big_int::__shift_up__( ubase_t * a, const int k ){

   const int blk = k / BLOCK_BIT;
   const int bit = k % BLOCK_BIT;

   if ( blk > 0 ){
      for ( int i = NUM_BLOCK - 1; i >= blk; i-- ){ a[ i ] = a[ i - blk ]; }
      for ( int i = blk - 1;       i >= 0;   i-- ){ a[ i ] = 0;            }
   }

   if ( bit > 0 ){
      for ( int i = NUM_BLOCK - 1; i > blk; i-- )  //  del bit               set bit       in case bit is present downward
      {  for ( int j = BLOCK_BIT - 1; j >= bit; j-- ){ a[  i  ] &= ~( 1U << j ); a[  i  ] |= ( ( ( a[   i   ] >> (             j - bit ) ) & 1U ) << j ); }
         for ( int j =       bit - 1; j >= 0;   j-- ){ a[  i  ] &= ~( 1U << j ); a[  i  ] |= ( ( ( a[ i - 1 ] >> ( BLOCK_BIT + j - bit ) ) & 1U ) << j ); }
      }{ for ( int j = BLOCK_BIT - 1; j >= bit; j-- ){ a[ blk ] &= ~( 1U << j ); a[ blk ] |= ( ( ( a[  blk  ] >> (             j - bit ) ) & 1U ) << j ); }
         for ( int j =       bit - 1; j >= 0;   j-- ){ a[ blk ] &= ~( 1U << j ); }
      }
   }

}

void jans::big_int::__shift_down__( ubase_t * a, const int k ){

   const int blk = k / BLOCK_BIT;
   const int bit = k % BLOCK_BIT;
   const int upi = NUM_BLOCK - blk - 1;
   const int upj = BLOCK_BIT - bit;

   if ( blk > 0 ){
      for ( int i = 0;       i <= upi;      i++ ){ a[ i ] = a[ i + blk ]; }
      for ( int i = upi + 1; i < NUM_BLOCK; i++ ){ a[ i ] = 0;            }
   }

   if ( bit > 0 ){
      for ( int i = 0; i < upi; i++ ) //          del bit                   set bit       in case bit is present upward
      {  for ( int j = 0;   j < upj;       j++ ){ a[  i  ] &= ~( 1U << j ); a[  i  ] |= ( ( ( a[   i   ] >> ( j + bit ) ) & 1U ) << j ); }
         for ( int j = upj; j < BLOCK_BIT; j++ ){ a[  i  ] &= ~( 1U << j ); a[  i  ] |= ( ( ( a[ i + 1 ] >> ( j - upj ) ) & 1U ) << j ); }
      }{ for ( int j = 0;   j < upj;       j++ ){ a[ upi ] &= ~( 1U << j ); a[ upi ] |= ( ( ( a[  upi  ] >> ( j + bit ) ) & 1U ) << j ); }
         for ( int j = upj; j < BLOCK_BIT; j++ ){ a[ upi ] &= ~( 1U << j ); }
      }
   }

}

