/**
    REAL: An efficient REad ALigner for next generation sequencing reads.
    Copyright (C) 2011 German Tischler

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/


#if ! defined(POPCNT_HPP)
#define POPCNT_HPP

#include <types.hpp>

template<unsigned int intsize>
struct PopCnt4
{
	/**
	 * Population count (number of 1 bits)
	 * imported from http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
	 * @param x number
	 * @return number of 1 bits in x
	 **/
	static unsigned int popcnt4(u_int32_t v)
	{
		v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
		v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
		return ( ( (v + (v >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24; // count 
	}
};
template<unsigned int longsize>
struct PopCnt8
{
	/**
	 * Population count (number of 1 bits)
	 * imported from Simon Gog et. al: SDSL
	 * there refered to be designed by Knuth
	 * @param x number
	 * @return number of 1 bits in x
	 **/
	static unsigned int popcnt8(u_int64_t x)
	{
		x = x-( (x>>1) & 0x5555555555555555ull);
		x = (x & 0x3333333333333333ull) + ((x >> 2) & 0x3333333333333333ull);
		x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0full;
		return  (0x0101010101010101ull*x >> 56);
	}
};

#if defined(__GNUC__)
template<>
struct PopCnt4<2>
{
	/**
	 * Population count (number of 1 bits)
	 * using the gnu compilers builtin function for unsigned long types
	 * @param u number
	 * @return number of 1 bits in u
	 **/
	static unsigned int popcnt4(u_int32_t const u)
	{
		return __builtin_popcountl(u);
	}	
};
template<>
struct PopCnt4<4>
{
	/**
	 * Population count (number of 1 bits)
	 * using the gnu compilers builtin function for unsigned long types
	 * @param u number
	 * @return number of 1 bits in u
	 **/
	static unsigned int popcnt4(u_int32_t const u)
	{
		return __builtin_popcount(u);
	}	
};

template<>
struct PopCnt8<4>
{
	/**
	 * Population count (number of 1 bits)
	 * using the gnu compilers builtin function for unsigned long types
	 * @param u number
	 * @return number of 1 bits in u
	 **/
	static unsigned int popcnt8(u_int64_t const u)
	{
		return (__builtin_popcountl(u>>32)) + (__builtin_popcountl(u&0x00000000FFFFFFFFULL));			
	}	
};
template<>
struct PopCnt8<8>
{
	/**
	 * Population count (number of 1 bits)
	 * using the gnu compilers builtin function for unsigned long types
	 * @param u number
	 * @return number of 1 bits in u
	 **/
	static unsigned int popcnt8(u_int64_t const u)
	{
		return __builtin_popcountl(u);
	}	
};
#endif

#endif
