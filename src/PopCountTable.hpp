/**
    REAL: An efficient REad ALigner for next generation sequencing reads.
    Copyright (C) 2010 German Tischler, Solon Pissis

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

#if ! defined(POPCOUNTTABLE_HPP)
#define POPCOUNTTABLE_HPP

#include "types.hpp"
#include "popcnt.hpp"

namespace toollib
{
	/**
	 * Lookup table for population count in 16 bit numbers. The table uses 2^19 bytes.
	 **/
	struct PopCountTable
	{
		private:
		u_int8_t const * const table;
		static u_int8_t * generateTable();
		
		public:
		/**
		 * compute population count of the i+1 most significant bits in m
		 * @param m
		 * @param i
		 * @return population count
		 **/
		unsigned int operator()(unsigned short const m) const { return table[m]; }

		/**
		 * constructor
		 **/
		PopCountTable();
		/**
		 * destructor
		 **/
		~PopCountTable();
	};
	
	struct PopCountBase
	{
		static PopCountTable const pctable;	
	};

	template<bool sse4>	
	struct PopCount : public PopCountBase
	{

		/**
		 * compute population count (number of 1 bits) of 32 bit number val
		 * @param val
		 * @return population count of val
		 **/
		static inline u_int32_t popcount4(u_int32_t val)
		{
			if ( sse4 )
				return PopCnt4<sizeof(unsigned int)>::popcnt4(val);
			else
				return pctable(val >> 16) + pctable(val & 0xFFFF);
		}

		/**
		 * compute population count (number of 1 bits) of the i+1 most significant bits of the 32 bit number val
		 * @param val
		 * @param i index of last bit to consider (starting from MSB)
		 * @return population count
		 **/
		static inline u_int32_t popcount4(u_int32_t val, unsigned int const i)
		{
			return popcount4(val >> (31-i));
		}

		/**
		 * compute population count (number of 1 bits)
		 * @param val
		 * @return population count
		 **/
		static inline u_int64_t popcount8(u_int64_t val)
		{
			if ( sse4 )
				return PopCnt8<sizeof(unsigned long)>::popcnt8(val);
			else
				return popcount4( val >> 32 ) + popcount4( val & 0xFFFFFFFFu);
		}

		/**
		 * compute population count (number of 1 bits) of the i+1 most significant bits of the 64 bit number val
		 * @param val
		 * @param i index of last bit to consider (starting from MSB)
		 * @return population count
		 **/
		static inline u_int64_t popcount8(u_int64_t val, unsigned int const i)
		{
			return popcount8(val >> (63-i));
		}


		static inline unsigned int bitcountpair(u_int32_t const n)
		{
			return popcount8(((n >> 1) | n) & 0x55555555ul);
		}

		static inline unsigned int bitcountpair(u_int64_t const n)
		{
			return popcount8(((n >> 1) | n) & 0x5555555555555555ull);
		}

		static inline unsigned int diffcountpair(u_int32_t const a, u_int32_t const b)
		{
			return bitcountpair( a ^ b );
		}

		static inline unsigned int diffcountpair(u_int64_t const a, u_int64_t const b)
		{
			return bitcountpair( a ^ b );
		}
	};
}
#endif
