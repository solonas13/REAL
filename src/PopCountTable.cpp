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

#if defined(HAVE_CONFIG_H)
#include "real_config.hpp"
#endif

#include "PopCountTable.hpp"
#include <cassert>

namespace toollib
{
	PopCountTable::PopCountTable() : table(generateTable()) {}
	PopCountTable::~PopCountTable() { delete [] table; }	
	u_int8_t * PopCountTable::generateTable()
	{
		u_int8_t * table = 0;
		
		try
		{
			table = new u_int8_t [ 1u<<16 ];
			
			for ( u_int32_t m = 0; m < (1u<<16); ++m )
			{
				u_int8_t c = 0;
				
				for ( u_int32_t t = m; t; t >>= 1 )
					if ( t & 1 )
						++c;
				
				table[m] = c;
			}

			return table;
		}
		catch(...)
		{
			delete [] table;
			table = 0;
			throw;
		}
	}
}
