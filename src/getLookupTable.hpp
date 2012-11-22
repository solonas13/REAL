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

#if ! defined(GETLOOKUPTABLE_HPP)
#define GETLOOKUPTABLE_HPP

#include <cstdlib>
#include "AutoArray.hpp"

template<typename mask_type>
AutoArray< size_t > getLookupTable(
        mask_type const * Alist0,
        size_t const masks,
        unsigned int const s0shift,
        unsigned int const histsize
        )
{
        AutoArray< size_t > Ahist0(2*histsize,false);
        std::fill(Ahist0.get(), Ahist0.get() + 2 * histsize, 0);

        u_int32_t lasthash = 0;
        size_t lastindex = 0;
        for ( size_t i = 0; i < masks; ++i )
        {
                u_int32_t const hash = (Alist0[i].sign >> s0shift);
                
                if ( hash != lasthash )
                {
                        Ahist0[2*lasthash + 0] = lastindex, Ahist0[2*lasthash + 1] = i;
                        lastindex = i, lasthash = hash;
                }
        }
        Ahist0[2*lasthash + 0] = lastindex, Ahist0[2*lasthash + 1] = masks;
                
        return Ahist0;
}

template<typename projector_type>
AutoArray< u_int32_t > getLookupTableProjected(
        u_int32_t const * const Alist0,
        size_t const masks,
        unsigned int const s0shift,
        unsigned int const histsize,
        projector_type const & projector
        )
{
        AutoArray< u_int32_t > Ahist0(2*histsize,false);
        std::fill(Ahist0.get(), Ahist0.get() + 2 * histsize, 0);

        u_int32_t lasthash = 0;
        size_t lastindex = 0;
        for ( size_t i = 0; i < masks; ++i )
        {
                u_int32_t const hash = (projector(Alist0[i]) >> s0shift);
                
                if ( hash != lasthash )
                {
                        Ahist0[2*lasthash + 0] = lastindex;
                        Ahist0[2*lasthash + 1] = i;
                        lastindex = i;
                        lasthash = hash;
                }
        }
        Ahist0[2*lasthash + 0] = lastindex, Ahist0[2*lasthash + 1] = masks;
                
        return Ahist0;
}
#endif
