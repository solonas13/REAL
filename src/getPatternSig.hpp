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

#if ! defined(GETPATTERNSIG_HPP)
#define GETPATTERNSIG_HPP

#include "types.hpp"

#include <stdexcept>

template<unsigned int k>
u_int32_t getSignature(u_int8_t const * p)
{
    u_int32_t m = 0;
    
    switch ( *(p++) )
    {
        case 'A':
            break;
        case 'C':
            m |= 1;
            break;
        case 'G':
            m |= 2;
            break;
        case 'T':
            m |= 3;
            break;
        default:
            m |= 4;
            break;
    }

    for ( unsigned int i = 1; i < k; ++i )
    {
        m <<= 3;

        switch ( *(p++) )
        {
            case 'A':
                break;
            case 'C':
                m |= 1;
                break;
            case 'G':
                m |= 2;
                break;
            case 'T':
                m |= 3;
                break;
            default:
                m |= 4;
                break;
        }    
    }
    
    return m;
}

inline u_int32_t getPartSig(
    u_int8_t const * base,
    unsigned int const partsyms
    )
{
    switch ( partsyms )
    {
        case 0: return 0;
        case 1: return getSignature<1>( base );
        case 2: return getSignature<2>( base );
        case 3: return getSignature<3>( base );
        case 4: return getSignature<4>( base );
        default: throw std::runtime_error("getPartSig called for unhandled length.");
    }
}

inline u_int32_t getPartSig(
    u_int8_t const ** L,
    unsigned int const i,
    unsigned int const fullparts,
    unsigned int const siglen,
    unsigned int const partsyms
    )
{
    switch ( partsyms )
    {
        case 0: return 0;
        case 1: return getSignature<1>( L[ 2*i+1 ] + fullparts * siglen );
        case 2: return getSignature<2>( L[ 2*i+1 ] + fullparts * siglen );
        case 3: return getSignature<3>( L[ 2*i+1 ] + fullparts * siglen );
        case 4: return getSignature<4>( L[ 2*i+1 ] + fullparts * siglen );
        default: throw std::runtime_error("getPartSig called for unhandled length.");
    }
}
#endif
