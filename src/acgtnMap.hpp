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

#if ! defined(ACGTNMAP_HPP)
#define ACGTNMAP_HPP

#include <string>
#include <algorithm>

namespace toollib
{
        inline static char remapChar(char c)
        {
                switch ( c )
                {
                        case 0: return 'A'; break;
                        case 1: return 'C'; break;
                        case 2: return 'G'; break;
                        case 3: return 'T'; break;
                        default: return 'N'; break;
                }
        }

        inline static char mapChar(char c)
        {
                switch ( c )
                {
                        case 'A': return 0; break;
                        case 'C': return 1; break;
                        case 'G': return 2; break;
                        case 'T': return 3; break;
                        default: return 4; break;
                }
        }

        inline std::string remapString(std::string a)
        {
                for ( unsigned int i = 0; i < a.size(); ++i )
                        a[i] = remapChar(a[i]);
                return a;
        }

        inline char invertN(char a)
        {
                switch ( a )
                {
                        case 0: return 3;
                        case 1: return 2;
                        case 2: return 1;
                        case 3: return 0;
                        default: return 4;
                }
        }
}
#endif
