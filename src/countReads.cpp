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

#include <stdexcept>
#include <cassert>
#include "countReads.hpp"
#include "FastReaderBase.hpp"

u_int64_t countLength(
    std::string const & textfilename,
    std::vector < std::pair < std::string, u_int64_t > > & ranges
    )
{
    FastReaderBase frb(textfilename);
    
    int c;
    bool ignorerestofline = false;
    size_t cnt = 0;
    
    std::string id;
    size_t idcnt = 0;
    
    while ( (c=frb.getNextCharacter()) != -1 )
    {
        switch ( c )
        {
            case '>':
            {
                ignorerestofline = true;
                idcnt = cnt;
                id.resize(0);
                break;
            }
            case '\n':
            {
                if ( ignorerestofline )
                {
                    std::cerr << "idcnt=" << idcnt << " id=" << id << std::endl;    
                    ranges.push_back ( std::pair < std::string, u_int64_t > (id, idcnt) );
                }
                ignorerestofline = false;
                break;
            }
            default:
            {
                if ( ! ignorerestofline )
                {
                    switch ( c )
                    {
                        case 'A': case 'C': case 'G': case 'T': case 'N': cnt++; break;
                    }
                }
                else
                {
                    id += static_cast<char>(c);
                }
            }
            break;
        }
    }
    
    ranges.push_back ( std::pair < std::string, u_int64_t > ("terminal", cnt) );
    
    return cnt;
}

AutoArray<u_int8_t> readFile(std::string const & textfilename, size_t const s)
{
        AutoArray<u_int8_t> AA(s,false);
        u_int8_t * A = AA.get();

    FastReaderBase frb(textfilename);
    
    int c;
    bool ignorerestofline = false;
    
    while ( (c=frb.getNextCharacter()) != -1 )
    {
        switch ( c )
        {
            case '>':
                ignorerestofline = true;
                break;
            case '\n':
                ignorerestofline = false;
                break;
            default:
            {
                if ( ! ignorerestofline )
                {
                    switch ( c )
                    {
                        case 'A': *(A++) = 0; break;
                        case 'C': *(A++) = 1; break;
                        case 'G': *(A++) = 2; break;
                        case 'T': *(A++) = 3; break;
                        case 'N': *(A++) = 4; break;
                    }
                }            
            }
            break;
        }
    }

    return AA;
}
