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

#if ! defined(GETTEXT_HPP)
#define GETTEXT_HPP

#include <memory>
#include <string>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "AutoTextArray.hpp"
#include "countReads.hpp"

template<bool sse4>
std::auto_ptr < AutoTextArray<sse4> > getText(std::string const & textname, std::vector < std::pair < std::string, u_int64_t > > & ranges)
{
        std::cerr << "Computing length of file " << textname << "...";
        double befl = clock();
        size_t const flen = countLength(textname.c_str(),ranges);
        double aftl = clock();
        std::cerr << "done." << " clocks " << (aftl-befl)/CLOCKS_PER_SEC << " length " << flen << std::endl;

        std::cerr << "Reading file " << textname << "...";
        befl = clock();
        AutoArray<u_int8_t> cont = readFile(textname.c_str(), flen);
        aftl = clock();
        std::cerr << "done." << " clocks " << (aftl-befl)/CLOCKS_PER_SEC << std::endl;
        
        std::auto_ptr < AutoTextArray<sse4> > ATA ( new AutoTextArray<sse4>(cont.get(), flen));
        
        #if 0
        for ( size_t zzz = 0; zzz < flen; ++zzz )
                assert ( (*ATA)[zzz] == cont[zzz] );
        #endif
        
        cont.release();

        return ATA;
}
#endif
