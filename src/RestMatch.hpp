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

#if ! defined(RESTMATCH_HPP)
#define RESTMATCH_HPP

#include <sys/types.h>

#include "PopCountTable.hpp"

template<bool sse4>
struct RestMatch
{
	static unsigned int getRestLen(unsigned int const readlen, unsigned int const patlen)
	{
                return readlen - patlen;	
	}
	
	static unsigned int getRestWords(unsigned int const readlen, unsigned int const patlen)
	{
                return (getRestLen(readlen,patlen) + 31)/32;
	}
	
	template<typename text_iterator>
	static unsigned int computeDistance(
		u_int64_t const * const restwordarray,
		unsigned int const fullrestwords,
		unsigned int const fracrestsyms,
		text_iterator const & T,
		unsigned int o
	)
	{
		unsigned int dist = 0;
		
		// T.isDontCareFree(fullrestwords*32 + fracrestwords);
		
		for ( unsigned int i = 0; i < fullrestwords; ++i, o += 32 )
		{
			u_int64_t const textword = T.getTextWord(o,32);
			u_int64_t const patword = restwordarray[i];
			unsigned int const d = toollib::PopCount<sse4>::diffcountpair(patword, textword);
			dist += d;

			#if 0
			std::cerr << "Dist of text=" << maskToString(textword,fracrestsyms)
				<< " to pat=" << maskToString(patword,fracrestsyms)
				<< " is " << d << std::endl;
			#endif
		}
		
		if ( fracrestsyms )
		{
			u_int64_t const textword = T.getTextWord(o,fracrestsyms);
			u_int64_t const patword = restwordarray[fullrestwords];
			unsigned int const d = toollib::PopCount<sse4>::diffcountpair(patword, textword);
			dist += d;
			
			#if 0
			std::cerr << "Dist of text=" << maskToString(textword,fracrestsyms)
				<< " to pat=" << maskToString(patword,fracrestsyms)
				<< " is " << d << std::endl;
			#endif
		}
		
		return dist;
	}

	// match offset: offset at which we match the seed
	static unsigned int getMatchOffset(unsigned int const readlen, unsigned int const patlen, bool const inv)
	{
                unsigned int const restlen = getRestLen(readlen,patlen);
                unsigned int const matchoffset = (inv ? restlen : 0);
                return matchoffset;	
	}
	
	static unsigned int getFullRestWords(unsigned int const readlen, unsigned int const patlen)
	{
                unsigned int const restlen = getRestLen(readlen,patlen);
                unsigned int const fullrestwords = restlen / 32;
                return fullrestwords;	        
	}
	
	static unsigned int getFracRestSyms(unsigned int const readlen, unsigned int const patlen, unsigned int const fullrestwords)
	{
                unsigned int const restlen = getRestLen(readlen,patlen);
                unsigned int const fracrestsyms = restlen - fullrestwords * 32;
                return fracrestsyms;
	}

        // offset between start of seed and start of rest	
	static int getTextRestOffset(unsigned int const readlen, unsigned int const patlen, bool const inv)
	{
                unsigned int const restlen = getRestLen(readlen,patlen);
                int const textrestoffset = inv ? (- static_cast<int>(restlen)) : static_cast<int>(patlen);
                return textrestoffset;	
	}

	template<typename pat_iterator>
	static void fillRestWordArray(
		// straight read
	        pat_iterator r, 
	        // pointer to array
	        u_int64_t * const restwordarray, 
	        // reverse and invert
	        bool const inv, 
	        // length of seed
	        unsigned int const seedlen,
	        // length of read
	        unsigned int const readlen,
	        // full rest words
	        unsigned int const fullrestwords,
	        // left over symbols
	        unsigned int const fracrestsyms)
	{
                if ( inv )
                {
                	r += readlen;

			for ( unsigned int i = 0; i < fullrestwords; ++i )
			{
				u_int64_t w = 0;
				
				for ( unsigned int j = 0; j < 32; ++j )
				{
					w <<= 2;
					switch ( *(--r) )
					{
						case 'A': w |= 3; break;
						case 'C': w |= 2; break;
						case 'G': w |= 1; break;
						case 'T': w |= 0; break;
					}
				}
				
				restwordarray[i] = w;
			}
			if ( fracrestsyms )
			{
				u_int64_t w = 0;
				
				for ( unsigned int j = 0; j < fracrestsyms; ++j )
				{
					w <<= 2;
					switch ( *(--r) )
					{
						case 'A': w |= 3; break;
						case 'C': w |= 2; break;
						case 'G': w |= 1; break;
						case 'T': w |= 0; break;
					}
				}
				
				restwordarray[fullrestwords] = w;
			}
                }
                else
                {                
                        r += seedlen;
                
			for ( unsigned int i = 0; i < fullrestwords; ++i )
			{
				u_int64_t w = 0;
				
				for ( unsigned int j = 0; j < 32; ++j )
				{
					w <<= 2;
					switch ( *(r++) )
					{
						case 'A': w |= 0; break;
						case 'C': w |= 1; break;
						case 'G': w |= 2; break;
						case 'T': w |= 3; break;
					}
				}
				
				restwordarray[i] = w;
			}
			if ( fracrestsyms )
			{
				u_int64_t w = 0;
				
				for ( unsigned int j = 0; j < fracrestsyms; ++j )
				{
					w <<= 2;
					switch ( *(r++) )
					{
						case 'A': w |= 0; break;
						case 'C': w |= 1; break;
						case 'G': w |= 2; break;
						case 'T': w |= 3; break;
					}
				}
		
				restwordarray[fullrestwords] = w;
			}
		}
	}

	template<typename pat_iterator>
	static void fillRestWordArrayMapped(
		// straight read
	        pat_iterator r, 
	        // pointer to array
	        u_int64_t * const restwordarray, 
	        // length of seed
	        unsigned int const seedlen,
	        // full rest words
	        unsigned int const fullrestwords,
	        // left over symbols
	        unsigned int const fracrestsyms)
	{
                r += seedlen;
        
                for ( unsigned int i = 0; i < fullrestwords; ++i )
                {
                        u_int64_t w = 0;
                        
                        for ( unsigned int j = 0; j < 32; ++j )
                        {
                                w <<= 2;
                                switch ( *(r++) )
                                {
                                        case 0: w |= 0; break;
                                        case 1: w |= 1; break;
                                        case 2: w |= 2; break;
                                        case 3: w |= 3; break;
                                }
                        }
                        
                        restwordarray[i] = w;
                }
                if ( fracrestsyms )
                {
                        u_int64_t w = 0;
                        
                        for ( unsigned int j = 0; j < fracrestsyms; ++j )
                        {
                                w <<= 2;
                                switch ( *(r++) )
                                {
                                        case 0: w |= 0; break;
                                        case 1: w |= 1; break;
                                        case 2: w |= 2; break;
                                        case 3: w |= 3; break;
                                }
                        }
        
                        restwordarray[fullrestwords] = w;
                }
	}

	template<typename pat_iterator>
	static void fillRestWordArrayReverseMapped(
		// straight read
	        pat_iterator r, 
	        // pointer to array
	        u_int64_t * const restwordarray, 
	        // length of read
	        unsigned int const readlen,
	        // full rest words
	        unsigned int const fullrestwords,
	        // left over symbols
	        unsigned int const fracrestsyms)
	{
                r += readlen;

                for ( unsigned int i = 0; i < fullrestwords; ++i )
                {
                        u_int64_t w = 0;
                        
                        for ( unsigned int j = 0; j < 32; ++j )
                        {
                                w <<= 2;
                                switch ( *(--r) )
                                {
                                        case 0: w |= 3; break;
                                        case 1: w |= 2; break;
                                        case 2: w |= 1; break;
                                        case 3: w |= 0; break;
                                }
                        }
                        
                        restwordarray[i] = w;
                }
                if ( fracrestsyms )
                {
                        u_int64_t w = 0;
                        
                        for ( unsigned int j = 0; j < fracrestsyms; ++j )
                        {
                                w <<= 2;
                                switch ( *(--r) )
                                {
                                        case 0: w |= 3; break;
                                        case 1: w |= 2; break;
                                        case 2: w |= 1; break;
                                        case 3: w |= 0; break;
                                }
                        }
                        
                        restwordarray[fullrestwords] = w;
                }
	}

};
#endif
