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

#if ! defined(PATTERN_HPP)
#define PATTERN_HPP

#include <string>
#include <sstream>
#include "acgtnMap.hpp"

struct FastAReader;
struct FastQReader;

struct PatternBase
{
        friend struct FastAReader;
        friend struct FastQReader;
        
	unsigned int patlen;
	size_t patid;
	
	char const * mapped;
	char const * transposed;

	PatternBase() : patlen(0), patid(0), mapped(0), transposed(0) {}

	unsigned int getQuality(unsigned int const /* i */) const
	{
		return 30;
	}
	
	bool isDontCareFree() const
	{
	        for ( unsigned int i = 0; i < patlen; ++i )
	                if ( mapped[i] > 3 )
	                        return false;
                return true;
	}
	
	unsigned int getPatternLength() const
	{
	        return patlen;
	}
	
	size_t getPatID() const
	{
	        return patid;
	}
	inline void computeMapped()
	{}
	std::string getStringId() const
	{
	        std::ostringstream ostr;
	        ostr << patid;
	        ostr.flush();
	        return ostr.str();
	}
};

struct PatternQualityBase : public PatternBase
{
	char const * quality;

	PatternQualityBase() : PatternBase(), quality(0) {}

	unsigned int getQuality(unsigned int const i) const
	{
		return quality[i];
	}	
};

struct Pattern : public PatternBase
{
        friend struct FastAReader;
        friend struct FastQReader;
        
        public:
	std::string sid;

	private:	
	std::string spattern;
	char const * pattern;
	
	public:	
	std::string smapped;
	std::string stransposed;
	
	Pattern() : PatternBase() {}
	
	void computeMapped()
	{
	        smapped.resize(patlen);
	        stransposed.resize(patlen);

	        std::string::iterator       smapped_a = smapped.begin();
	        std::string::iterator const smapped_e = smapped.end();
	        std::string::reverse_iterator stransposed_a = stransposed.rbegin();
	        char const * tpattern = pattern;
	        while ( smapped_a != smapped_e )
	        {
	                switch ( *(tpattern++) )
	                {
	                        case 'A': *(smapped_a++) = toollib::mapChar('A'); *(stransposed_a++) = toollib::invertN(toollib::mapChar('A')); break;
	                        case 'C': *(smapped_a++) = toollib::mapChar('C'); *(stransposed_a++) = toollib::invertN(toollib::mapChar('C')); break;
	                        case 'G': *(smapped_a++) = toollib::mapChar('G'); *(stransposed_a++) = toollib::invertN(toollib::mapChar('G')); break;
	                        case 'T': *(smapped_a++) = toollib::mapChar('T'); *(stransposed_a++) = toollib::invertN(toollib::mapChar('T')); break;
	                        default:  *(smapped_a++) = toollib::mapChar('N'); *(stransposed_a++) = toollib::invertN(toollib::mapChar('N')); break;
	                }
                }
                
                mapped = smapped.c_str();
                transposed = stransposed.c_str();
	}
	std::string const & getStringId() const
	{
	        return sid;
	}
};
#endif
