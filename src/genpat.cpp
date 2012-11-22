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

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <ctime>
#include "getText.hpp"
#include "StreamTextFile.hpp"

u_int64_t rand64()
{
	return
		(static_cast<u_int64_t>(rand() & 0xFF) << 0) |
		(static_cast<u_int64_t>(rand() & 0xFF) << 8) |
		(static_cast<u_int64_t>(rand() & 0xFF) << 16) |
		(static_cast<u_int64_t>(rand() & 0xFF) << 24) |
		(static_cast<u_int64_t>(rand() & 0xFF) << 32) |
		(static_cast<u_int64_t>(rand() & 0xFF) << 40) |
		(static_cast<u_int64_t>(rand() & 0xFF) << 48) |
		(static_cast<u_int64_t>(rand() & 0xFF) << 56);		
}

char change(char c)
{
	char tc = c;
	
	while ( tc == c )
	{
		switch ( rand64() % 4 )
		{
			case 0: tc = 'A'; break;
			case 1: tc = 'C'; break;
			case 2: tc = 'G'; break;
			case 3: tc = 'T'; break;
			default: break;
		}
	}
	
	return tc;
}

// filename, numpat, maxerr
int main(int argc, char * argv[])
{
	srand(time(0));

	if ( argc < 6 )
	{
		std::cerr << "usage: " << argv[0] << " <filename> <numpat> <patlen> <errprob> <1|0[FASTQ|FASTA]>" << std::endl;
		return EXIT_FAILURE;
	}
	
	std::string filename = argv[1];
	unsigned int const numpat = atoi(argv[2]);
	unsigned int const patlen = atoi(argv[3]);
	double const errprob = atof(argv[4]);
	bool const fastq = atoi(argv[5]);
		
	bool const sse4 = false;
	std::vector < std::pair < std::string, u_int64_t > > ranges;
	std::auto_ptr < AutoTextArray<sse4> > AATA = getText<sse4>(filename,ranges);
	AutoTextArray<sse4> const & ATA = *(AATA.get());
	
	u_int64_t numpos = (ATA.getN() >= patlen) ? (ATA.getN() - patlen + 1) : 0;
			
	std::cerr << "Number of positions is " << numpos << std::endl;	
	
	if ( ! numpos )
		return EXIT_SUCCESS;
	
	std::vector<u_int64_t> vpos;
	for ( u_int64_t i = 0; i < numpat; ++i )
		vpos.push_back(rand64() % numpos);
		
	std::sort(vpos.begin(),vpos.end());
	
	for ( unsigned int i = 0; i < vpos.size(); ++i )
	{
	        // u_int64_t const word = ATA.getTextWord(vpos[i], patlen);
	        std::string s(patlen,' ');
	        for ( unsigned int j = 0; j < patlen; ++j )
	                switch ( ATA[ vpos[i] + j] )
	                {
	                        case 0: s[j] = 'A'; break;
	                        case 1: s[j] = 'C'; break;
	                        case 2: s[j] = 'G'; break;
	                        case 3: s[j] = 'T'; break;
                                default: s[j] = 'N'; break;
                        }
                        
		bool inv = rand64() & 1;	
		
		if ( inv )
                        SignatureConstruction<unsigned int>::reverseString(s);		        
                
                std::string const r = s;
		
		std::ostringstream commentstr;
			
		commentstr << (fastq?"@":">");
		commentstr << "p" << vpos[i];
		
		if ( inv )
			commentstr << "_inv";
		
		for ( unsigned int j = 0; j < patlen; ++j )
		{
			if ( ((rand64() % 65536) / 65536.0) <= errprob )
			{
				commentstr << "_" << j << s[j];
				char t = change(s[j]);
				s[j] = t;
				commentstr << t;
				// std::cerr << "Modify position " << j << std::endl;
			}
		}
		
		if ( fastq )
			commentstr << " length=" << patlen;
		
		if ( !fastq )
		{
			std::cout << commentstr.str() << std::endl;
			std::cout << s << std::endl;
		}
		else
		{
			std::cout << commentstr.str() << std::endl;
			std::cout << s << std::endl;
			std::cout << "+" << std::endl;
			
			for ( unsigned int j = 0; j < patlen; ++j )
			        if ( s[j] == r[j] )
        				std::cout << "D";
                                else
                                        std::cout << "*";
			std::cout << std::endl;
		}
		
		if ( i % 1000000 == 0 )
			std::cerr << i << std::endl;
	}
		
	return EXIT_SUCCESS;
}
