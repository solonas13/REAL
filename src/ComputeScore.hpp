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

#if ! defined(COMPUTESCORE_HPP)
#define COMPUTESCORE_HPP

#include "AutoTextArray.hpp"
#include "Scoring.hpp"

template<bool sse4, typename pattern_type, bool scores>
struct ComputeScore
{

};

template<bool sse4, typename pattern_type>
struct ComputeScore<sse4,pattern_type,false>
{
        static float computeScore(
                bool const,
                AutoTextArray<sse4> const &,
                pattern_type const &,
                Scoring const &,
                unsigned int const,
                unsigned int const
                )
        {
                return 1.0f;
        }
};

template<bool sse4, typename pattern_type>
struct ComputeScore<sse4,pattern_type,true>
{
        static float computeScore(
                bool const inverted,
                AutoTextArray<sse4> const & ATA,
                pattern_type const & pattern,
                Scoring const & scoring,
                unsigned int const pos,
                unsigned int const patl
                )
        {
                if ( inverted )
                {
                        double rawscore = 1.0f;
                        char const * mapped = inverted ? pattern.transposed : pattern.mapped;
                        unsigned int j = patl;

                        unsigned int wordpos = pos/32;
                        u_int64_t const firstword = ATA.getTextWord( wordpos );
                        u_int64_t const mask = 0x3ull;
                        unsigned int const firstrest = 32-(pos % 32);
                        unsigned int firstloopstop = std::min(firstrest,patl);

                        unsigned int shift = 2*(firstrest-1);
                        unsigned int i = 0;
                        unsigned int rest = patl;
                        
                        for ( ; i < firstloopstop; ++i, shift-=2 )
                        {
                                unsigned int const refbase = ((firstword>>shift)&mask);
                                unsigned int const patbase = mapped[i];
                                unsigned int const quality = pattern.getQuality(--j);
                                double const localrawscore = scoring.getRawLogScoreTable(refbase,patbase,quality);
                                rawscore += localrawscore;
                                // assert ( ((firstword>>shift)&mask) == ATA[pos+i] );
                        }

                        rest -= firstloopstop;
                        
                        while ( rest >= 32 )
                        {
                                shift = 62;
                                unsigned int const loopstop = i+32;
                                u_int64_t const word = ATA.getTextWord(++wordpos);
                                
                                for ( ; i < loopstop; ++i, shift -= 2 )
                                {
                                        unsigned int const refbase = ((word>>shift)&mask);
                                        unsigned int const patbase = mapped[i];
                                        unsigned int const quality = pattern.getQuality(--j);
                                        double const localrawscore = scoring.getRawLogScoreTable(refbase,patbase,quality);
                                        rawscore += localrawscore;
                                        // assert ( ((word>>shift)&mask) == ATA[pos+i] );
                                }
                                
                                rest -= 32;
                        }
                        
                        if ( rest )
                        {
				shift = 62;
				unsigned int const loopstop = i+rest;
				u_int64_t const word = ATA.getTextWord(++wordpos);
					
				for ( ; i < loopstop; ++i, shift -= 2 )
				{
					unsigned int const refbase = ((word>>shift)&mask);
					unsigned int const patbase = mapped[i];
					unsigned int const quality = pattern.getQuality(--j);
					double const localrawscore = scoring.getRawLogScoreTable(refbase,patbase,quality);
					rawscore += localrawscore;
					// assert ( ((word>>shift)&mask) == ATA[pos+i] );
				}
                        }                        

                        return rawscore;
                }
                else
                {
                        double rawscore = 1.0f;
                        char const * mapped = inverted ? pattern.transposed : pattern.mapped;

                        unsigned int wordpos = pos/32;
                        u_int64_t const firstword = ATA.getTextWord( wordpos );
                        u_int64_t const mask = 0x3ull;
                        unsigned int const firstrest = 32-(pos % 32);
                        unsigned int firstloopstop = std::min(firstrest,patl);

                        unsigned int shift = 2*(firstrest-1);
                        unsigned int i = 0;
                        unsigned int rest = patl;
                        
                        for ( ; i < firstloopstop; ++i, shift-=2 )
                        {
                                unsigned int const refbase = ((firstword>>shift)&mask);
                                unsigned int const patbase = mapped[i];
                                unsigned int const quality = pattern.getQuality(i);
                                double const localrawscore = scoring.getRawLogScoreTable(refbase,patbase,quality);
                                rawscore += localrawscore;
                                // assert ( ((firstword>>shift)&mask) == ATA[pos+i] );
                        }

                        rest -= firstloopstop;
                        
                        while ( rest >= 32 )
                        {
                                shift = 62;
                                unsigned int const loopstop = i+32;
                                u_int64_t const word = ATA.getTextWord(++wordpos);
                                
                                for ( ; i < loopstop; ++i, shift -= 2 )
                                {
                                        unsigned int const refbase = ((word>>shift)&mask);
                                        unsigned int const patbase = mapped[i];
                                        unsigned int const quality = pattern.getQuality(i);
                                        double const localrawscore = scoring.getRawLogScoreTable(refbase,patbase,quality);
                                        rawscore += localrawscore;
                                        // assert ( ((word>>shift)&mask) == ATA[pos+i] );
                                }
                                
                                rest -= 32;
                        }
                        
                        if ( rest )
                        {                        
                                shift = 62;
                                unsigned int const loopstop = i+rest;
                                u_int64_t const word = ATA.getTextWord(++wordpos);
                                        
                                for ( ; i < loopstop; ++i, shift -= 2 )
                                {
                                        unsigned int const refbase = ((word>>shift)&mask);
                                        unsigned int const patbase = mapped[i];
                                        unsigned int const quality = pattern.getQuality(i);
                                        double const localrawscore = scoring.getRawLogScoreTable(refbase,patbase,quality);
                                        rawscore += localrawscore;
                                        // assert ( ((word>>shift)&mask) == ATA[pos+i] );
                                }
                        }                        

                        return rawscore;
                }
        }
};
#endif
