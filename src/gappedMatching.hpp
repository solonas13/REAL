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

#if ! defined(GAPPEDMATCHINHG_HPP)
#define GAPPEDMATCHING_HPP

#include <iostream>
#include "AutoTextArray.hpp"
#include "acgtnMap.hpp"

inline unsigned int delta(char const a, char const b)
{
        if ( a==b )
                return 0;
        else
                return 1;
}

struct MatrixAccessor
{
        unsigned int * const M;
        unsigned int const columns;
        
        MatrixAccessor(unsigned int * const rM, unsigned int const rcolumns)
        : M(rM), columns(rcolumns) {}
        
        inline unsigned int       & operator()(unsigned int i, unsigned int j)       { return M[i*columns+j]; }
        inline unsigned int const & operator()(unsigned int i, unsigned int j) const { return M[i*columns+j]; }
};

template<bool sse4, typename iterator>
void gappedMatch(
        unsigned int const patl,
        unsigned int const rpos,
        bool const inverted,
        AutoTextArray<sse4> const & ATA,
	// gaps
	unsigned int * const PG,
	unsigned int * const PH,
	unsigned int const gapwindow,
	unsigned int const restlength,
	iterator const & pattern
        )
{
        // length of seed
        unsigned int const seedl = patl - restlength;

        #if 0        
        std::cerr << "seed position: " << rpos << std::endl;
        #endif
        
        // place where we start scanning the text
        u_int32_t const pos = inverted ? rpos : (rpos + seedl );
        
        #if 0
        std::cerr << "Scan position: " << pos << std::endl;
        #endif

        unsigned int const scanlength =
                inverted
                ?
                std::min ( gapwindow, pos )
                :
                std::min ( gapwindow, static_cast<u_int32_t>(ATA.getN())-pos )
        ;
        
        #if 0
        std::cerr << "scan length " << scanlength << std::endl;
        #endif
        
        if ( 
                (scanlength >= restlength)
                &&
                (inverted ? 
                        ATA.isDontCareFree(pos-scanlength,scanlength) 
                        :
                        ATA.isDontCareFree(pos,scanlength))
        )
        {
                MatrixAccessor G(PG, restlength+1);
                MatrixAccessor H(PH, restlength+1);
        
                for ( unsigned int i = 0; i < scanlength+1; ++i )
                        for ( unsigned int j = 0; j < restlength+1; ++j )
                                G(i,j) = H(i,j) = 0;
                                
                for( unsigned int i = 0; i < scanlength+1; i++)
                        for( unsigned int j = 0; j < restlength+1; j++)
                        {
                                H(i,0) = i;
                                H(0,j) = j;
                        }

                /*DP algorithm*/
                if ( inverted )
                {
                        for(unsigned int i = 1; i < scanlength+1; i++)				//Matrix G and H calculation 
                                for(unsigned int j = 1; j < restlength+1; j++)
                                {                                                        
                                        unsigned int const del = delta( ATA[pos - i], pattern[seedl + j - 1] );
                                        
                                        G(i,j) = G (i-1,j-1) + del;
                                        H(i,j) = 0;

                                        if ( del )
                                        {
                                                unsigned int const rep = G(i,j);
                                                unsigned int const gap = (i<j) ? G(i,i) : G(j,j);
                                                
                                                if(gap < rep)
                                                        H(i,j) = std::abs (static_cast<int>(i) - static_cast<int>(j));

                                                G(i,j) = std::min(rep, gap);
                                        }
                                }
                }
                else
                {
                        for(unsigned int i = 1; i < scanlength+1; i++)				//Matrix G and H calculation 
                                for(unsigned int j = 1; j < restlength+1; j++)
                                {                                                        
                                        unsigned int const del = delta( ATA[pos + i - 1], pattern[seedl + j - 1] );
                                        
                                        G(i,j) = G (i-1,j-1) + del;
                                        H(i,j) = 0;

                                        if ( del )
                                        {
                                                unsigned int const rep = G(i,j);
                                                unsigned int const gap = (i<j) ? G(i,i) : G(j,j);
                                                
                                                if(gap < rep)
                                                        H(i,j) = std::abs (static_cast<int>(i) - static_cast<int>(j));

                                                G(i,j) = std::min(rep, gap);
                                        }
                                }
                }

                // these should be parameters
                unsigned int const maxg = 8;
                unsigned int const maxmis = 4;

                /*figure out where to start from*/
                unsigned int start = 0;
                unsigned int mini = maxmis;		//give priority to least mismatches
                unsigned int min_g = maxg;
                unsigned int const jj = restlength;

                for (unsigned int i = 0; i < scanlength + 1; i++)
                {
                        if ( i < restlength )
                        {
                                if ( G(i,jj) <= maxmis && restlength - i <= maxg )
                                {
                                        if ( G(i,jj) <= mini && restlength - i <= min_g )
                                        {
                                                mini = G(i,jj); 
                                                min_g = restlength - i;
                                                start = i;
                                        }
                                }
                        }
                        else if ( i > restlength )
                        {
                                if ( G(i,jj) <= maxmis && i - restlength <= maxg )
                                {
                                        if ( G(i,jj) <= mini && i - restlength <= min_g )
                                        {
                                                mini = G(i,jj);
                                                min_g = i - restlength; 
                                                start = i;
                                        }
                                }
                        }
                        else if ( i == restlength )
                        {
                                if ( G(i,jj) <= mini ) 
                                {
                                        mini = G(i,jj); 
                                        min_g = restlength - i;	//gap = 0
                                        start = i;
                                }
                        }
                }

                /*Trace back*/
                int jjj = restlength;
                unsigned int gap_length = 0;
                unsigned int gap_pos = 0;
                unsigned int total_mis = G(start,jjj);
                int where = 0;

                for ( int i = start ; i >= 0; )
                {
                        if ( H(i,jjj) == 0 )
                        {
                                i--; jjj--;
                        }
                        else			//here must go only once
                        {
                                if ( i > jjj )	
                                {
                                        where = 1;	//gap is in the text
                                        gap_length = H(i,jjj);
                                        gap_pos = jjj;	
                                        i -= H(i,jjj);
                                }
                                else		
                                {
                                        where = -1;	//gap is in the pattern
                                        gap_length = H(i,jjj);
                                        gap_pos = i;	
                                        jjj -= H(i,jjj);
                                }
                        }
                }
                
                #if 0
                std::cerr << "Pattern ";
                for ( unsigned int i = 0; i < restlength; ++i )
                        std::cerr << toollib::remapChar(pattern[ seedl + i ]);
                std::cerr << std::endl;
                
                std::cerr << "Text ";
                for ( unsigned int i = 0; i < scanlength; ++i )
                        std::cerr << toollib::remapChar ( ATA[pos + (inverted?-(i+1):i)] );
                std::cerr << std::endl;
                
                std::cerr << "H:" << std::endl;
                for ( unsigned int i = 0; i < scanlength+1; ++i )
                {
                        for ( unsigned int j = 0; j < restlength+1; ++j )
                                std::cerr << H(i,j) << " ";
                        std::cerr << std::endl;
                }        
                std::cerr << "G:" << std::endl;
                for ( unsigned int i = 0; i < scanlength+1; ++i )
                {
                        for ( unsigned int j = 0; j < restlength+1; ++j )
                                std::cerr << G(i,j) << " ";
                        std::cerr << std::endl;
                }   
                
                std::cerr << "start " << start << std::endl; 
                #endif    
                
                std::cerr << "Number of mismatches " << total_mis << std::endl;
                
                if ( gap_length )
                {
                        std::cerr << "Length of gap " << gap_length << std::endl;
                        if ( where == 1 )
                                std::cerr << "Gap in pattern at position " << gap_pos << std::endl;
                        else
                                std::cerr << "Gap in text at position " << gap_pos << std::endl;
                }
                else
                {
                        std::cerr << "gap length is zero." << std::endl;
                }
        }
}
#endif
