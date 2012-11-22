#if ! defined(MATCH_HPP)
#define MATCH_HPP

#include <sys/types.h>
#include "RangeVector.hpp"
#include "AutoTextArray.hpp"
#include "RestWordBuffer.hpp"
#include "Scoring.hpp"
#include "ComputeScore.hpp"

#include "UniqueMatchInfo.hpp"

/*The following function gives the total scoring of an alignment in constant time
Note: double current_score is the value of G[i][m], i.e. the score of an alignment without the gap penalties
*/
inline double total_scoring
( 
        unsigned int gap, 
        double current_score, 
        double gap_open_penalty, 
        double gap_extend_penalty, 
        double gap_open_offset_penalty 
)
{
        if ( gap % 3 == 0 )                             //The mod3 has some biological meaning which I am not aware
                return ( current_score + ( gap * gap_extend_penalty ) + gap_open_penalty + gap_open_offset_penalty );    
        else
                return ( current_score + ( gap * gap_extend_penalty ) + gap_open_penalty );      
}


/*The following function computes the limits of the i-th coordinate for the matrix G in constant time*/
inline unsigned int i_limits
( 
        unsigned int n, 
        unsigned int m, 
        unsigned int* up, 
        unsigned int* down, 
        unsigned int MAXgap 
)
{
        if ( (int) m - (int) MAXgap < 0 )       (* up )    = 0;
        else                                    (* up )    = m - MAXgap;
        if ( m + MAXgap > n )                   (* down )  = n;
        else                                    (* down )  = m + MAXgap;
        return 0;
}

/*The following function computes the limits of the j-th coordinate for matrix G and H in constant time*/
inline unsigned int j_limits
( 
        unsigned int i, 
        unsigned int m, 
        unsigned int* left, 
        unsigned int* right, 
        unsigned int MAXgap 
)
{
	if ( (int) i - (int) MAXgap > 0 )	(* left )   = i - MAXgap;
	else					(* left )   = 1;
	if ( i + MAXgap > m ) 			(* right )  = m;
	else					(* right )  = i + MAXgap;
	return 0;
}

template<typename data_type>
struct AgmMatrix
{
        u_int64_t const n;
        u_int64_t const m;
        AutoArray< data_type > AG; 
        AutoArray< data_type * > G;
        
        AgmMatrix(unsigned int const rn, unsigned int const rm)
        : n(rn), m(rm), AG((m+1) * (n+1)), G(n+1)
        {
                for ( u_int64_t i = 0; i < n+1; ++i )
                        G[i] = AG.get() + i*(m+1);
                for ( u_int64_t i = 0; i < (n+1)*(m+1); ++i )
                        AG[i] = 0;
        }
        
        data_type * operator[](unsigned int i)
        {
                return G[i];
        }
        data_type const * operator[](unsigned int i) const
        {
                return G[i];
        }
};

template<typename data_type>
std::ostream & operator<<(std::ostream & out, AgmMatrix<data_type> const & matrix)
{
        out << "<AgmMatrix>(" << std::endl;
        
        for ( u_int64_t i = 0; i < matrix.n+1; ++i )
        {
                for ( u_int64_t j = 0; j < matrix.m+1; ++j )
                {
                         out << matrix[i][j];
                         if ( j+1 < (matrix.m+1) )
                                 out << "\t";       
                }
                std::cerr << std::endl;
        }
        
        out << "</AgmMatrix>)";
        
        return out;
}

/*The following function gives the position of the gap in O(m) time*/
inline unsigned int backtracing 
( 
        AgmMatrix<unsigned int> & H,
        unsigned int n,
        unsigned int m, 
        unsigned int start, 
        unsigned int where, 
        unsigned int* gap_pos 
)
{
        int i, j;
        ( *gap_pos ) = 0;

        if ( where == 1 || where == 2 )
        {
                i = start; j = m;       //we start backtracing from the last column
        }
        else
        {
                i = n; j = start;       //we start backtracing from the last row
        }
        while ( i >= 0 )
        {
                if ( H[i][j] == 0 )
                {
                        --i; --j;
                }
                else                            
                {
                        if ( i > j )    
                                ( *gap_pos ) = j;       
                        else            
                                ( *gap_pos ) = i;
                        break;  
                }
        }
        return 0;
}

/*The following function computes the optimal alignment using matrix G in O(2*MAXgap+1) time
Note:   double MINscore is the minimum score of any alignment used for initialisation.
        double gap_open_penalty, double gap_extend_penalty, double gap_open_offset_penalty are arguments given by the user to represent the gap penalty.
*/

inline unsigned int opt_solution 
( 
        AgmMatrix<double> & G,
        unsigned int n, 
        unsigned int m, 
        unsigned int MAXgap,
        double MINscore, 
        double gap_open_penalty, 
        double gap_extend_penalty, 
        double gap_open_offset_penalty, 
        double* MAXscore, 
        unsigned int* MINgap, 
        unsigned int* where, 
        unsigned int* start 
)
{
        double score = MINscore;                                //Initialisation
                                
        unsigned int up = 0;
        unsigned int down = 0;
        i_limits( n, m, &up, &down, MAXgap );                   // computes the i coordinates for matrix G for the last column

        for ( unsigned int i = up ; i <= down ; i++ )
        {
                double temp_score = 0;
                if ( i < m )
                {
                        if ( G[i][m] >= MINscore && m - i <= MAXgap )
                        {
                                temp_score = total_scoring ( m - i, G[i][m], gap_open_penalty, gap_extend_penalty, gap_open_offset_penalty );
                                if ( temp_score > score )
                                {
                                        score = temp_score;
                                        ( *MAXscore ) = score; 
                                        ( *MINgap ) = m - i;
                                        ( *where ) = 1;         //where: gap is in the text and start backtracing from the last column
                                        ( *start ) = i;         //backtrace from cell G[start,m]
                                }
                        }
                }
                else if ( i > m )
                {
                        if ( G[i][m] >= MINscore && i - m <= MAXgap )
                        {
                                temp_score = total_scoring( i - m, G[i][m], gap_open_penalty, gap_extend_penalty, gap_open_offset_penalty );
                                if (  temp_score > score )
                                {
                                        score = temp_score;
                                        ( *MAXscore ) = score; 
                                        ( *MINgap ) = i - m;
                                        ( *where ) = 2;         //where: gap is in the pattern and start backtracing from last column
                                        ( *start ) = i;         //backtrace from cell G[start,m]
                                }
                        }
                }
                else if ( i == m )
                {
                        if ( G[i][m] >= MINscore )
                        {
                                temp_score = total_scoring( 0, G[i][m], gap_open_penalty, gap_extend_penalty, gap_open_offset_penalty );
                                if (  temp_score > score ) // mgap = 0 
                                {
                                        score = temp_score;
                                        ( *MAXscore ) = score; 
                                        ( *MINgap ) = 0;
                                        ( *where ) = 0;         //there is no gap
                                        ( *start ) = m;         //no need to backtrace
                                }
                        }
                }
        }

        if ( m + MAXgap > n )
        {
                #if 0
                std::cerr << "It should not go here!!!" << std::endl;
                getchar();
                #endif
                unsigned int left = 0;
                unsigned int right = 0;
                j_limits ( n, m, &left, &right, MAXgap );       // computes the j coordinates for matrix G for the last row

                for ( unsigned int j = left ; j < right ; j++ )
                {
                        double temp_score = 0;
                        if ( G[n][j] >= MINscore && n - j <= MAXgap )
                        {
                                temp_score = total_scoring( n - j, G[n][j], gap_open_penalty, gap_extend_penalty, gap_open_offset_penalty );
                                if (  temp_score > score )
                                {
                                        score = temp_score;
                                        ( *MAXscore ) = score; 
                                        ( *MINgap ) = n - j;
                                        ( *where ) = 3;         //where: gap is in the pattern and start backtracing from last row
                                        ( *start ) = j;         //backtrace from cell G[n,start]
                                }
                        }
                }
        }        

        return 0;
}

/*The following function is the algorithm for computing the dynamic-programming matrix G[0..n+1,0..m+1] 
in Theta((2*MAXgap+1)*m) time and Theta(2*n*m)space. 
Note:	unsigned int MAXgap is an argument given by the user to represent the maximum allowed gap
 	The following function also computes matrix H[0..n+1,0..m+1] for the backtracing.
*/
template<typename text_type, typename pattern_type, typename quality_type>
inline unsigned int agm 
( 
        AgmMatrix<double> & G,
        AgmMatrix<unsigned int> & H,
	text_type const & t,
	u_int64_t const ri,
	unsigned int n,
	pattern_type const & pattern, 
	u_int64_t const rj,
	unsigned int m, 
	unsigned int MAXgap,
	quality_type const & quality,
	Scoring const & scoring
)
{
	for( unsigned int i = 1; i < n + 1 ; i++ )		// matrix G and H calculation
	{
		if( i < MAXgap + 1)	H[i][0] = i;		// matrix H initialisation

		unsigned int left = 0;
		unsigned int right = 0;
		j_limits ( i, m, &left, &right, MAXgap );	// computes the j coordinates for matrix G and H	 

		for( unsigned int j = left ; j <= right ; j++ )
		{
			double gap = 0;
			double mis = 0;

			if( j < MAXgap + 1)	H[0][j] = j;	// matrix H initialisation

			if( i < j )
			{
				mis = G[i-1][j-1] + scoring.getRawLogScoreTable ( t[i+ri], pattern[j+rj], quality.getQuality(j+rj) ); 
				gap = G[i][i];
				G[i][j] = std::max ( mis, gap );

				#if 0
				std::cerr << "i=" << i << " j=" << j << " mis=" << mis << " gap=" << gap << std::endl;
				#endif
				
				if ( gap > mis )        H[i][j] = j - i;
				else		        H[i][j] = 0;
			}
			else if ( i > j )
			{
				mis = G[i-1][j-1] + scoring.getRawLogScoreTable ( t[i+ri], pattern[j+rj], quality.getQuality(j+rj) );
				gap = G[j][j];
				G[i][j] = std::max ( mis, gap );

				#if 0
				std::cerr << "i=" << i << " j=" << j << " mis=" << mis << " gap=" << gap << std::endl;
				#endif

				if ( gap > mis )	H[i][j] = i - j;
				else			H[i][j] = 0;
			}
			else if ( i == j )
			{
				G[i][j] = G[i-1][j-1] + scoring.getRawLogScoreTable ( t[i+ri], pattern[j+rj], quality.getQuality(j+rj) );    
                                H[i][j] = 0;
			}
		}
	}	
	return 0;
}


template<bool sse4, typename a_mask_type, typename b_mask_type, typename iterator, typename pattern_type, bool scores, typename updater>
void match(
        a_mask_type const * const list_a,
        b_mask_type const * const list_b,
        unsigned int const seedmaxk,
        unsigned int const totalmaxk,
        typename a_mask_type::signature_type const s_a,
        typename b_mask_type::signature_type const s_b,
        bool const inverted,
        u_int64_t const fileid,
        unsigned int const sampleshift,
        size_t const * const lookup,
        RangeVector<sse4> const & RV,
        /* rest matching */
        AutoTextArray<sse4> const & ATA,
        RestWordBuffer<sse4> const & RWB,
	/* scoring */
        pattern_type const & pattern,
        Scoring const & scoring,
        float const epsilon,
        //
        typename updater::info_type & info
        )
{
        unsigned int matchoffset; 
        int textrestoffset;
        u_int64_t const * restwordarray;
        
        if ( inverted )
        {
                matchoffset = RWB.reversematchoffset;
                textrestoffset = RWB.reversetextrestoffset;
                restwordarray = RWB.Breverse;
        }
        else
        {
                matchoffset = RWB.straightmatchoffset;
                textrestoffset = RWB.straighttextrestoffset;
                restwordarray = RWB.Bstraight;
        }

        unsigned int const prefix = s_a >> sampleshift;
        size_t const low = lookup[2*prefix + 0];
        size_t const high = lookup[2*prefix + 1];

        std::pair< a_mask_type const *, a_mask_type const * > eq = 
                std::equal_range (list_a + low, list_a + high, a_mask_type(s_a));

        for ( a_mask_type const * p = eq.first; p != eq.second; ++p )
        {
        	// compute number of errors in seed
        	unsigned int const seedk = toollib::PopCount<sse4>::diffcountpair(s_b,list_b[p->ptr].sign);
                
                if( seedk <= seedmaxk )
                {
                	// seed position
                        u_int32_t const rpos = p->getPos(list_b);
                        
                        if ( rpos >= matchoffset )
                        {
                        	// suspected pattern position
                                u_int32_t const pos = rpos - matchoffset;
   
                                if ( RV.isPositionValid(pos,RWB.patl) && ATA.isDontCareFree(pos,RWB.patl) )
                                {
                                        u_int32_t const restpos = rpos + textrestoffset;
                                        
                                        unsigned int const restk = RestMatch<sse4>::computeDistance ( restwordarray, RWB.fullrestwords, RWB.fracrestsyms, ATA, restpos );
                                        unsigned int const totalk = seedk + restk;

                                        if ( totalk <= totalmaxk )
                                        {
                                                unsigned int const fragid = RV.positionToRange(pos);
                                                double const score = ComputeScore<sse4,pattern_type,scores>::computeScore(inverted,ATA,pattern,scoring,pos,RWB.patl);
                                                
                                                updater::update(inverted,fileid,pos,totalk,score,epsilon,fragid,info);
                                        }
                                }
                        }
                }
        }        
}

#include <map>

struct GapInfo
{
        unsigned int MINgap;    //the minimum gap achieved 
        unsigned int where;     //where is the gap
        unsigned int start;     //where to start backtracing from
        unsigned int gap_pos;   //the position of the gap
};

template<bool sse4, typename a_mask_type, typename b_mask_type, typename iterator, typename pattern_type, bool scores, typename updater>
void matchGaps(
        a_mask_type const * const list_a,
        b_mask_type const * const list_b,
        unsigned int const seedmaxk,
        typename a_mask_type::signature_type const s_a,
        typename b_mask_type::signature_type const s_b,
        bool const inverted,
        // u_int64_t const fileid,
        unsigned int const sampleshift,
        size_t const * const lookup,
        RangeVector<sse4> const & RV,
        /* rest matching */
        AutoTextArray<sse4> const & ATA,
        RestWordBuffer<sse4> const & RWB,
	/* scoring */
        pattern_type const & pattern,
        Scoring const & scoring,
        float const /* epsilon */,
        //
        typename updater::info_type & info,
        unsigned int patid,
        std::map<unsigned int, GapInfo> & gapinfos
        )
{
        unsigned int matchoffset; 
        int textrestoffset;
        u_int64_t const * restwordarray;
        
        if ( inverted )
        {
                matchoffset = RWB.reversematchoffset;
                textrestoffset = RWB.reversetextrestoffset;
                restwordarray = RWB.Breverse;
        }
        else
        {
                matchoffset = RWB.straightmatchoffset;
                textrestoffset = RWB.straighttextrestoffset;
                restwordarray = RWB.Bstraight;
        }

        unsigned int const prefix = s_a >> sampleshift;
        size_t const low = lookup[2*prefix + 0];
        size_t const high = lookup[2*prefix + 1];

        std::pair< a_mask_type const *, a_mask_type const * > eq = 
                std::equal_range (list_a + low, list_a + high, a_mask_type(s_a));
        
        if ( updater::matchGaps(info) )
        {
                for ( a_mask_type const * p = eq.first; p != eq.second; ++p )
                {
                        // compute number of errors in seed
                        unsigned int const seedk = toollib::PopCount<sse4>::diffcountpair(s_b,list_b[p->ptr].sign);
                        
                        if( seedk <= seedmaxk )
                        {
                                // seed position
                                u_int32_t const rpos = p->getPos(list_b);
                                
                                if ( 
                                        RV.isPositionValid(rpos,RWB.seedl) 
                                        && 
                                        ATA.isDontCareFree(rpos,RWB.seedl) 
                                )
                                {
                                        unsigned int const range = RV.positionToRange(rpos);
                                        // u_int64_t const plow = RV.ranges[range].second;
                                        u_int64_t const phigh = RV.ranges[range+1].second;
                                
                                        if ( ! inverted )
                                        {
                                                double const seedscore = ComputeScore<sse4,pattern_type,scores>::computeScore(inverted,ATA,pattern,scoring,rpos,RWB.seedl);
                                        
                                                // rest of text
                                                u_int64_t const n = std::min ( static_cast<u_int64_t>(phigh - RWB.seedl - rpos), static_cast<u_int64_t>(2 * RWB.patl) );
                                                u_int64_t const m = RWB.patl - RWB.seedl;
                                                
                                                if ( n && m && ATA.isDontCareFree(rpos + RWB.seedl,n) )
                                                {
                                                        #if 0
                                                        std::cerr << "n=" << n << " m=" << m << std::endl;
                                                        #endif
                                                
                                                        AgmMatrix<double> G(n,m);
                                                        AgmMatrix<unsigned int> H(n,m);
                                                        double MAXscore;        //the returned scored
                                                        unsigned int MINgap;    //the minimum gap achieved 
                                                        unsigned int where;     //where is the gap
                                                        unsigned int start;     //where to start backtracing from
                                                        unsigned int gap_pos;   //the position of the gap

                                                        agm(G,H,ATA, rpos + RWB.seedl - 1 , n, pattern.mapped, RWB.seedl-1, RWB.patl-RWB.seedl, 3 , pattern, scoring );
                                                        opt_solution ( G, n, m, 3 /*MAXgap*/, -100.00/*MINscore*/, -1.0/*gap_open_penalty*/, -1.0/*gap_extend_penalty*/,  -1.0/*gap_open_offset_penalty*/,  &MAXscore, &MINgap, &where, &start );
                                                        
                                                        backtracing ( H, n, m, start, where, &gap_pos );

                                                        if ( where == 1 )
                                                        {
                                                                //gap is in the text[gap_pos]
                                                        }                
                                                        else if ( where == 2 || where == 3 )
                                                        {
                                                                //gap is in the pattern[gap_pos]
                                                        }

                                                        #if 0
                                                        std::cerr << G;
                                                        std::cerr << H;
                                                        #endif
                                                        
                                                        if ( MINgap )
                                                        {
                                                                double const completescore = seedscore + MAXscore;
                                                              
                                                                // no previous match  
                                                                if ( info.getState() == UniqueMatchInfoBase::NoMatch )
                                                                {
                                                                        info.setState( UniqueMatchInfoBase::Gapped );
                                                                        
                                                                        info.setScore( completescore  );
                                                                        info.setPosition(rpos);
                                                                        
                                                                        GapInfo gapinfo;
                                                                        gapinfo.MINgap = MINgap;
                                                                        gapinfo.where = where;
                                                                        gapinfo.start = start;
                                                                        gapinfo.gap_pos = gap_pos;
                                                                        gapinfos[patid] = gapinfo;
                                                                }
                                                                else if ( info.getState() == UniqueMatchInfoBase::Gapped )
                                                                {
                                                                        // better match
                                                                        if ( completescore > info.getScore() + 1e-6 )
                                                                        {
                                                                                // larger score, use
                                                                                info.setScore( completescore  );                                                                
                                                                                info.setPosition(rpos);
                                                                                
                                                                                GapInfo gapinfo;
                                                                                gapinfo.MINgap = MINgap;
                                                                                gapinfo.where = where;
                                                                                gapinfo.start = start;
                                                                                gapinfo.gap_pos = gap_pos;
                                                                                gapinfos[patid] = gapinfo;
                                                                        }
                                                                        else if ( completescore < info.getScore() - 1e-6 )
                                                                        {
                                                                                // ignore, smaller score
                                                                        }
                                                                        else
                                                                        {
                                                                                if ( gapinfos.find(patid) != gapinfos.end() )
                                                                                        gapinfos.erase(gapinfos.find(patid));
                                                                                        
                                                                                #if 0
                                                                                // info.setState(UniqueMatchInfoBase::NonUnique);
                                                                                #endif
                                                                        }
                                                                }
                                                                
                                                                // std::cerr << "Score = " << MAXscore << " Gap = " << MINgap << " Where = " << where << " Gap_Pos = " << gap_pos <<std::endl;
                                                        }
                                                        else
                                                        {
                                                                // std::cerr << "No gap found." << std::endl;        
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
}
#endif
