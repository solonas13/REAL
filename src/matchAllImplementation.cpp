#if ! defined(MATCHALLIMPLEMENTATION_HPP)
#define MATCHALLIMPLEMENTATION_HPP

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

#include "real_config.hpp"

#include "PatternIdReader.hpp"
#include "MatcherBase.hpp"
#include "match.hpp"
#include "ListSetBlockReader.hpp"
#include "SignatureConstruction.hpp"
#include "RestWordBuffer.hpp"
#include "RangeVector.hpp"
#include "types.hpp"
#include "UniqueMatchInfo.hpp"
#include "matchAll.hpp"
#include "countReads.hpp"
#include "Mask.hpp"
#include "u_sort.hpp"
#include "getHistSize.hpp"
#include "getRadixSortTemp.hpp"
#include "AutoTextArray.hpp"
#include "RestMatch.hpp"
#include "ComputeScore.hpp"
#include "getFileList.hpp"
#include "getText.hpp"
#include "MapTextFile.hpp"
#include "getFileID.hpp"
#include "RealTimeClock.hpp"
#include "acgtnMap.hpp"
#include "getLookupTable.hpp"
#include "Lock.hpp"
#include "getNumLists.hpp"

#include "SynchronousReader.hpp"
#include "AsynchronousReader.hpp"
#include "FastAReader.hpp"
#include "FastQReader.hpp"
#include "FastFileDecoder.hpp"
#include "FastQualityFileDecoder.hpp"
#include "AsynchronousWriter.hpp"

#include <memory>
#include <fstream>

#include "real_config.hpp"

#include "RestWordBuffer.hpp"
#include "RangeVector.hpp"
#include "types.hpp"
#include "UniqueMatchInfo.hpp"
#include "matchUnique.hpp"
#include "countReads.hpp"
#include "Mask.hpp"
#include "u_sort.hpp"
#include "getHistSize.hpp"
#include "getRadixSortTemp.hpp"
#include "AutoTextArray.hpp"
#include "RestMatch.hpp"
#include "ComputeScore.hpp"
#include "getFileList.hpp"
#include "SignatureConstruction.hpp"
#include "getText.hpp"
#include "MapTextFile.hpp"
#include "getFileID.hpp"
#include "RealTimeClock.hpp"
#include "acgtnMap.hpp"
#include "getLookupTable.hpp"
#include "Lock.hpp"

#include "SynchronousReader.hpp"
#include "AsynchronousReader.hpp"
#include "FastAReader.hpp"
#include "FastQReader.hpp"
#include "FastFileDecoder.hpp"
#include "FastQualityFileDecoder.hpp"
#include "AsynchronousWriter.hpp"

#include <memory>
#include <fstream>

struct MatchPosAndError
{
        u_int64_t pos;
        u_int64_t file;
        u_int64_t frag;
        bool inverted;
        unsigned int k;
        double score;
        
        MatchPosAndError() {}
        MatchPosAndError(
                bool const rinverted,
                u_int64_t const rfileid,
                u_int64_t const rpos,
                unsigned int const rk,
                double const rscore,
                u_int64_t const rfrag)
        : pos(rpos), 
          file(rfileid),
          frag(rfrag),
          inverted(rinverted), k(rk), score(rscore) {}
};

inline bool operator< ( MatchPosAndError const & A, MatchPosAndError const & B)
{
        if ( A.k != B.k )
                return A.k < B.k;
        else if ( A.pos != B.pos )
                return A.pos < B.pos;
        else if ( A.file != B.file )
                return A.file < B.file;
        else if ( A.frag != B.frag )
                return A.frag < B.frag;
        else if ( A.score != B.score )
                return A.score < B.score;
        else
                return A.inverted < B.inverted;
}

inline bool operator!= ( MatchPosAndError const & A, MatchPosAndError const & B)
{
        return 
                (A<B) || (B<A);
}

inline bool operator== ( MatchPosAndError const & A, MatchPosAndError const & B)
{
        return !(A!=B);

}

inline void unifyMatches(std::vector<MatchPosAndError> & V)
{
	std::sort(V.begin(),V.end());
			
	std::vector < MatchPosAndError > Vu;
	Vu.push_back(V[0]);
			
	for ( unsigned int i = 1; i < V.size(); ++i )
	        if ( V[i] != V[i-1] )
		        Vu.push_back(V[i]);
        Vu.swap(V);
}

struct VectorUpdater
{
        typedef std::vector < MatchPosAndError > info_type;

        static bool matchGaps(info_type const & /* info */)
        {
                return false;
        }

        static inline void update(
                bool const inverted,
                size_t const fileid,
                unsigned int const pos,
                unsigned int const totalk,
                float const score,
                float const,
                unsigned int const fragid,
                info_type & V
        )
        {
                V.push_back(MatchPosAndError(inverted,fileid,pos,totalk,score,fragid));
        }                                                                                                                
};

template<bool scores>
std::pair < u_int64_t, u_int64_t > getReadMemory(u_int64_t const num_reads)
{
        // output number
        std::cerr << "number of reads " << num_reads << std::endl;
	// compute memory used for reads
	u_int64_t const read_memory = num_reads * sizeof(UniqueMatchInfo<scores>);
	std::cerr << "memory used for reads " << (read_memory/(1024*1024)) << " mb" << std::endl;
	//
	return std::pair<u_int64_t, u_int64_t>(num_reads,read_memory);
}

static u_int64_t getNonListMem(std::pair < u_int64_t, u_int64_t > const & read_info)
{
        // lookup table parameters 
	unsigned int const histsize = getHistSize();
        unsigned int const num_lists = getNumLists();
        // compute memory used for lookup tables
        size_t const lookupmem = num_lists * histsize * 2* sizeof(size_t);
        // output
        // std::cerr << "Using " << lookupmem << " bytes of memory for lookup tables." << std::endl;
        // memory we are using for structure which are not the lists
        size_t const non_list_mem = read_info.second + lookupmem;
        
        return non_list_mem;
}

template<typename signature_type, typename ptr_type>
static double getListElemSize()
{
        return
                8*
                (
                (getNumLists()>>1) * sizeof(BaseMask<signature_type,ptr_type>) + 
                ((getNumLists()>>1) + getRadixSortTemp()) * sizeof(Mask<signature_type,ptr_type>)
                )
                ;
}

template<typename signature_type, typename ptr_type>
static u_int64_t getNListMax(
        RealOptions const & opts,
        std::pair < u_int64_t, u_int64_t > const & read_info,
        u_int64_t const ATAsize
)
{
        if ( getNonListMem(read_info) + ATAsize > opts.usemem )
                throw std::bad_alloc();

        // number of elements per list
        return static_cast<u_int64_t>( (8.0*(opts.usemem - getNonListMem(read_info) - ATAsize) ) / getListElemSize<signature_type,ptr_type>());
}

template<
        typename signature_type, bool sse4,
        typename ptr_type, typename pattern_type,
        bool scores
>
struct AllMatcher : public MatcherBase<signature_type,sse4,ptr_type,pattern_type,scores>
{
        typedef MatcherBase<signature_type,sse4,ptr_type,pattern_type,scores> base;

        AllMatcher(
                u_int64_t const n_list, 
                RealOptions const & ropts, 
                SignatureConstruction<signature_type> const & rSC,
                MapTextFile<signature_type,sse4> & rMTF,
                Scoring const & rscoring,
                AutoTextArray<sse4> const & rATA,
                RangeVector<sse4> const & rRV
        ) : base(n_list,ropts,rSC,rMTF,rscoring,rATA,rRV)
        {
        }

        inline bool match(
                pattern_type const & pattern, 
                u_int64_t const fi,
                RestWordBuffer<sse4> & RWB,
                u_int64_t & handled,
                std::vector < MatchPosAndError > & localmatches
                ) const
        {
                std::vector < MatchPosAndError > V;
                
                unsigned int const patl = pattern.getPatternLength();

                if ( patl < static_cast<unsigned int>(base::opts.seedl) )
                {
                        std::cerr << "Skipping pattern " << toollib::remapString(std::string(pattern.mapped,pattern.mapped+patl)) << " shorter than seed length." << std::endl;
                        handled++;
                        return false;
                }
                
                bool pattern_dontcarefree = true;                                                
                for ( unsigned int i = 0; i < patl; ++i )
                        if ( pattern.mapped[i] > 3 )
                                pattern_dontcarefree = false;
                
                if ( ! pattern_dontcarefree )
                {
                        handled++;
                        return false;
                }
                
                RWB.setup(patl);
        
                u_int32_t m[4];

                if ( base::SC.signatureMapped(pattern.mapped,&m[0]) )
                {
                        // fill straight rest word array (pattern bits beyond seed)
                        RWB.setupStraight(pattern.mapped);

                        // straight pattern fragment sets
                        signature_type const straight[] = { base::SC.s0(m[0],m[1]), base::SC.s1(m[0],m[2]), base::SC.s2(m[0],m[3]), base::SC.s3(m[1],m[2]), base::SC.s4(m[1],m[3]), base::SC.s5(m[2],m[3]) };
                        
                        typedef Mask<signature_type,ptr_type> a_mask_type;
                        typedef BaseMask<signature_type,ptr_type> b_mask_type;

			float const epsilon = base::opts.getFilterValue(patl);

                        ::match< sse4, a_mask_type, b_mask_type, char const *, pattern_type, scores,VectorUpdater>(
                                base::fulllists[0],
                                base::baselists[2],
                                base::opts.seedkmax,
                                base::opts.totalkmax,
                                straight[0],
                                straight[5],
                                false,
                                fi,
                                base::SC.s0shift(getSampleBits()),
                                base::lookup[0],
                                base::RV,
                                base::ATA,
                                RWB,
                                pattern,base::scoring,
                                epsilon,
                                localmatches
                                );


                        ::match<sse4, a_mask_type, b_mask_type, char const *, pattern_type, scores,VectorUpdater>(base::fulllists[1],base::baselists[1],base::opts.seedkmax,base::opts.totalkmax,straight[1],straight[4],false,fi,base::SC.s1shift(getSampleBits()),base::lookup[1],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,localmatches);
                        ::match<sse4, a_mask_type, b_mask_type, char const *, pattern_type, scores,VectorUpdater>(base::fulllists[2],base::baselists[0],base::opts.seedkmax,base::opts.totalkmax,straight[2],straight[3],false,fi,base::SC.s2shift(getSampleBits()),base::lookup[2],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,localmatches);
                        ::match<sse4, b_mask_type, a_mask_type, char const *, pattern_type, scores,VectorUpdater>(base::baselists[0],base::fulllists[2],base::opts.seedkmax,base::opts.totalkmax,straight[3],straight[2],false,fi,base::SC.s3shift(getSampleBits()),base::lookup[3],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,localmatches);
                        ::match<sse4, b_mask_type, a_mask_type, char const *, pattern_type, scores,VectorUpdater>(base::baselists[1],base::fulllists[1],base::opts.seedkmax,base::opts.totalkmax,straight[4],straight[1],false,fi,base::SC.s4shift(getSampleBits()),base::lookup[4],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,localmatches);
                        ::match<sse4, b_mask_type, a_mask_type, char const *, pattern_type, scores,VectorUpdater>(base::baselists[2],base::fulllists[0],base::opts.seedkmax,base::opts.totalkmax,straight[5],straight[0],false,fi,base::SC.s5shift(getSampleBits()),base::lookup[5],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,localmatches);

                        // compute reverse transpose signature seed
                        u_int32_t im[4] = { 0,0,0,0 };
                        
                        base::SC.reverseMappedSignature ( pattern.mapped, &im[0] );
                        // reverse pattern fragment sets
                        signature_type const reverse[] = { base::SC.s0(im[0],im[1]), base::SC.s1(im[0],im[2]), base::SC.s2(im[0],im[3]), base::SC.s3(im[1],im[2]), base::SC.s4(im[1],im[3]), base::SC.s5(im[2],im[3]) };

                        // fill reverse rest word array (pattern bits beyond seed)
                        RWB.setupReverse(pattern.mapped);

                        ::match<sse4, a_mask_type, b_mask_type, char const *, pattern_type, scores,VectorUpdater>(base::fulllists[0],base::baselists[2],base::opts.seedkmax,base::opts.totalkmax,reverse[0],reverse[5],true,fi,base::SC.s0shift(getSampleBits()),base::lookup[0],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,localmatches);
                        ::match<sse4, a_mask_type, b_mask_type, char const *, pattern_type, scores,VectorUpdater>(base::fulllists[1],base::baselists[1],base::opts.seedkmax,base::opts.totalkmax,reverse[1],reverse[4],true,fi,base::SC.s1shift(getSampleBits()),base::lookup[1],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,localmatches);
                        ::match<sse4, a_mask_type, b_mask_type, char const *, pattern_type, scores,VectorUpdater>(base::fulllists[2],base::baselists[0],base::opts.seedkmax,base::opts.totalkmax,reverse[2],reverse[3],true,fi,base::SC.s2shift(getSampleBits()),base::lookup[2],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,localmatches);
                        ::match<sse4, b_mask_type, a_mask_type, char const *, pattern_type, scores,VectorUpdater>(base::baselists[0],base::fulllists[2],base::opts.seedkmax,base::opts.totalkmax,reverse[3],reverse[2],true,fi,base::SC.s3shift(getSampleBits()),base::lookup[3],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,localmatches);
                        ::match<sse4, b_mask_type, a_mask_type, char const *, pattern_type, scores,VectorUpdater>(base::baselists[1],base::fulllists[1],base::opts.seedkmax,base::opts.totalkmax,reverse[4],reverse[1],true,fi,base::SC.s4shift(getSampleBits()),base::lookup[4],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,localmatches);
                        ::match<sse4, b_mask_type, a_mask_type, char const *, pattern_type, scores,VectorUpdater>(base::baselists[2],base::fulllists[0],base::opts.seedkmax,base::opts.totalkmax,reverse[5],reverse[0],true,fi,base::SC.s5shift(getSampleBits()),base::lookup[5],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,localmatches);
                }	

                handled++;

                return true;
        }
};


template<bool sse4, typename signature_type, typename reader_type, bool scores>
void EnumerateAllMatches<sse4,signature_type,reader_type,scores>::doMatching(RealOptions const & opts)
{
        typedef typename reader_type::pattern_type pattern_type;

#if defined(_OPENMP)
        unsigned int const numthreads = omp_get_max_threads();
        std::cerr << "Starting matching for all hits using " << numthreads << " matching threads." << std::endl;
#else
        unsigned int const numthreads = 1;
#endif
        
	u_int64_t const numpat = reader_type::countPatterns(opts.patternfilename);
	
        // get read information
        std::pair < u_int64_t, u_int64_t > const read_info = getReadMemory<scores>(numpat);
        
        // get filenames of input files
        std::vector< std::string > filenames;
        std::string const & textfilesuffix = ".fa";
        getFileList(opts.textfilename,filenames,textfilesuffix);
        
        SignatureConstruction<signature_type> const SC(opts.seedl,opts.nu);
        
        // double const cps = CLOCKS_PER_SEC;

        int const qualityOffset = opts.qualityOffset ? opts.qualityOffset : reader_type::getOffset(opts.patternfilename);
        
        if ( ! qualityOffset )
                throw std::runtime_error("Unable to automatically detect FastQ quality format.");
        
        Scoring const scoring(opts.similarity, opts.gc, opts.trans, opts.err, opts.gcmut_bias);

        std::auto_ptr < AsynchronousWriter  > output;
        
        if ( opts.outputfilename == "-" )
        {
		output = std::auto_ptr < AsynchronousWriter  >( new AsynchronousWriter(STDOUT_FILENO,16) );
	}
	else
	{
		output = std::auto_ptr < AsynchronousWriter  >( new AsynchronousWriter(opts.outputfilename,16) );
	}

        for ( unsigned int fi = 0; fi < filenames.size(); ++fi )
        {
                bool const lastfile = (fi+1)==filenames.size();

                std::cerr << "Processing file " << filenames[fi] << (lastfile?" (last processed file)":"")<< std::endl;
                
                std::vector < std::pair < std::string, u_int64_t > > ranges;
                std::auto_ptr < AutoTextArray<sse4> > AATA = getText<sse4>(filenames[fi],ranges);
                
                AutoTextArray<sse4> const & ATA = *(AATA.get());
                RangeVector<sse4> RV(ranges);
                
                if ( ATA.getN() < static_cast<unsigned int>(opts.seedl) )
                {
                        std::cerr << "File " << filenames[fi] << " is too small for seed length, skipping it." << std::endl;
                        continue;
                }
                
                std::cerr << ATA.size() << " bytes allocated for text." << std::endl;

                // pointer type in masks
                typedef u_int32_t ptr_type;

                // number of element per list
                size_t const n_list_max = getNListMax<signature_type,ptr_type>(opts,read_info,ATA.size() + (ATA.getN()/8) );
        
                std::cerr << "We have space for " << n_list_max << " list elements." << std::endl;
        
                if ( ! n_list_max )
                        throw std::bad_alloc();

                // open text file
                MapTextFile<signature_type,sse4> MTF(ATA, opts.seedl, opts.nu);
                
                u_int64_t const filesize = ATA.getN();
                std::cerr << "File size = " << filesize << std::endl;
                u_int64_t const expblocks = (filesize - opts.seedl + 1 + (n_list_max-1)) / n_list_max;
                std::cerr << "Expected number of blocks = " << expblocks << std::endl;
                u_int64_t const n_list = (filesize + (expblocks-1)) / expblocks;
                std::cerr << "Using n_list = " << n_list << std::endl;

                unsigned int masks = 0;

                AllMatcher<signature_type,sse4,ptr_type,pattern_type,scores> AM(n_list,opts,SC,MTF,scoring,ATA,RV);
                
                /**
                 * read blocks of input file
                 **/
                while ( (masks = AM.readNextBlock()) > 0 )
                {        
                        PatternIdReader<reader_type> patfile(opts.patternfilename, qualityOffset, numthreads);
                
                        u_int64_t totalhandled = 0;

                        toollib::Lock stderrlock;

#if defined(_OPENMP)                        
                        #pragma omp parallel
#endif
                        {
                                RestWordBuffer<sse4> RWB(opts.seedl);
                                std::pair < typename reader_type::block_type *, typename reader_type::idblock_type * > block;
                                std::auto_ptr < std::ostringstream > tempostr(new std::ostringstream);

                                while ( patfile.getBlock(block) )
                                {
                                        u_int64_t handled = 0;
                                
                                        for ( u_int64_t z = 0; z < block.first->blocksize; ++z )
                                        {
                                                pattern_type const & pattern = block.first->getPattern(z);
                                                std::vector < MatchPosAndError > localmatches;
                                                AM.match(pattern, fi, RWB, handled,localmatches);

                                                if ( localmatches.size() )
                                                {
                                                        unifyMatches(localmatches);
                                                        
                                                        for ( unsigned int ll = 0; ll < localmatches.size(); ++ll )
                                                        {
                                                                MatchPosAndError const & M = localmatches[ll];

                                                                (*tempostr)
                                                                        << block.second->ids[z] << "\t" 
                                                                        << (M.inverted ? 
                                                                                toollib::remapString(
                                                                                        std::string(
                                                                                                pattern.transposed,
                                                                                                pattern.transposed + pattern.getPatternLength()) ) : 
                                                                                toollib::remapString(
                                                                                        std::string(pattern.mapped, pattern.mapped+pattern.getPatternLength() )
                                                                                        )
                                                                                ) 
                                                                                << "\t";

                                                                if ( scores )
                                                                        (*tempostr) << ComputeScore<sse4,pattern_type,scores>::computeScore(M.inverted,ATA,pattern,scoring,M.pos,pattern.getPatternLength());
                                                                
                                                                (*tempostr) << "\t"
                                                                        << 1 << "\t"
                                                                        << "a" << "\t"
                                                                        << pattern.getPatternLength() << "\t"
                                                                        << (M.inverted ? "-" : "+") << "\t"
                                                                        << RV.positionToId(M.pos) << "\t"
                                                                        << M.pos - ranges[RV.positionToRange(M.pos)].second +1 << "\t"
                                                                        << /* type of hit << */ "\t"
                                                                        << M.k
                                                                        << std::endl;                                                                        
                                                                        
                                                                if ( tempostr->str().size() > 16384 )
                                                                {
                                                                        std::string tempstring = tempostr->str();
                                                                        output->write ( tempstring.begin(), tempstring.end() );
                                                                        tempostr = std::auto_ptr < std::ostringstream > ( new std::ostringstream );
                                                                }
                                                        }
                                                }

                                        }
                                        
                                        patfile.returnBlock(block);

#if defined(_OPENMP)
                                        #pragma omp critical
#endif
                                        {
                                                totalhandled += handled;
                                                std::cerr << "\r                                                              \r" << (static_cast<double>(totalhandled)/numpat) << std::flush;
                                        }
                                }
                        }
                        std::cerr << std::endl;
                }
                std::cerr << "All done." << std::endl;
        }
}
#endif
