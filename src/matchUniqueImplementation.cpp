#if ! defined(MATCHUNIQUEIMPLEMENTATION_HPP)
#define MATCHUNIQUEIMPLEMENTATION_HPP

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

#include "BinRadixSort.hpp"
#include "PatternIdReader.hpp"
#include "MatcherBase.hpp"
#include "match.hpp"
#include "ListSetBlockReader.hpp"
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
#include "RangeSet.hpp"
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

/**
 * memory:
 *
 * text
 * range vector + rank dictionary
 * lookup table
 * Unique match info
 * temp mem for radix sort
 * lists
 **/
template<bool scored_matching>
struct UpdateUniqueInfo
{

};

template<>
struct UpdateUniqueInfo<false>
{
        static bool const scores = false;

        typedef UniqueMatchInfo<scores> info_type;
        
        static bool matchGaps(info_type const & info)
        {
                return 
                	(info.getState() == UniqueMatchInfoBase::NoMatch)
                	||
                	(info.getState() == UniqueMatchInfoBase::Gapped)
                	;
        }

        static void update(
                bool const inverted,
                size_t const fileid,
                unsigned int const pos,
                unsigned int const totalk,
                float const,
                float const,
                unsigned int const fragid,
                info_type & info
                )
        {
                switch ( info.getState() )
                {
                        case UniqueMatchInfoBase::NoMatch:
                        case UniqueMatchInfoBase::Gapped:
                        {
                                info.setState(inverted ? UniqueMatchInfoBase::Reverse : UniqueMatchInfoBase::Straight);
                                info.setPosition(pos);
                                info.setFileid(fileid);
                                info.setErrors(totalk);
                                info.setFragment(fragid);
                                break;
                        }
                        case UniqueMatchInfoBase::Straight:
                        case UniqueMatchInfoBase::Reverse:
                        {
                                // found a match with less errors
                                if ( totalk < info.getErrors() )
                                {
                                        info.setState(inverted ? UniqueMatchInfoBase::Reverse : UniqueMatchInfoBase::Straight);
                                        info.setPosition(pos);
                                        info.setFileid(fileid);
                                        info.setErrors(totalk);                
                                        info.setFragment(fragid);
                                }
                                else if ( 
                                        totalk == info.getErrors() 
                                        && 
                                        ( (pos != info.getPosition()) || (fileid != info.getFileid()) || (fragid != info.getFragment()) )
                                )
                                {
                                        info.setState(UniqueMatchInfoBase::NonUnique);
                                }
                                // more errors than already known, ignore
                                else
                                {
                                }
                                break;
                        }
                        case UniqueMatchInfoBase::NonUnique:
                        {
                                // less errors than seen before, take
                                if ( totalk < info.getErrors() )
                                {
                                        info.setState(inverted ? UniqueMatchInfoBase::Reverse : UniqueMatchInfoBase::Straight);
                                        info.setPosition(pos);
                                        info.setFileid(fileid);
                                        info.setErrors(totalk); 
                                        info.setFragment(fragid);               
                                }
                                break;
                        }
                }
        }
};

template<>
struct UpdateUniqueInfo<true>
{
        static bool const scores = true;

        typedef UniqueMatchInfo<scores> info_type;

        static bool matchGaps(info_type const & info)
        {
        	return
                	(info.getState() == UniqueMatchInfoBase::NoMatch)
                	||
                	(info.getState() == UniqueMatchInfoBase::Gapped)
                	;
        }

        static void update(
                bool const inverted,
                size_t const fileid,
                unsigned int const pos,
                unsigned int const totalk,
                float const score,
                float const epsilon,
                unsigned int const fragid,
                info_type & info
                )
        {
                switch ( info.getState() )
                {
                        case UniqueMatchInfoBase::Gapped:
                        case UniqueMatchInfoBase::NoMatch:
                        {
                                info.setState(inverted ? UniqueMatchInfoBase::Reverse : UniqueMatchInfoBase::Straight);
                                info.setPosition(pos);
                                info.setFileid(fileid);
                                info.setErrors(totalk);
                                info.setScore(score);
                                info.setFragment(fragid);
                                break;
                        }
                        case UniqueMatchInfoBase::Straight:
                        case UniqueMatchInfoBase::Reverse:
                        {
                                // greater score?
                                if ( (score > info.getScore()+epsilon) )
                                {
                                        info.setState(inverted ? UniqueMatchInfoBase::Reverse : UniqueMatchInfoBase::Straight);
                                        info.setPosition(pos);
                                        info.setFileid(fileid);
                                        info.setErrors(totalk);                
                                        info.setScore(score);                                                                
                                        info.setFragment(fragid);
                                }
                                // same score, different position
                                else if ( 
                                        (score > info.getScore()-epsilon) 
                                        && 
                                        ( (pos != info.getPosition()) || (fileid != info.getFileid()) || (fragid != info.getFragment()) )
                                )
                                {
                                        info.setState(UniqueMatchInfoBase::NonUnique);
                                }
                                // smaller score, ignore
                                else
                                {
                                
                                }
                                break;
                        }
                        case UniqueMatchInfoBase::NonUnique:
                        {
                                // greater score?
                                if ( score > info.getScore()+epsilon )
                                {
                                        info.setState(inverted ? UniqueMatchInfoBase::Reverse : UniqueMatchInfoBase::Straight);
                                        info.setPosition(pos);
                                        info.setFileid(fileid);
                                        info.setErrors(totalk);                
                                        info.setScore(score);   
                                        info.setFragment(fragid);                                                          
                                }
                                break;
                        }
                }

        }

};

template<bool scores, typename output_type, typename pattern_type>
void printMatchUnlocked(
        UniqueMatchInfo<scores> const & info,
        output_type & output,
        pattern_type const & pattern,
        std::string const & id,
        unsigned int const patl,
        u_int64_t & unique,
        RangeSet const & RS
        )
{
        switch ( info.getState() )
        {
                case UniqueMatchInfoBase::Straight:
                        {
                                output
                                        << id << "\t" 
                                        << toollib::remapString(
                                                std::string(
                                                        pattern.mapped,
                                                        pattern.mapped + patl)) << "\t";
                                
                                if ( scores )
                                        output << (info.getScore());
                                
                                output
                                        <<  "\t"
                                        << 1 << "\t"
                                        << "a" << "\t"
                                        << patl << "\t"
                                        << "+" << "\t"
                                        << RS.ranges[info.getFileid()][info.getFragment()].first << "\t"
                                        // << fileids[info.getFileid()] << "\t"
                                        << info.getPosition() - RS.ranges[info.getFileid()][info.getFragment()].second  + 1 << "\t"
                                        << /* type of hit << */ "\t"
                                        << info.getErrors();
                                output << "\n";

                                unique++;                                
                        }
                        break;
                case UniqueMatchInfoBase::Reverse:
                        {
                                output
                                        << id << "\t" 
                                        << toollib::remapString(
                                                std::string(
                                                        pattern.transposed,pattern.transposed+patl)) << "\t";
                                if ( scores )
                                        output << (info.getScore());
                                        
                                output << "\t"
                                        << 1 << "\t"
                                        << "a" << "\t"
                                        << patl << "\t"
                                        << "-" << "\t"
                                        << RS.ranges[info.getFileid()][info.getFragment()].first << "\t"
                                        //<< fileids[info.getFileid()] << "\t"
                                        << info.getPosition() - RS.ranges[info.getFileid()][info.getFragment()].second + 1 << "\t"
                                        << /* type of hit << */ "\t"
                                        << info.getErrors();
                                output << "\n";

                                unique++;
                        }
                        break;
                default:
                        break;        	        
        }                                                        
}

template<bool scores, typename output_type, typename pattern_type>
void printMatch(
        UniqueMatchInfo<scores> const & info,
        output_type & output,
        toollib::Lock & stderrlock,
        pattern_type const & pattern,
        std::string const & id,
        unsigned int const patl,
        u_int64_t & unique,
        RangeSet const & RS
        )
{
        switch ( info.getState() )
        {
                case UniqueMatchInfoBase::Straight:
                case UniqueMatchInfoBase::Reverse:
                        stderrlock.lock();
                        printMatchUnlocked(info,output,pattern,id,patl,unique,RS);
                        stderrlock.unlock();
                        break;
                default:
                        break;
        }                                                        
}

template<
        typename signature_type, bool sse4,
        typename ptr_type, typename pattern_type,
        bool scores
>
struct UniqueMatcher : public MatcherBase<signature_type,sse4,ptr_type,pattern_type,scores>
{
        typedef MatcherBase<signature_type,sse4,ptr_type,pattern_type,scores> base;

        UniqueMatcher(
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
                UniqueMatchInfo<scores> & info, 
                u_int64_t const fi,
                RestWordBuffer<sse4> & RWB,
                u_int64_t & handled) const
        {
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
                        
                        typedef Mask<signature_type, ptr_type> mask_type;
                        typedef BaseMask<signature_type, ptr_type> base_mask_type;

			float const epsilon = base::opts.getFilterValue(patl);
			
			// std::cerr << "patl=" << patl << " epsilon=" << epsilon << std::endl;
                        
                        // try to uniquely match straight pattern
                        ::match<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(
                                base::fulllists[0],base::baselists[2],
                                base::opts.seedkmax,
                                base::opts.totalkmax,
                                straight[0],straight[5],
                                false /* straight */,
                                fi,
                                base::SC.s0shift(getSampleBits()),
                                base::lookup[0],
                                base::RV,
                                base::ATA,
                                RWB,
                                pattern, 
                                base::scoring,
                                epsilon,
                                info
                        );

                        bool const uni0s = (info.getState() == UniqueMatchInfoBase::Straight) && (info.getErrors() == 0);
                        
                        if ( (! uni0s) || scores )
                        {
                                ::match<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::fulllists[1],base::baselists[1],base::opts.seedkmax,base::opts.totalkmax,straight[1],straight[4],false,fi,base::SC.s1shift(getSampleBits()),base::lookup[1],base::RV,base::ATA,RWB
                                        ,pattern,base::scoring,epsilon,info
                                );
                                ::match<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::fulllists[2],base::baselists[0],base::opts.seedkmax,base::opts.totalkmax,straight[2],straight[3],false,fi,base::SC.s2shift(getSampleBits()),base::lookup[2],base::RV,base::ATA,RWB
                                        ,pattern,base::scoring,epsilon,info
                                );
                                ::match<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[0],base::fulllists[2],base::opts.seedkmax,base::opts.totalkmax,straight[3],straight[2],false,fi,base::SC.s3shift(getSampleBits()),base::lookup[3],base::RV,base::ATA,RWB
                                        ,pattern,base::scoring,epsilon,info
                                );
                                ::match<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[1],base::fulllists[1],base::opts.seedkmax,base::opts.totalkmax,straight[4],straight[1],false,fi,base::SC.s4shift(getSampleBits()),base::lookup[4],base::RV,base::ATA,RWB
                                        ,pattern,base::scoring,epsilon,info
                                );
                                ::match<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[2],base::fulllists[0],base::opts.seedkmax,base::opts.totalkmax,straight[5],straight[0],false,fi,base::SC.s5shift(getSampleBits()),base::lookup[5],base::RV,base::ATA,RWB
                                        ,pattern,base::scoring,epsilon,info
                                );
                        }

                        // compute reverse transpose signature seed
                        u_int32_t im[4] = { 0,0,0,0 };
                        
                        base::SC.reverseMappedSignature ( pattern.mapped, &im[0] );
                        // reverse pattern fragment sets
                        signature_type const reverse[] = { base::SC.s0(im[0],im[1]), base::SC.s1(im[0],im[2]), base::SC.s2(im[0],im[3]), base::SC.s3(im[1],im[2]), base::SC.s4(im[1],im[3]), base::SC.s5(im[2],im[3]) };

                        // fill reverse rest word array (pattern bits beyond seed)
                        RWB.setupReverse(pattern.mapped);

                        // try to uniquely match reverse pattern
                        ::match<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::fulllists[0],base::baselists[2],base::opts.seedkmax,base::opts.totalkmax,reverse[0],reverse[5],true,fi,base::SC.s0shift(getSampleBits()),base::lookup[0],base::RV,base::ATA,RWB
                                ,pattern,base::scoring,epsilon,info
                        );

                        bool const uni0r = (info.getState() == UniqueMatchInfoBase::Reverse) && (info.getErrors() == 0);

                        if ( (! uni0r) || scores )
                        {
                                ::match<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::fulllists[1],base::baselists[1],base::opts.seedkmax,base::opts.totalkmax,reverse[1],reverse[4],true,fi,base::SC.s1shift(getSampleBits()),base::lookup[1],base::RV,base::ATA,RWB
                                        ,pattern,base::scoring,epsilon,info
                                );
                                ::match<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::fulllists[2],base::baselists[0],base::opts.seedkmax,base::opts.totalkmax,reverse[2],reverse[3],true,fi,base::SC.s2shift(getSampleBits()),base::lookup[2],base::RV,base::ATA,RWB
                                        ,pattern,base::scoring,epsilon,info
                                );
                                ::match<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[0],base::fulllists[2],base::opts.seedkmax,base::opts.totalkmax,reverse[3],reverse[2],true,fi,base::SC.s3shift(getSampleBits()),base::lookup[3],base::RV,base::ATA,RWB
                                        ,pattern,base::scoring,epsilon,info
                                );
                                ::match<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[1],base::fulllists[1],base::opts.seedkmax,base::opts.totalkmax,reverse[4],reverse[1],true,fi,base::SC.s4shift(getSampleBits()),base::lookup[4],base::RV,base::ATA,RWB
                                        ,pattern,base::scoring,epsilon,info
                                );
                                ::match<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[2],base::fulllists[0],base::opts.seedkmax,base::opts.totalkmax,reverse[5],reverse[0],true,fi,base::SC.s5shift(getSampleBits()),base::lookup[5],base::RV,base::ATA,RWB
                                        ,pattern,base::scoring,epsilon,info
                                );
                        }

                        #if 0
                        if ( lastfilelastblock )
                                printMatch<scores>(info,output,stderrlock,pattern,patl,unique,RS);
                        #endif
                }	

                handled++;

                return true;
        }
        inline bool matchGaps(
                pattern_type const & pattern, 
                UniqueMatchInfo<scores> & info, 
                RestWordBuffer<sse4> & RWB,
                std::map<unsigned int, GapInfo> & gapinfos,
                u_int64_t & handled) const
        {
                if ( info.getState() == UniqueMatchInfoBase::NoMatch || info.getState() == UniqueMatchInfoBase::Gapped )
                {
                        unsigned int const patl = pattern.getPatternLength();
                        unsigned int const patid = pattern.getPatID();

                        if ( patl < static_cast<unsigned int>(base::opts.seedl) )
                        {
                                std::cerr << "Skipping pattern " << toollib::remapString(std::string(pattern.mapped,pattern.mapped+patl)) << "shorter than seed length." << std::endl;
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
                                // straight pattern fragment sets
                                signature_type const straight[] = { base::SC.s0(m[0],m[1]), base::SC.s1(m[0],m[2]), base::SC.s2(m[0],m[3]), base::SC.s3(m[1],m[2]), base::SC.s4(m[1],m[3]), base::SC.s5(m[2],m[3]) };
                                
                                typedef Mask<signature_type, ptr_type> mask_type;
                                typedef BaseMask<signature_type, ptr_type> base_mask_type;
                                
                                // compute reverse transpose signature seed
                                u_int32_t im[4] = { 0,0,0,0 };
                                
                                base::SC.reverseMappedSignature ( pattern.mapped, &im[0] );
                                // reverse pattern fragment sets
                                signature_type const reverse[] = { base::SC.s0(im[0],im[1]), base::SC.s1(im[0],im[2]), base::SC.s2(im[0],im[3]), base::SC.s3(im[1],im[2]), base::SC.s4(im[1],im[3]), base::SC.s5(im[2],im[3]) };

                                float const epsilon = base::opts.getFilterValue(patl);

                                ::matchGaps<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::fulllists[0],base::baselists[2],base::opts.seedkmax,straight[0],straight[5],false /* straight */,base::SC.s0shift(getSampleBits()),base::lookup[0],base::RV,base::ATA,RWB,pattern, base::scoring,epsilon,info,patid,gapinfos);
                                ::matchGaps<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::fulllists[1],base::baselists[1],base::opts.seedkmax,straight[1],straight[4],false,base::SC.s1shift(getSampleBits()),base::lookup[1],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,info,patid,gapinfos);
                                ::matchGaps<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::fulllists[2],base::baselists[0],base::opts.seedkmax,straight[2],straight[3],false,base::SC.s2shift(getSampleBits()),base::lookup[2],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,info,patid,gapinfos);
                                ::matchGaps<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[0],base::fulllists[2],base::opts.seedkmax,straight[3],straight[2],false,base::SC.s3shift(getSampleBits()),base::lookup[3],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,info,patid,gapinfos);
                                ::matchGaps<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[1],base::fulllists[1],base::opts.seedkmax,straight[4],straight[1],false,base::SC.s4shift(getSampleBits()),base::lookup[4],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,info,patid,gapinfos);
                                ::matchGaps<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[2],base::fulllists[0],base::opts.seedkmax,straight[5],straight[0],false,base::SC.s5shift(getSampleBits()),base::lookup[5],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,info,patid,gapinfos);

                                ::matchGaps<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::fulllists[0],base::baselists[2],base::opts.seedkmax,reverse[0],reverse[5],true,base::SC.s0shift(getSampleBits()),base::lookup[0],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,info,patid,gapinfos);
                                ::matchGaps<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::fulllists[1],base::baselists[1],base::opts.seedkmax,reverse[1],reverse[4],true,base::SC.s1shift(getSampleBits()),base::lookup[1],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,info,patid,gapinfos);
                                ::matchGaps<sse4, mask_type, base_mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::fulllists[2],base::baselists[0],base::opts.seedkmax,reverse[2],reverse[3],true,base::SC.s2shift(getSampleBits()),base::lookup[2],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,info,patid,gapinfos);
                                ::matchGaps<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[0],base::fulllists[2],base::opts.seedkmax,reverse[3],reverse[2],true,base::SC.s3shift(getSampleBits()),base::lookup[3],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,info,patid,gapinfos);
                                ::matchGaps<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[1],base::fulllists[1],base::opts.seedkmax,reverse[4],reverse[1],true,base::SC.s4shift(getSampleBits()),base::lookup[4],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,info,patid,gapinfos);
                                ::matchGaps<sse4, base_mask_type, mask_type, char const *,pattern_type,scores,UpdateUniqueInfo<scores> >(base::baselists[2],base::fulllists[0],base::opts.seedkmax,reverse[5],reverse[0],true,base::SC.s5shift(getSampleBits()),base::lookup[5],base::RV,base::ATA,RWB,pattern,base::scoring,epsilon,info,patid,gapinfos);
                                
                        }
                }	

                handled++;

                return true;
        }
};

template<bool sse4, typename signature_type>
struct GappedFactorArrayConstIterator
{
	typedef signature_type value_type;
	typedef signature_type const * pointer;
	typedef signature_type const & reference;
	typedef int64_t difference_type;
	typedef std::random_access_iterator_tag iterator_category;


	AutoTextArrayGappedComparator<sse4> const * data;
	u_int32_t const * A;
	u_int64_t i;
		
	GappedFactorArrayConstIterator() : data(0), A(0), i(0) {}
	GappedFactorArrayConstIterator(
		AutoTextArrayGappedComparator<sse4> const * rdata, 
		u_int32_t const * rA,
		u_int64_t const ri
		) 
	: data(rdata), A(rA), i(ri) {}
	
	u_int64_t operator*() const
	{
		return (*data)(A[i]);
	}
	
	u_int64_t operator[](u_int64_t j) const
	{
		return (*data)(A[i+j]);
	}

	GappedFactorArrayConstIterator<sse4,signature_type> & operator++()
	{
		i++;
		return *this;
	}
	GappedFactorArrayConstIterator<sse4,signature_type> operator++(int)
	{
		GappedFactorArrayConstIterator<sse4,signature_type> temp = *this;
		i++;
		return temp;
	}
	GappedFactorArrayConstIterator<sse4,signature_type> & operator--()
	{
		i--;
		return *this;
	}
	GappedFactorArrayConstIterator<sse4,signature_type> operator--(int)
	{
		GappedFactorArrayConstIterator<sse4,signature_type> temp = *this;
		i--;
		return temp;
	}

	GappedFactorArrayConstIterator<sse4,signature_type> & operator+=(u_int64_t j)
	{
		i += j;
		return *this;
	}
	GappedFactorArrayConstIterator<sse4,signature_type> & operator-=(u_int64_t j)
	{
		i -= j;
		return *this;
	}
	
	bool operator<(GappedFactorArrayConstIterator<sse4,signature_type> const & I) const
	{
		return i < I.i;
	}
	bool operator==(GappedFactorArrayConstIterator<sse4,signature_type> const & I) const
	{
		return 
			(data == I.data) 
			&&
			(A == I.A)
			&&
			(i == I.i);
	}
	bool operator!=(GappedFactorArrayConstIterator<sse4,signature_type> const & I) const
	{
		return ! ( (*this) == I );
	}
};

template<bool sse4, typename signature_type>
inline GappedFactorArrayConstIterator<sse4,signature_type> operator+ ( GappedFactorArrayConstIterator<sse4,signature_type> const & I, u_int64_t j )
{
	GappedFactorArrayConstIterator<sse4,signature_type> J = I;
	J += j;
	return J;
}

template<bool sse4, typename signature_type>
inline GappedFactorArrayConstIterator<sse4,signature_type> operator- ( GappedFactorArrayConstIterator<sse4,signature_type> const & I, u_int64_t j )
{
	GappedFactorArrayConstIterator<sse4,signature_type> J = I;
	J -= j;
	return J;
}

template<bool sse4, typename signature_type>
inline typename GappedFactorArrayConstIterator<sse4,signature_type>::difference_type operator- ( GappedFactorArrayConstIterator<sse4,signature_type> const & I, GappedFactorArrayConstIterator<sse4,signature_type> const & J )
{
        return
                static_cast < typename GappedFactorArrayConstIterator<sse4,signature_type>::difference_type > (I.i) -
                static_cast < typename GappedFactorArrayConstIterator<sse4,signature_type>::difference_type > (J.i);
}

template<bool sse4, typename signature_type>
struct GappedFactorArray
{
        AutoTextArray<sse4> const & ATA;
        RangeVector<sse4> const & RV;
        RealOptions const & opts; 
        AutoTextArrayGappedComparator<sse4> comp;
        u_int64_t const n_list_max;
        AutoArray<u_int32_t> A;
        AutoArray<u_int32_t> lookup;
        u_int64_t ii;
        SignatureConstruction<signature_type> const & SC;
        Scoring const & scoring;
        
        GappedFactorArrayConstIterator<sse4,signature_type> begin;
        GappedFactorArrayConstIterator<sse4,signature_type> end;
        
        u_int64_t getMaxListSize(u_int64_t const uniqueinfomemory) const
        {
                u_int64_t const numgfa = 3;
        
                u_int64_t const textmemory = ATA.size();
                u_int64_t const rangevectormemory = RV.size();
                u_int64_t const lookupmemory = numgfa * 2 * getHistSize() * sizeof(u_int32_t);
                u_int64_t const trlmemory = textmemory + rangevectormemory + lookupmemory;
                u_int64_t const nonlistmemory = trlmemory + uniqueinfomemory;
                
                #if 0
                std::cerr << "text memory: " << textmemory << std::endl;
                std::cerr << "range vector memory: " << rangevectormemory << std::endl;
                std::cerr << "lookup memory: " << lookupmemory << std::endl;
                std::cerr << "trlmemory: " << trlmemory << std::endl;
                std::cerr << "uniqueinfomemory: " << uniqueinfomemory << std::endl;
                std::cerr << "nonlistmemory: " << nonlistmemory << std::endl;
                std::cerr << "usemem: " << opts.usemem << std::endl;
                #endif

                if ( nonlistmemory > opts.usemem )
                {
                        std::cerr << "Insufficient memory." << std::endl;
                        throw std::bad_alloc();
                }

                u_int64_t const restmemory = opts.usemem - nonlistmemory;
                u_int64_t const elementmemory = numgfa * sizeof(u_int32_t);
                u_int64_t const rn_list_max = (restmemory / elementmemory);

                #if 0
                std::cerr << "restmemory: " << restmemory << std::endl;
                std::cerr << "elementmemory: " << elementmemory << std::endl;
                #endif
                std::cerr << "n_list_max: " << rn_list_max << std::endl;

                return rn_list_max;        
        }
        
        GappedFactorArray(
                AutoTextArray<sse4> const & rATA, 
                RangeVector<sse4> const & rRV,
                RealOptions const & ropts,
                u_int64_t const gapsize,
                u_int64_t const uniqueinfomemory,
                SignatureConstruction<signature_type> const & rSC,
                Scoring const & rscoring
        )
        : ATA(rATA), RV(rRV), opts(ropts), 
          comp(ATA,opts.seedl/4,gapsize*(opts.seedl/4),opts.seedl/4),
          n_list_max(getMaxListSize(uniqueinfomemory)),
          A(n_list_max,false),
          ii(0),
          SC(rSC),
          scoring(rscoring)
        {
        }
        
        u_int64_t getNextBlock()
        {
                u_int64_t const c = ATA.fillDontCareFreeGappedTextWordPositions(A.get(),n_list_max,ii,comp.p,comp.g,comp.q);
                
		std::cerr << "[B" << c;
                BinRadixSort<u_int32_t *, AutoTextArrayGappedComparator<sse4> >::binradixsort(
                        A.get(),A.get()+c,static_cast<signature_type>(1) << (2*2*(opts.seedl/4)-1),comp);
		std::cerr << "]";

                lookup = getLookupTableProjected(
                        A.get(), 
                        c, 
                        SC.s0shift(getSampleBits()), 
                        getHistSize(), 
                        comp);

                #if 0
                for ( u_int64_t i = 0; i < c; ++i )
                {
                        u_int64_t k = A[i];
                
                        for ( u_int64_t j = 0; j < comp.p; ++j ) std::cerr << ATA[ k++ ];
                        std::cerr << "*";
                        for ( u_int64_t j = 0; j < comp.g; ++j ) std::cerr << ATA[ k++ ];
                        std::cerr << "*";
                        for ( u_int64_t j = 0; j < comp.q; ++j ) std::cerr << ATA[ k++ ];
                                
                        std::cerr << "\t" << ATA.getGappedTextWord(A[i],comp.p,comp.g,comp.q);
                                
                        std::cerr << std::endl;
                }
                
                std::cerr << "------------------------------------------" << std::endl;
                #endif

                #if 0
                for ( unsigned int i = 0; i < getHistSize(); ++i )
                {
                        u_int32_t const low = lookup[2*i];
                        u_int32_t const high = lookup[2*i+1];
                        
                        std::cerr << "+++++++++++++" << std::endl;
                        for ( unsigned int j = low; j != high; ++j )
                        {
                                u_int64_t k = A[j];
                        
                                for ( u_int64_t j = 0; j < comp.p; ++j ) std::cerr << ATA[ k++ ];
                                std::cerr << "*";
                                for ( u_int64_t j = 0; j < comp.g; ++j ) std::cerr << ATA[ k++ ];
                                std::cerr << "*";
                                for ( u_int64_t j = 0; j < comp.q; ++j ) std::cerr << ATA[ k++ ];                		        
                                
                                std::cerr << std::endl;
                        }
                        std::cerr << "-------------" << std::endl;
                }
                #endif

                begin = GappedFactorArrayConstIterator<sse4,signature_type>(&comp,A.get(),0);
		end = GappedFactorArrayConstIterator<sse4,signature_type>(&comp,A.get(),c);
	
                return c;
        }

        template<typename pattern_type, bool scores, typename updater>
	void match(
		signature_type const s_a,
		signature_type const s_b,
		u_int32_t const seedoffset,
	        AutoTextArrayGappedComparator<sse4> const & ocomp,
	        u_int32_t const ooffset,
	        RestWordBuffer<sse4> & RWB,
	        bool const inverted,
	        pattern_type const & pattern,
	        u_int64_t const fileid,
	        //
	        typename updater::info_type & info,
	        //
	        float const epsilon
		) const
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

		signature_type const prefix = s_a >> SC.s0shift(getSampleBits());
		
		u_int32_t const low = lookup[2*prefix + 0];
		u_int32_t const high = lookup[2*prefix + 1];
		
		std::pair< GappedFactorArrayConstIterator<sse4,signature_type>, GappedFactorArrayConstIterator<sse4,signature_type> > const eq = 
			std::equal_range (begin+low, begin+high, s_a);

		for ( u_int64_t p = eq.first.i; p != eq.second.i; ++p )
		{
			u_int64_t const fragpos = A[p];
			
			if ( fragpos >= seedoffset )
			{
				// seed position
				u_int64_t const seedpos = fragpos-seedoffset;
				
				// position of rest of seed
				u_int64_t const opos = seedpos + ooffset;
				
				// compute number of errors in seed
				unsigned int const seedk = toollib::PopCount<sse4>::diffcountpair(static_cast<signature_type>(s_b),static_cast<signature_type>(ocomp(opos)));

        			if( seedk <= opts.seedkmax )
        			{
                                        if ( seedpos >= matchoffset )
                                        {
                                                // suspected pattern position
                                                u_int32_t const pos = seedpos - matchoffset;
                   
                                                if ( RV.isPositionValid(pos,RWB.patl) && ATA.isDontCareFree(pos,RWB.patl) )
                                                {
                                                        u_int32_t const restpos = seedpos + textrestoffset;
                                                        
                                                        unsigned int const restk = RestMatch<sse4>::computeDistance ( restwordarray, RWB.fullrestwords, RWB.fracrestsyms, ATA, restpos );
                                                        unsigned int const totalk = seedk + restk;

                                                        if ( totalk <= opts.totalkmax )
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
	}
};

template<bool sse4, typename signature_type>
struct GappedFactorArraySet
{
	std::auto_ptr < GappedFactorArray<sse4,signature_type> > FA0;
        std::auto_ptr < GappedFactorArray<sse4,signature_type> > FA1;
        std::auto_ptr < GappedFactorArray<sse4,signature_type> > FA2;

        AutoTextArrayGappedComparator<sse4> const ocomp23;
	AutoTextArrayGappedComparator<sse4> const ocomp03;
        AutoTextArrayGappedComparator<sse4> const ocomp13;
        
        RealOptions const & opts;
        SignatureConstruction<signature_type> const & SC;
        
        GappedFactorArraySet(
                AutoTextArray<sse4> const & ATA, 
                RangeVector<sse4> const & RV,
                RealOptions const & ropts,
                u_int64_t const uniqueinfomemory,
                SignatureConstruction<signature_type> const & rSC,
                Scoring const & scoring
        )
        :
        	FA0(new GappedFactorArray<sse4,signature_type> (ATA,RV,ropts,0/*gapsize*/,uniqueinfomemory,rSC,scoring)),
        	FA1(new GappedFactorArray<sse4,signature_type> (ATA,RV,ropts,1/*gapsize*/,uniqueinfomemory,rSC,scoring)),
        	FA2(new GappedFactorArray<sse4,signature_type> (ATA,RV,ropts,2/*gapsize*/,uniqueinfomemory,rSC,scoring)),
        	ocomp23(ATA,ropts.seedl/4,0,ropts.seedl/4),
        	ocomp03(ATA,ropts.seedl/4,2*(ropts.seedl/4),ropts.seedl/4),
        	ocomp13(ATA,ropts.seedl/4,ropts.seedl/4,ropts.seedl/4),
        	opts(ropts),
        	SC(rSC)
	{}
	
	bool getNextBlock()
	{
	        toollib::RealTimeClock rtc; rtc.start();
	        std::cerr << "Getting next block...";

	        std::cerr << "(";
		u_int64_t const fa0 = FA0->getNextBlock();
		std::cerr << ")";

	        std::cerr << "(";
		u_int64_t const fa1 = FA1->getNextBlock();
		std::cerr << ")";

	        std::cerr << "(";
		u_int64_t const fa2 = FA2->getNextBlock();
		std::cerr << ")";

		std::cerr << "done, time " << rtc.getElapsedSeconds() << " size " << fa0 << std::endl;
		return fa0 || fa1 || fa2;
	}
	
	void reset()
	{
		FA0.reset();
		FA1.reset();
		FA2.reset();
	}
	
	template<typename pattern_type, bool scores, typename info_type>
	bool match(
		pattern_type const & pattern,
		u_int64_t const fi,
		RestWordBuffer<sse4> & RWB,
		info_type & info,
		u_int64_t & handled
		)
	{
		unsigned int const patl = pattern.getPatternLength();

		if ( patl < static_cast<unsigned int>(opts.seedl) )
		{
			std::cerr << "Skipping pattern " << toollib::remapString(std::string(pattern.mapped,pattern.mapped+patl)) << "shorter than seed length." << std::endl;
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

		if ( SC.signatureMapped(pattern.mapped,&m[0]) )
		{
			// fill straight rest word array (pattern bits beyond seed)
			RWB.setupStraight(pattern.mapped);

			// straight pattern fragment sets
			signature_type const straight[] = { 
				SC.s0(m[0],m[1]), // SC.s0(m[0],m[1]), 
				SC.s1(m[0],m[2]), 
				SC.s2(m[0],m[3]), 
				SC.s3(m[1],m[2]), 
				SC.s4(m[1],m[3]), 
				SC.s5(m[2],m[3]) // SC.s5(m[2],m[3]) 
			};
			
			// std::cerr << "Seed " << SC.maskToString(m[0],m[1],m[2],m[3]) << std::endl;

			float const epsilon = opts.getFilterValue(patl);

			FA0->match<pattern_type,scores,UpdateUniqueInfo<scores> >(straight[0],straight[5], 0 /* seedoffset */, ocomp23, 2*(opts.seedl/4) /* seed rest offset */, RWB, false,pattern,fi,info,epsilon);

                        bool const uni0s = (info.getState() == UniqueMatchInfoBase::Straight) && (info.getErrors() == 0);
                        
                        if ( (! uni0s) || scores )
                        {
				FA0->match<pattern_type,scores,UpdateUniqueInfo<scores> >(straight[3],straight[2], opts.seedl/4 /* seedoffset */, ocomp03, 0 /* seed rest offset */, RWB, false,pattern,fi,info,epsilon);
				FA0->match<pattern_type,scores,UpdateUniqueInfo<scores> >(straight[5],straight[0], 2*(opts.seedl/4) /* seedoffset */, ocomp23, 0 /* seed rest offset */, RWB, false,pattern,fi,info,epsilon);							
				FA1->match<pattern_type,scores,UpdateUniqueInfo<scores> >(straight[1],straight[4], 0 /* seedoffset */, ocomp13, 1*(opts.seedl/4) /* seed rest offset */, RWB, false,pattern,fi,info,epsilon);
				FA1->match<pattern_type,scores,UpdateUniqueInfo<scores> >(straight[4],straight[1], opts.seedl/4 /* seedoffset */, ocomp13, 0 /* seed rest offset */, RWB, false,pattern,fi,info,epsilon);
				FA2->match<pattern_type,scores,UpdateUniqueInfo<scores> >(straight[2],straight[3], 0 /* seedoffset */, ocomp23, 1*(opts.seedl/4) /* seed rest offset */, RWB, false,pattern,fi,info,epsilon);
			}

			// compute reverse transpose signature seed
			u_int32_t im[4] = { 0,0,0,0 };

			SC.reverseMappedSignature ( pattern.mapped, &im[0] );
			// reverse pattern fragment sets
			signature_type const reverse[] = { 
				SC.s0(im[0],im[1]), SC.s1(im[0],im[2]), 
				SC.s2(im[0],im[3]), SC.s3(im[1],im[2]),
				SC.s4(im[1],im[3]), SC.s5(im[2],im[3]) };

			// fill reverse rest word array (pattern bits beyond seed)
			RWB.setupReverse(pattern.mapped);

			FA0->match<pattern_type,scores,UpdateUniqueInfo<scores> >(reverse[0],reverse[5], 0 /* seedoffset */, ocomp23, 2*(opts.seedl/4) /* seed rest offset */, RWB, true,pattern,fi,info,epsilon);

                        bool const uni0r = (info.getState() == UniqueMatchInfoBase::Reverse) && (info.getErrors() == 0);

                        if ( (! uni0r) || scores )
                        {
				FA0->match<pattern_type,scores,UpdateUniqueInfo<scores> >(reverse[3],reverse[2], opts.seedl/4 /* seedoffset */, ocomp03, 0 /* seed rest offset */, RWB, true,pattern,fi,info,epsilon);							
				FA0->match<pattern_type,scores,UpdateUniqueInfo<scores> >(reverse[5],reverse[0], 2*(opts.seedl/4) /* seedoffset */, ocomp23, 0 /* seed rest offset */, RWB, true,pattern,fi,info,epsilon);							
				FA1->match<pattern_type,scores,UpdateUniqueInfo<scores> >(reverse[1],reverse[4], 0 /* seedoffset */, ocomp13, 1*(opts.seedl/4) /* seed rest offset */, RWB, true,pattern,fi,info,epsilon);
				FA1->match<pattern_type,scores,UpdateUniqueInfo<scores> >(reverse[4],reverse[1], opts.seedl/4 /* seedoffset */, ocomp13, 0 /* seed rest offset */, RWB, true,pattern,fi,info,epsilon);
				FA2->match<pattern_type,scores,UpdateUniqueInfo<scores> >(reverse[2],reverse[3], 0 /* seedoffset */, ocomp23, 1*(opts.seedl/4) /* seed rest offset */, RWB, true,pattern,fi,info,epsilon);
			}

			handled++;

			return true;
		}        
	
		return false;
	}
};

template<bool sse4, typename signature_type, typename reader_type, bool scores>
void EnumerateUniqueMatches<sse4,signature_type,reader_type,scores>::doMatching(RealOptions const & opts)
{
        typedef typename reader_type::pattern_type pattern_type;
        
#if defined(_OPENMP)
        unsigned int const numthreads = omp_get_max_threads();
        std::cerr << "Starting matching for unique best hits using " << numthreads << " matching threads." << std::endl;
#else
        unsigned int const numthreads = 1;
#endif
        
	u_int64_t const numpat = reader_type::countPatterns(opts.patternfilename);
	
        // allocate memory for matches
        AutoArray < UniqueMatchInfo<scores> > uniqueinfo(numpat);
        
        // get filenames of input files
        std::vector< std::string > filenames;
        std::string const & textfilesuffix = ".fa";
        getFileList(opts.textfilename,filenames,textfilesuffix);
        
        SignatureConstruction<signature_type> const SC(opts.seedl,opts.nu);

        u_int64_t unique = 0;
        
        // double const cps = CLOCKS_PER_SEC;

        int const qualityOffset = opts.qualityOffset ? opts.qualityOffset : reader_type::getOffset(opts.patternfilename);
        
        if ( ! qualityOffset )
                throw std::runtime_error("Unable to automatically detect FastQ quality format.");
        
        Scoring const scoring(opts.similarity, opts.gc, opts.trans, opts.err, opts.gcmut_bias);
        RangeSet RS;

        for ( unsigned int fi = 0; fi < filenames.size(); ++fi )
        {
                bool const lastfile = (fi+1)==filenames.size();

                std::cerr << "Processing file " << filenames[fi] << (lastfile?" (last processed file)":"")<< std::endl;

                std::vector < std::pair < std::string, u_int64_t > > ranges;
                
                std::auto_ptr < AutoTextArray<sse4> > AATA = getText<sse4>(filenames[fi],ranges);
                AutoTextArray<sse4> const & ATA = *(AATA.get());
                RangeVector<sse4> RV(ranges);                
                RS.addRange(ranges);

                u_int64_t const filesize = ATA.getN();
                std::cerr << "File size = " << filesize << std::endl;
                
                if ( ATA.getN() < static_cast<unsigned int>(opts.seedl) )
                {
                        std::cerr << "File " << filenames[fi] << " is too small for seed length, skipping it." << std::endl;
                        continue;
                }

                if ( ranges.size() > UniqueMatchInfo<scores>::getMaxFragmentsPerFile() )
                {
                        std::cerr << "Number of fragments " << ranges.size() << " in file is larger than limit " << UniqueMatchInfo<scores>::getMaxFragmentsPerFile() << " we can handle, skipping it." << std::endl;
                        continue;
                }

		// #define GAPPED_FACTOR_ARRAY
                #if defined(GAPPED_FACTOR_ARRAY)
                toollib::Lock handledlock;
                
                if ( opts.seedl > 32 )
                	throw std::runtime_error("Cannot handle seeds of length > 32 with gapped factor array.");

	        GappedFactorArraySet<sse4,signature_type> GFAS(ATA,RV,opts,uniqueinfo.size(),SC,scoring);
	        
	        clock_t matchclocks = 0;
                
                while ( (GFAS.getNextBlock()) )
                {
                        reader_type patfile(opts.patternfilename, qualityOffset);
                        typename reader_type::stream_data_type ad(patfile,FastFileDecoderBase::default_blocksize,std::max(4u,3*numthreads));
                        typename reader_type::stream_reader_type ar(ad);
                        u_int64_t totalhandled = 0;

#if defined(_OPENMP)
#pragma omp parallel
#endif
			{                
                                typename reader_type::block_type * patternblock = 0;
                                RestWordBuffer<sse4> RWB(opts.seedl);

                                while ( (patternblock = ar.getBlock()) )
                                {
                                        u_int64_t handled = 0;
                                
                                        clock_t const bef = clock();
                                        for ( u_int64_t z = 0; z < patternblock->blocksize; ++z )
                                        {
                                                pattern_type const & pattern = patternblock->getPattern(z);
                                                GFAS.match< pattern_type, scores, UniqueMatchInfo<scores> >(pattern,fi,RWB,
                                                	uniqueinfo[pattern.getPatID()],handled
                                                	);
                                        }
                                        clock_t const aft = clock();
                                       
                                       handledlock.lock();
                                       totalhandled += handled;
                                       matchclocks += (aft-bef);
                                       handledlock.unlock(); 
                                       std::cerr << "\r                                  \r" << static_cast<double>(totalhandled)/static_cast<double>(numpat);

                                       ar.returnBlock(patternblock);                                       
				}
			}
                }
                std::cerr << "\r                                  \r" << 1 << std::endl;
                GFAS.reset();
                std::cerr << "Total matching time " << 
                	static_cast<double>(matchclocks)/CLOCKS_PER_SEC << std::endl;

                #else

                // open text file
                MapTextFile<signature_type,sse4> MTF(ATA, opts.seedl, opts.nu);

                // pointer type in masks
                typedef u_int32_t ptr_type;
                
                u_int64_t const textmemory = ATA.size();
                u_int64_t const rangevectormemory = RV.size();
                u_int64_t const lookupmemory = 2 * getNumLists() * getHistSize() * sizeof(size_t);
                u_int64_t const trlmemory = textmemory + rangevectormemory + lookupmemory;
                u_int64_t const uniqueinfomemory = uniqueinfo.size();
                u_int64_t const nonlistmemory = trlmemory + uniqueinfomemory;
                
                if ( nonlistmemory > opts.usemem )
                {
                        std::cerr << "Insufficient memory." << std::endl;
                        throw std::bad_alloc();
                }
                
                u_int64_t const restmemory = opts.usemem - nonlistmemory;
                u_int64_t const elementmemory = ((getNumLists()>>1) * sizeof(BaseMask<signature_type,ptr_type>) + ((getNumLists()>>1) + getRadixSortTemp()) * sizeof(Mask<signature_type,ptr_type>));
                u_int64_t const n_list_max = restmemory / elementmemory;
               
                std::cerr << "memory         : " << std::endl;
                std::cerr << " -to use       : " << opts.usemem << std::endl;
                std::cerr << " -text         : " << textmemory << std::endl;
                std::cerr << " -range vector : " << rangevectormemory << std::endl;
                std::cerr << " -uniqueinfo   : " << uniqueinfomemory << std::endl;
                std::cerr << " -lookupmemory : " << lookupmemory << std::endl;
                std::cerr << "sum            : " << nonlistmemory << std::endl;
                std::cerr << "rest           : " << restmemory << std::endl;
                std::cerr << " -elementmemory: " << elementmemory << std::endl;

                // number of element per list        
                std::cerr << "We have space for " << n_list_max << " list elements." << std::endl;
        
                if ( ! n_list_max )
                        throw std::bad_alloc();

                u_int64_t const expblocks = (filesize - opts.seedl + 1 + (n_list_max-1)) / n_list_max;
                std::cerr << "Expected number of blocks = " << expblocks << std::endl;
                u_int64_t const n_list = (filesize + (expblocks-1)) / expblocks;
                std::cerr << "Using n_list = " << n_list << std::endl;

                unsigned int masks = 0;

                UniqueMatcher<signature_type,sse4,ptr_type,pattern_type,scores> UM(n_list,opts,SC,MTF,scoring,ATA,RV);
                
                /**
                 * read blocks of input file
                 **/
                while ( (masks = UM.readNextBlock()) > 0 )
                {
#if 0
                        bool const lastblock = (!UM.have_next);
                        bool const lastfilelastblock = lastfile && lastblock;
#endif
                                
                        reader_type patfile(opts.patternfilename, qualityOffset);
                        typename reader_type::stream_data_type ad(patfile,FastFileDecoderBase::default_blocksize,std::max(4u,3*numthreads));
                        typename reader_type::stream_reader_type ar(ad);

                        unsigned int totalhandled = 0;

                        toollib::Lock stderrlock;

#if defined(_OPENMP)                        
                        #pragma omp parallel
#endif
                        {
                                typename reader_type::block_type * patternblock = 0;
                                RestWordBuffer<sse4> RWB(opts.seedl);
                                
                                while ( (patternblock = ar.getBlock()) )
                                {
                                        u_int64_t handled = 0;

                                        for ( u_int64_t z = 0; z < patternblock->blocksize; ++z )
                                        {
                                                pattern_type const & pattern = patternblock->getPattern(z);
                                                UM.match(pattern, uniqueinfo[pattern.getPatID()], fi, RWB, handled);
                                        }
                                        
                                        ar.returnBlock(patternblock);

#if defined(_OPENMP)
                                        #pragma omp critical
#endif
                                        {
                                                totalhandled += handled;
                                                std::cerr << "\r                                                              \r" << (static_cast<double>(totalhandled)/numpat) << "\t" << ar.getFillState() << std::flush;
                                        }
                                }
                        }
                        std::cerr << std::endl;
                }
                std::cerr << "All done." << std::endl;
                #endif
        }

        if ( opts.gaps )
        {
		std::map<unsigned int, GapInfo> gapinfos;

		for ( unsigned int fi = 0; fi < filenames.size(); ++fi )
		{
			bool const lastfile = (fi+1)==filenames.size();

			std::cerr << "Processing file " << filenames[fi] << (lastfile?" (last processed file)":"")<< std::endl;

			std::vector < std::pair < std::string, u_int64_t > > ranges;
			
			std::auto_ptr < AutoTextArray<sse4> > AATA = getText<sse4>(filenames[fi],ranges);
			AutoTextArray<sse4> const & ATA = *(AATA.get());
			RangeVector<sse4> RV(ranges);                
			// RS.addRange(ranges);

			u_int64_t const filesize = ATA.getN();
			std::cerr << "File size = " << filesize << std::endl;
			
			if ( ATA.getN() < static_cast<unsigned int>(opts.seedl) )
			{
				std::cerr << "File " << filenames[fi] << " is too small for seed length, skipping it." << std::endl;
				continue;
			}

			if ( ranges.size() > UniqueMatchInfo<scores>::getMaxFragmentsPerFile() )
			{
				std::cerr << "Number of fragments " << ranges.size() << " in file is larger than limit " << UniqueMatchInfo<scores>::getMaxFragmentsPerFile() << " we can handle, skipping it." << std::endl;
				continue;
			}
			
			// open text file
			MapTextFile<signature_type,sse4> MTF(ATA, opts.seedl, opts.nu);

			// pointer type in masks
			typedef u_int32_t ptr_type;
			
			u_int64_t const textmemory = ATA.size();
			u_int64_t const rangevectormemory = RV.size();
			u_int64_t const lookupmemory = 2 * getNumLists() * getHistSize() * sizeof(size_t);
			u_int64_t const trlmemory = textmemory + rangevectormemory + lookupmemory;
			u_int64_t const uniqueinfomemory = uniqueinfo.size();
			u_int64_t const nonlistmemory = trlmemory + uniqueinfomemory;
			
			if ( nonlistmemory > opts.usemem )
			{
				std::cerr << "Insufficient memory." << std::endl;
				throw std::bad_alloc();
			}
			
			u_int64_t const restmemory = opts.usemem - nonlistmemory;
			u_int64_t const elementmemory = ((getNumLists()>>1) * sizeof(BaseMask<signature_type,ptr_type>) + ((getNumLists()>>1) + getRadixSortTemp()) * sizeof(Mask<signature_type,ptr_type>));
			u_int64_t const n_list_max = restmemory / elementmemory;
		       
			std::cerr << "memory         : " << std::endl;
			std::cerr << " -to use       : " << opts.usemem << std::endl;
			std::cerr << " -text         : " << textmemory << std::endl;
			std::cerr << " -range vector : " << rangevectormemory << std::endl;
			std::cerr << " -uniqueinfo   : " << uniqueinfomemory << std::endl;
			std::cerr << " -lookupmemory : " << lookupmemory << std::endl;
			std::cerr << "sum            : " << nonlistmemory << std::endl;
			std::cerr << "rest           : " << restmemory << std::endl;
			std::cerr << " -elementmemory: " << elementmemory << std::endl;

			// number of element per list        
			std::cerr << "We have space for " << n_list_max << " list elements." << std::endl;
		
			if ( ! n_list_max )
				throw std::bad_alloc();

			u_int64_t const expblocks = (filesize - opts.seedl + 1 + (n_list_max-1)) / n_list_max;
			std::cerr << "Expected number of blocks = " << expblocks << std::endl;
			u_int64_t const n_list = (filesize + (expblocks-1)) / expblocks;
			std::cerr << "Using n_list = " << n_list << std::endl;

			unsigned int masks = 0;

			UniqueMatcher<signature_type,sse4,ptr_type,pattern_type,scores> UM(n_list,opts,SC,MTF,scoring,ATA,RV);
			
			/**
			 * read blocks of input file
			 **/
			while ( (masks = UM.readNextBlock()) > 0 )
			{
	#if 0
				bool const lastblock = (!UM.have_next);
				bool const lastfilelastblock = lastfile && lastblock;
	#endif
					
				reader_type patfile(opts.patternfilename, qualityOffset);
				typename reader_type::stream_data_type ad(patfile,FastFileDecoderBase::default_blocksize,std::max(4u,3*numthreads));
				typename reader_type::stream_reader_type ar(ad);

				unsigned int totalhandled = 0;

				toollib::Lock stderrlock;

	#if defined(_OPENMP)                        
				#pragma omp parallel
	#endif
				{
					typename reader_type::block_type * patternblock = 0;
					RestWordBuffer<sse4> RWB(opts.seedl);
					
					while ( (patternblock = ar.getBlock()) )
					{
						u_int64_t handled = 0;
					
						for ( u_int64_t z = 0; z < patternblock->blocksize; ++z )
						{
							pattern_type const & pattern = patternblock->getPattern(z);

							UniqueMatchInfo<scores> & info = uniqueinfo[pattern.getPatID()];

							UM.matchGaps(pattern, info, RWB, gapinfos, handled);
						}
							
						
						ar.returnBlock(patternblock);

	#if defined(_OPENMP)
						#pragma omp critical
	#endif
						{
							totalhandled += handled;
							std::cerr << "\r                                                              \r" << (static_cast<double>(totalhandled)/numpat) << "\t" << ar.getFillState() << std::flush;
						}
					}
				}
				std::cerr << std::endl;
			}
			std::cerr << "All done." << std::endl;
		}
	}

        PatternIdReader<reader_type> pir (opts.patternfilename, qualityOffset, numthreads);
        std::pair < typename reader_type::block_type *, typename reader_type::idblock_type * > block;
        
        /**
         * output
         **/
        std::auto_ptr < AsynchronousWriter  > aw;
        
        if ( opts.outputfilename == "-" )
        {
		aw = std::auto_ptr < AsynchronousWriter  >( new AsynchronousWriter(STDOUT_FILENO,16) );
	}
	else
	{
		aw = std::auto_ptr < AsynchronousWriter  >( new AsynchronousWriter(opts.outputfilename,16) );
	}
        
        while (  pir.getBlock(block) )
        {
                for ( u_int64_t i = 0; i < block.first->blocksize; ++i )
                {
                        pattern_type const & pattern = block.first->getPattern(i);
                        std::string const & id = block.second->ids[i];
                        UniqueMatchInfo<scores> & info = uniqueinfo[pattern.getPatID()];
                        
                        std::ostringstream ostr;
                        printMatchUnlocked<scores>(info,ostr,pattern,id,pattern.getPatternLength(),unique,RS);
                        std::string const os = ostr.str();
                        aw->write(os.begin(), os.end());

                        #if defined(GAPPED_MATCHING)
                        if ( (info.getState() == UniqueMatchInfoBase::Gapped) && ( gapinfos.find(pattern.getPatID() )  != gapinfos.end() ) )
                        {
                                GapInfo const & gapinfo = gapinfos.find(pattern.getPatID())->second;
                        
                                std::cerr << "Found gapped match, score " << info.getScore()
                                        << " position=" << info.getPosition()
                                        << " seed length=" << opts.seedl
                                        << " MINgap=" << gapinfo.MINgap
                                        << " where=" << gapinfo.where
                                        << " start=" << gapinfo.start
                                        << " gap_pos=" << gapinfo.gap_pos
                                        << std::endl;
                        }
                        #endif
                }

                pir.returnBlock(block);                
        }

        std::cerr << "unique: " << unique << std::endl;
}
#endif
