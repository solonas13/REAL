#if ! defined(RESTWORDBUFFER_HPP)
#define RESTWORDBUFFER_HPP

#include "AutoArray.hpp"
#include "RestMatch.hpp"

template<bool sse4>
struct RestWordBuffer
{
        AutoArray<u_int64_t> ABstraight;
        u_int64_t * Bstraight;
        AutoArray<u_int64_t> ABreverse;
        u_int64_t * Breverse;

        unsigned int const seedl;
        unsigned int maxpatl;
        unsigned int patl;

        unsigned int numrestwords;
        unsigned int fullrestwords;
        unsigned int fracrestsyms;
        int straighttextrestoffset;
        int reversetextrestoffset;
        unsigned straightmatchoffset;
        unsigned reversematchoffset;

        RestWordBuffer(unsigned int const rseedl)
        : ABstraight(), Bstraight(0), ABreverse(), Breverse(0), seedl(rseedl), maxpatl(0), patl(0)
        {
        
        }
        
        void setup(unsigned int npatl)
        {
                if ( npatl > maxpatl )
                {
                        ABstraight = AutoArray<u_int64_t>(RestMatch<sse4>::getRestWords(npatl,seedl));
                        Bstraight = ABstraight.get();
                        ABreverse = AutoArray<u_int64_t>(RestMatch<sse4>::getRestWords(npatl,seedl));
                        Breverse = ABreverse.get();
                        maxpatl = npatl;
                }
                if ( npatl != patl )
                {
                        patl = npatl;
                        numrestwords = RestMatch<sse4>::getRestWords(patl,seedl);
                        fullrestwords = RestMatch<sse4>::getFullRestWords(patl, seedl);
                        fracrestsyms = RestMatch<sse4>::getFracRestSyms(patl, seedl, fullrestwords);
                        straighttextrestoffset = RestMatch<sse4>::getTextRestOffset(patl, seedl, false);
                        reversetextrestoffset = RestMatch<sse4>::getTextRestOffset(patl, seedl, true);
                        straightmatchoffset = RestMatch<sse4>::getMatchOffset(patl, seedl, false);
                        reversematchoffset = RestMatch<sse4>::getMatchOffset(patl, seedl, true);                        
                }
        }
        
        template<typename iterator_type>
        void setupStraight(iterator_type mapped)
        {
                // fill straight rest word array (pattern bits beyond seed)
                RestMatch<sse4>::fillRestWordArrayMapped(
                        mapped, 
                        Bstraight, 
                        seedl, 
                        fullrestwords, 
                        fracrestsyms);
        }
        
        template<typename iterator_type>
        void setupReverse(iterator_type mapped)
        {
                // fill reverse rest word array (pattern bits beyond seed)
                RestMatch<sse4>::fillRestWordArrayReverseMapped(
                        mapped, 
                        Breverse, 
                        patl, 
                        fullrestwords, 
                        fracrestsyms);
        }
};
#endif
