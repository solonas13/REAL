#if ! defined(MATCHERBASE_HPP)
#define MATCHERBASE_HPP

#include "ListSetBlockReader.hpp"
#include "AutoTextArray.hpp"
#include "RangeVector.hpp"
#include "Scoring.hpp"

template<
        typename signature_type, bool sse4,
        typename ptr_type, typename pattern_type,
        bool scores
>
struct MatcherBase : public ListSetBlockReader<signature_type,sse4>
{
        typedef ListSetBlockReader<signature_type,sse4> base;

        Scoring const & scoring;
        AutoTextArray<sse4> const & ATA;
        RangeVector<sse4> const & RV;

        MatcherBase(
                u_int64_t const n_list, 
                RealOptions const & ropts, 
                SignatureConstruction<signature_type> const & rSC,
                MapTextFile<signature_type,sse4> & rMTF,
                Scoring const & rscoring,
                AutoTextArray<sse4> const & rATA,
                RangeVector<sse4> const & rRV
        ) : ListSetBlockReader<signature_type,sse4>(n_list,ropts,rSC,rMTF), 
            scoring(rscoring), ATA(rATA), RV(rRV)
        {
        
        }

};
#endif
