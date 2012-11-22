#if ! defined(LISTSET_HPP)
#define LISTSET_HPP

#include "AutoArray.hpp"
#include "SignatureConstruction.hpp"
#include "Mask.hpp"
#include "RealOptions.hpp"
#include "u_sort.hpp"
#include "getHistSize.hpp"
#include "MapTextFile.hpp"
#include "RealTimeClock.hpp"

template<typename signature_type>
struct ListSet
{
        static unsigned int const numfulllists = 3;
        static unsigned int const numbaselists = 3;

        // allocate lists
        AutoArray< Mask<signature_type, unsigned int> > Afulllists[numfulllists];
        AutoArray< BaseMask<signature_type, unsigned int> > Abaselists[numbaselists];

        Mask<signature_type, unsigned int> const * fulllists[numfulllists];
        BaseMask<signature_type, unsigned int> const * baselists[numfulllists];
        
        u_int64_t n_list;
        RealOptions const & opts;
        SignatureConstruction<signature_type> const & SC;

        AutoArray< size_t > Alookup[numfulllists + numbaselists];
        size_t * lookup[numfulllists + numbaselists];
        
        ListSet(u_int64_t const rn_list, RealOptions const & ropts, SignatureConstruction<signature_type> const & rSC) : n_list(rn_list), opts(ropts), SC(rSC)
        {
                for ( unsigned int i = 0; i < numfulllists; ++i )
                        Afulllists[i] = AutoArray< Mask<signature_type, unsigned int> > ( n_list, false);
                for ( unsigned int i = 0; i < numbaselists; ++i )
                        Abaselists[i] = AutoArray< BaseMask<signature_type, unsigned int> > ( n_list, false);
        }
        
        void sort(u_int64_t const masks)
        {
                for ( unsigned int i = 0; i < numfulllists; ++i )
                        MaskSort< signature_type, unsigned int>::sort(Afulllists[i], Abaselists[numbaselists-i-1], masks, opts.sort_threads);

                for ( unsigned int i = 0; i < numfulllists; ++i )
                        fulllists[i] = Afulllists[i].get();
                for ( unsigned int i = 0; i < numbaselists; ++i )
                        baselists[i] = Abaselists[i].get();

                // lookup table parameters 
                unsigned int const histsize = getHistSize();

                Alookup[0] = getLookupTable(fulllists[0], masks, SC.s0shift(getSampleBits()), histsize);
                Alookup[1] = getLookupTable(fulllists[1], masks, SC.s1shift(getSampleBits()), histsize);
                Alookup[2] = getLookupTable(fulllists[2], masks, SC.s2shift(getSampleBits()), histsize);
                Alookup[3] = getLookupTable(baselists[0], masks, SC.s3shift(getSampleBits()), histsize);
                Alookup[4] = getLookupTable(baselists[1], masks, SC.s4shift(getSampleBits()), histsize);
                Alookup[5] = getLookupTable(baselists[2], masks, SC.s5shift(getSampleBits()), histsize);
                
                for ( unsigned int i = 0; i < sizeof(lookup)/sizeof(lookup[0]); ++i )
                        lookup[i] = Alookup[i].get();
        }
};
#endif
