#if ! defined(LISTSETBLOCKREADER_HPP)
#define LISTSETBLOCKREADER_HPP

#include "ListSet.hpp"

template<typename signature_type, bool sse4>
struct ListSetBlockReader : public ListSet<signature_type>
{
        typedef ListSet<signature_type> base;

        bool have_next;
        MapTextFile<signature_type,sse4> & MTF;
        
        ListSetBlockReader(
                u_int64_t const n_list, 
                RealOptions const & ropts, 
                SignatureConstruction<signature_type> const & rSC,
                MapTextFile<signature_type,sse4> & rMTF
        ) : ListSet<signature_type>(n_list,ropts,rSC), have_next(false), MTF(rMTF)
        {
        
        }

        u_int64_t readNextBlock()
        {
                u_int64_t const masks = MTF.readLists(
                        base::Afulllists[0].get(),
                        base::Afulllists[1].get(),
                        base::Afulllists[2].get(), 
                        base::Abaselists[0].get(), 
                        base::Abaselists[1].get(), 
                        base::Abaselists[2].get(), 
                        base::n_list, 
                        have_next);

                std::cerr << "Obtained " << masks << " fragments of size " << base::opts.seedl << std::endl;

                if ( masks )
                {
                        std::cerr << "Sorting fragments...";
                        
                        toollib::RealTimeClock rtc; rtc.start();
                        clock_t bef_sort = clock();

                        base::sort(masks);                        
                        clock_t aft_sort = clock();

                        std::cerr << "done, real time " << rtc.getElapsedSeconds() << " computation time " << static_cast<double>(aft_sort-bef_sort)/CLOCKS_PER_SEC << std::endl;
                }
                
                return masks;
        }
};
#endif
