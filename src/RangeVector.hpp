#if ! defined(RANGEVECTOR_HPP)
#define RANGEVECTOR_HPP

#include <vector>
#include <string>
#include "types.hpp"
#include "AutoArray.hpp"
#include "ERank222B.hpp"

template<bool sse4>
struct RangeVector
{
        std::vector < std::pair < std::string, u_int64_t > > const & ranges;
        AutoArray<u_int64_t> AB;
        u_int64_t * const B;
        Rank::ERank222B<sse4> R;
        
        u_int64_t size() const
        {
                 return AB.size() + sizeof(u_int64_t *) + R.size();
        }
        
        static AutoArray<u_int64_t> fillRange(
                std::vector < std::pair < std::string, u_int64_t > > const & ranges
                )
        {
                AutoArray<u_int64_t> AB ( (ranges.at(ranges.size()-1).second+1+63)/64 );

                for ( u_int64_t i = 0; i < (ranges.at(ranges.size()-1).second+1+63)/64; ++i )
                        AB[i] = 0;

                Rank::BitWriter8 W(AB.get());
                
                for ( unsigned int j = 0; j+1 < ranges.size(); ++j )
                {
                        u_int64_t len = ranges[j+1].second-ranges[j].second;
                        W.writeBit(true);
                        for ( u_int64_t i = 1; i < len; ++i )
                                W.writeBit(false);
                }
                W.writeBit(true);
                
                W.flush();
                
                return AB;
        }
        
        RangeVector(std::vector < std::pair < std::string, u_int64_t > > const & rranges)
        : ranges(rranges), AB ( fillRange(ranges) ), B ( AB.get() ),
          R(B, ((ranges.at(ranges.size()-1).second+1+63)/64)*64 )
        {
                #if 0
                for ( u_int64_t i = 0; i <= ranges.back().second; ++i )
                        std::cerr << Rank::getBit(B,i) << "(" << R.rank1(i)-1 << ")";
                std::cerr << std::endl;
                #endif
        }
        
        unsigned int positionToRange(u_int64_t const pos) const
        {
                return R.rank1(pos)-1;
        }
        bool isPositionValid(u_int64_t const pos, unsigned int const patl) const
        {
                unsigned int const range = positionToRange(pos);
                bool const ok = (pos + patl) <= ranges[range+1].second;

                #if 0                
                if ( ! ok )
                {
                        std::cerr << "position " << pos << " rejected"
                                << " " << ranges[range].second
                                << " " << ranges[range+1].second
                                << " patl=" << patl
                                << std::endl;
                }
                #endif
                
                return ok;
        }
        std::string const & positionToId(u_int64_t const pos) const
        {
                return ranges[positionToRange(pos)].first;
        }
};
#endif
