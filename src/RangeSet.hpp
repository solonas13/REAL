#if ! defined(RANGESET_HPP)
#define RANGESET_HPP

#include <vector>
#include <map>
#include <string>
#include <sys/types.h>

struct RangeSet
{
        std::vector < std::vector <  std::pair < std::string, u_int64_t > > > ranges;
        
        RangeSet() {}
        
        void addRange(std::vector < std::pair < std::string, u_int64_t > > const & rranges)
        {
                ranges.push_back ( rranges );
        }
        
        std::string getId(unsigned int const fileid, u_int64_t pos)
        {
                std::vector <  std::pair < std::string, u_int64_t > > const & lranges = ranges[fileid];
                
                unsigned int i = 0;
                while ( lranges[i+1].second <= pos )
                        ++i;
                        
                return lranges[i].first;
        }
};
#endif

