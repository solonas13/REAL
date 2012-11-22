#if ! defined(SPACETABLE_HPP)
#define SPACETABLE_HPP

#include "AutoArray.hpp"

struct SpaceTable
{
        protected:
        AutoArray<bool> spacetable;
        AutoArray<bool> nospacetable;

        SpaceTable() : spacetable(256), nospacetable(256)
        {
	        for ( unsigned int i = 0; i < 256; ++i )
	        {
	                spacetable[i] = isspace(i);
	                nospacetable[i] = !spacetable[i];
                }        
        }
};
#endif
