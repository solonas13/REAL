#if ! defined(FASTIDBLOCK_HPP)
#define FASTIDBLOCK_HPP

#include "types.hpp"
#include "AutoArray.hpp"
#include <string>

struct FastIDBlock
{
	unsigned int blockid;
	u_int64_t numid;
	u_int64_t blocksize;
	AutoArray < std::string > aids;
	std::string * ids;
	u_int64_t idbase;
	
	FastIDBlock ( )
	: numid(0), blocksize(0), aids(numid,false), ids(aids.get()), idbase(0) {}
	
	void setup(u_int64_t rnumid)
	{
		numid = rnumid;
		blocksize = 0;
		aids = AutoArray < std::string >(numid,false);
		ids = aids.get();
		idbase = 0;
	}
};
#endif
