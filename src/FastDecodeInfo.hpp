#if ! defined(FASTDECODEINFO_HPP)
#define FASTDECODEINFO_HPP

#include "types.hpp"
#include <ostream>

struct FastDecodeInfo
{
	unsigned int patlen;
	u_int64_t idbase;
	u_int64_t idnbase;
	u_int64_t patpos;
	u_int64_t patid;
	u_int64_t patnum;
	u_int64_t patnpos;
	u_int64_t patnid;
	u_int64_t patnnum;
	
	u_int64_t getNumPatterns() const
	{
		return patnum + patnnum;
	}
};

std::ostream & operator<< (std::ostream & out, FastDecodeInfo const & fdi);
#endif
