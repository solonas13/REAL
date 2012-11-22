#if ! defined(PATTERNBASEBLOCK_HPP)
#define PATTERNBASEBLOCK_HPP

#include "AutoArray.hpp"
#include "Pattern.hpp"

struct PatternBaseBlock
{
	unsigned int blockid;
	AutoArray<char> mappedBlock;
	AutoArray<char> transposedBlock;
	unsigned int numpat;
	AutoArray<PatternBase> Apatterns;
	PatternBase * patterns;
	unsigned int blocksize;
	unsigned int patlen;

	PatternBaseBlock()
	: blockid(0), numpat(0), patterns(0), blocksize(0), patlen(0)
	{
	}
	
	PatternBase const & getPattern(u_int64_t const i) const
	{
		return patterns[i];
	}
	
	void setup(
		AutoArray<char> rmappedBlock, 
		AutoArray<char> rtransposedBlock, 
		unsigned int const rnumpat,
		unsigned int const rpatlen	
		)
	{
		mappedBlock = rmappedBlock;
		transposedBlock = rtransposedBlock;
		numpat = rnumpat;
		Apatterns = AutoArray<PatternBase>(numpat,false);
		patterns = Apatterns.get();
		blocksize = 0;
		patlen = rpatlen;
	}
};
#endif
