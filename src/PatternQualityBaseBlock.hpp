#if ! defined(PATTERNQUALITYBASEBLOCK_HPP)
#define PATTERNQUALITYBASEBLOCK_HPP

#include "AutoArray.hpp"
#include "Pattern.hpp"

struct PatternQualityBaseBlock
{
	unsigned int blockid;
	AutoArray<char> mappedBlock;
	AutoArray<char> transposedBlock;
	AutoArray<char> quality;
	unsigned int numpat;
	AutoArray<PatternQualityBase> Apatterns;
	PatternQualityBase * patterns;
	unsigned int blocksize;
	unsigned int patlen;

	PatternQualityBaseBlock()
	: blockid(0), numpat(0), patterns(0), blocksize(0), patlen(0)
	{
	}
	
	PatternQualityBase const & getPattern(u_int64_t const i) const
	{
		return patterns[i];
	}
	
	void setup(
		AutoArray<char> rmappedBlock, 
		AutoArray<char> rtransposedBlock, 
		AutoArray<char> rquality,
		unsigned int const rnumpat,
		unsigned int const rpatlen	
		)
	{
		mappedBlock = rmappedBlock;
		transposedBlock = rtransposedBlock;
		quality = rquality;
		numpat = rnumpat;
		Apatterns = AutoArray<PatternQualityBase>(numpat,false);
		patterns = Apatterns.get();
		blocksize = 0;
		patlen = rpatlen;
	}
};
#endif
