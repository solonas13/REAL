#if ! defined(FASTIDDECODER_HPP)
#define FASTIDDECODER_HPP

#include "types.hpp"
#include "FastReaderBase.hpp"
#include "FastIDBlock.hpp"

struct FastIDDecoder
{
	typedef FastIDBlock block_type;

	u_int64_t const numpat;
	u_int64_t decoded;
	FastReaderBase FRB;
	unsigned int blocksize;
	u_int64_t idbase;
	
	FastIDDecoder(
		std::string const & filename,
		u_int64_t const rnumpat,
		u_int64_t offset,
		unsigned int rblocksize,
		u_int64_t ridbase)
	: numpat(rnumpat), decoded(0), FRB(filename,16,1024,offset), blocksize(rblocksize), idbase(ridbase)
	{
	
	}
	
	bool getNextId(std::string & s)
	{
		if ( decoded == numpat )
		{
			return false;
		}
		if ( !FRB.getString(s) )
		{
			return false;
		}
		return true;
	}
	
	u_int64_t fillPatternBlock(block_type & block, unsigned int const rblocksize)
	{
		assert ( blocksize == rblocksize );
		block.setup(blocksize);
		block.idbase = idbase + decoded;

		while ( 
			(decoded < numpat)
			&&
			(block.blocksize < block.numid) 
			&& 
			getNextId(block.ids[block.blocksize])
		)
		{
			decoded++;
			block.blocksize++;
		}
			
		return block.blocksize;
	}
};
#endif
