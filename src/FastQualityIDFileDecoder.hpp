#if ! defined(FASTQUALITYIDDECODER_HPP)
#define FASTQUALITYIDDECODER_HPP

#include "FastQualityDecoder.hpp"
#include "FastFileDecoderState.hpp"
#include "FastFileDecoderBase.hpp"

struct FastQualityIDFileDecoder : public FastFileDecoderBase
{
	FastQualityDecoder decoder;
	FastFileDecoderState state;
	std::auto_ptr < FastIDDecoder > subdecoder;
	unsigned int patlen;
	bool dc;
	bool active;
	typedef FastIDBlock block_type;
	
	void activate()
	{
		if ( active )
		{
			if ( dc )
				subdecoder = decoder.getIdNDecoder(patlen);
			else
				subdecoder = decoder.getIdDecoder(patlen);
		}	
	}
	
	FastQualityIDFileDecoder(std::string const & filename, unsigned int const blocksize = default_blocksize)
	: decoder(filename, blocksize),
	  state(decoder.getPatternLengthVector())
	{
		active = state.getNext(patlen,dc);
		activate();
	}

	u_int64_t fillPatternBlock(block_type & block, unsigned int const blocksize)
	{
		while ( active )
		{
		        u_int64_t const c = subdecoder->fillPatternBlock(block,blocksize);
			if ( c )
				return c;

			active = state.getNext(patlen,dc);
			
			activate();
		}
		return false;
	}
};
#endif
