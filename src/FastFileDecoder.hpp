#if ! defined(FASTFILEDECODER_HPP)
#define FASTFILEDECODER_HPP

#include "FastDecoder.hpp"
#include "FastFileDecoderState.hpp"
#include "AsynchronousReader.hpp"
#include "FastIDBlock.hpp"
#include "FastIDFileDecoder.hpp"
#include "FastFileDecoderBase.hpp"

struct FastFileDecoder : public FastFileDecoderBase
{
        typedef FastFileDecoder reader_type;
        typedef PatternBase pattern_type;
        typedef PatternBaseBlock block_type;
        typedef FastIDBlock idblock_type;
        typedef FastIDFileDecoder idfile_type;
        typedef AsynchronousFastReaderData<reader_type> stream_data_type;
        typedef AsynchronousStreamReader<stream_data_type> stream_reader_type;
        typedef AsynchronousFastIdData<reader_type> stream_id_type;
        typedef AsynchronousStreamReader<stream_id_type> stream_idreader_type;
        
	FastDecoder decoder;
	FastFileDecoderState state;
	std::auto_ptr < FastSubDecoder > subdecoder;
	unsigned int patlen;
	bool dc;
	bool active;
	
	void activate()
	{
		if ( active )
		{
			if ( dc )
				subdecoder = decoder.getNDecoder(patlen);
			else
				subdecoder = decoder.getDecoder(patlen);
		}	
	}
	
	static u_int64_t countPatterns(std::string const & filename)
	{
	        return FastDecoder::countPatterns(filename);
	}
	
	FastFileDecoder(std::string const & filename, int = 0, unsigned int const blocksize = default_blocksize)
	: decoder(filename,blocksize), 
	  state(decoder.getPatternLengthVector())
	{
		active = state.getNext(patlen,dc);
		activate();
	}

	u_int64_t fillPatternBlock(block_type & block, unsigned int const blocksize)
	{
		while ( active )
		{
			subdecoder->setupBlock(block,blocksize);
		
			if ( subdecoder->fillPatternBlock(block,blocksize) )
				return block.blocksize;
			
			active = state.getNext(patlen,dc);
			
			activate();
		}
		return false;
	}
	static int getOffset(std::string const /*inputfile*/)
	{
	        return -1;
	}	

};
#endif
