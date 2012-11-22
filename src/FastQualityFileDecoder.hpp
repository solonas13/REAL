#if ! defined(FASTQUALITYFILEDECODER_HPP)
#define FASTQUALITYFILEDECODER_HPP

#include "FastQualityDecoder.hpp"
#include "FastFileDecoderState.hpp"
#include "FastIDBlock.hpp"
#include "FastQualityIDFileDecoder.hpp"
#include "FastFileDecoderBase.hpp"

struct FastQualityFileDecoder : public FastFileDecoderBase
{
        typedef FastQualityFileDecoder reader_type;
        typedef PatternQualityBase pattern_type;
        typedef PatternQualityBaseBlock block_type;
        typedef FastIDBlock idblock_type;
        typedef FastQualityIDFileDecoder idfile_type;
        typedef AsynchronousFastReaderData<reader_type> stream_data_type;
        typedef AsynchronousFastIdData<reader_type> stream_id_type;
        typedef AsynchronousStreamReader<stream_data_type> stream_reader_type;
        typedef AsynchronousStreamReader<stream_id_type> stream_idreader_type;

	FastQualityDecoder decoder;
	FastFileDecoderState state;
	std::auto_ptr < FastQualitySubDecoder > subdecoder;
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
	        return FastQualityDecoder::countPatterns(filename);
	}

	
	FastQualityFileDecoder(std::string const & filename, int = 0, unsigned int blocksize = default_blocksize)
	: decoder(filename, blocksize), 
	  state(decoder.getPatternLengthVector())
	{
		active = state.getNext(patlen,dc);
		activate();
	}

	u_int64_t fillPatternBlock(PatternQualityBaseBlock & block, unsigned int const blocksize)
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
