#if ! defined(FASTAREADER_HPP)
#define FASTAREADER_HPP

#include "FastReaderBase.hpp"
#include "Pattern.hpp"
#include "SpaceTable.hpp"
#include "CharTermTable.hpp"
#include "CharBuffer.hpp"
#include "PatternBlock.hpp"
#include "AsynchronousReader.hpp"
#include "FastIDBlock.hpp"

struct FastAReader : public FastReaderBase, public SpaceTable
{
        typedef FastAReader reader_type;
        typedef Pattern pattern_type;
        typedef PatternBlock<pattern_type> block_type;
        typedef AsynchronousStreamReaderData<reader_type> stream_data_type;
        typedef AsynchronousStreamReader< AsynchronousStreamReaderData<reader_type> > stream_reader_type;

        typedef reader_type idfile_type;
        typedef FastIDBlock idblock_type;
        typedef AsynchronousIdData<reader_type> stream_id_type;
        typedef AsynchronousStreamReader<stream_id_type> stream_idreader_type;

	CharTermTable scanterm;
	CharTermTable newlineterm;

        CharBuffer idbuffer;	
        CharBuffer patbuffer;

	bool foundnextmarker;
	u_int64_t nextid;
	
	FastAReader(std::string const & filename, int = -1 /* qualityOffset*/, unsigned int numbuffers = 16, unsigned int bufsize = 16*1024)
	: FastReaderBase(filename,numbuffers,bufsize),
	  scanterm('>'), newlineterm('\n'),
	  foundnextmarker(false), nextid(0)
	{
		findNextMarker();
	}
	
	void findNextMarker()
	{
		int c;
		
		while ( !scanterm[(c=getNextCharacter())] )
		{
			// std::cerr << "Got " << static_cast<char>(c) << " without term, " << scanterm[c] << std::endl;
		}
		
		foundnextmarker = (c=='>');
	}
	
	u_int64_t countPatterns()
	{
	        u_int64_t count = 0;
	        
	        while ( skipPattern() )
	                ++count;

                return count;
	}

	static u_int64_t getPatternLength(std::string const & filename)
	{
	        FastAReader reader(filename);
	        pattern_type pat;
	        
	        if ( ! reader.getNextPatternUnlocked(pat) )
	                return 0;
                else
                        return pat.patlen;
	}
	
	static u_int64_t countPatterns(std::string const & filename)
	{
	        FastAReader reader(filename);
	        return reader.countPatterns();
	}

        bool skipPattern()
	{
		if ( ! foundnextmarker )
			return false;
			
		foundnextmarker = false;
	
		int c;

		while ( !newlineterm[c=getNextCharacter()] )
		{
                }
		
		if ( c < 0 )
			return false;
			
		while ( !scanterm[(c=getNextCharacter())] )
		{
                }

		foundnextmarker = (c=='>');

		return true;
	}

        bool getNextPatternUnlocked(pattern_type & pattern)
	{
		if ( ! foundnextmarker )
			return false;
			
		foundnextmarker = false;
	
		int c;

		idbuffer.reset();
		while ( !newlineterm[c=getNextCharacter()] )
			idbuffer.bufferPush(c);
		
		if ( c < 0 )
			return false;
			
		idbuffer.assign ( pattern.sid );

		patbuffer.reset();
		while ( !scanterm[(c=getNextCharacter())] )
			if ( nospacetable[c] )
				patbuffer.bufferPush(c);

                patbuffer.assign(pattern.spattern);
		pattern.pattern = pattern.spattern.c_str();
		pattern.patlen = pattern.spattern.size();
		pattern.patid = nextid++;
		
		foundnextmarker = (c=='>');

		return true;
	}

	unsigned int fillPatternBlock(
	        pattern_type * pat,
	        unsigned int const s)
	{
	        unsigned int i = 0;
	        
	        while ( i < s )
	        {
	                if ( getNextPatternUnlocked(pat[i]) )
	                        ++i;
                        else
                                break;
	        }
	        
	        return i;
	}
	
	unsigned int fillIdBlock(FastIDBlock & block, unsigned int const s)
	{
	        block.setup(s);
	        block.idbase = nextid;
	        
	        unsigned int i = 0;
	        
	        while ( i < block.numid )
	        {
	                pattern_type pat;
	                
	                if ( getNextPatternUnlocked(pat) )
	                        block.ids[i++] = pat.sid;
                        else
                                break;
	        }
	        
	        block.blocksize = i;
	        
	        return i;
	}
	
	static int getOffset(std::string const /*inputfile*/)
	{
	        return -1;
	}	
};
#endif
