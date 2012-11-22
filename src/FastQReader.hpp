#if ! defined(FASTQREADER_HPP)
#define FASTQREADER_HPP

#include "FastReaderBase.hpp"
#include "FASTQEntry.hpp"
#include "SpaceTable.hpp"
#include "CharTermTable.hpp"
#include "CharBuffer.hpp"
#include "PatternBlock.hpp"
#include "FastIDBlock.hpp"

struct FastQReader : public FastReaderBase, public SpaceTable
{
        typedef FastQReader reader_type;
        typedef FASTQEntry pattern_type;
        typedef PatternBlock<pattern_type> block_type;
        typedef AsynchronousStreamReaderData<reader_type> stream_data_type;
        typedef AsynchronousStreamReader< AsynchronousStreamReaderData<reader_type> > stream_reader_type;

        typedef reader_type idfile_type;
        typedef FastIDBlock idblock_type;
        typedef AsynchronousIdData<reader_type> stream_id_type;
        typedef AsynchronousStreamReader<stream_id_type> stream_idreader_type;

	CharTermTable atscanterm;
	CharTermTable plusscanterm;
	CharTermTable newlineterm;

        CharBuffer idbuffer;
        CharBuffer patbuffer;
        CharBuffer plusbuffer;
        CharBuffer qualbuffer;

	bool foundnextmarker;
	int qualityOffset;
	
	u_int64_t nextid;
	
	FastQReader(std::string const & filename, int rqualityOffset = 0, unsigned int numbuffers = 16, unsigned int bufsize = 16*1024)
	: FastReaderBase(filename,numbuffers,bufsize),
	  atscanterm('@'), plusscanterm('+'), newlineterm('\n'),
	  foundnextmarker(false),
	  qualityOffset(rqualityOffset), nextid(0)
	{
		findNextMarker();
	}
	
	void findNextMarker()
	{
		int c;
		
		while ( !atscanterm[(c=getNextCharacter())] )
		{
			// std::cerr << "Got " << static_cast<char>(c) << " without term, " << scanterm[c] << std::endl;
		}
		
		foundnextmarker = (c=='@');
	}
	
	u_int64_t countPatterns()
	{
	        pattern_type pattern;
	        u_int64_t count = 0;
	        
	        while ( skipPattern(pattern) )
	                ++count;

                return count;
	}
	
	static u_int64_t getPatternLength(std::string const & filename)
	{
	        FastQReader reader(filename);
	        pattern_type pat;
	        
	        if ( ! reader.getNextPatternUnlocked(pat) )
	                return 0;
                else
                        return pat.patlen;
	}

	static u_int64_t countPatterns(std::string const & filename)
	{
	        FastQReader reader(filename);
	        return reader.countPatterns();
	}

        bool skipPattern(pattern_type & pattern)
	{
	        
		if ( ! foundnextmarker )
			return false;
			
		foundnextmarker = false;
	
		int c;

		while ( !newlineterm[c=getNextCharacter()] )
		{}
		
		if ( c < 0 )
			return false;
			
		patbuffer.reset();
		while ( !plusscanterm[(c=getNextCharacter())] )
			if ( nospacetable[c] )
				patbuffer.bufferPush(c);

		pattern.patlen = patbuffer.length;
		
		while ( !newlineterm[c=getNextCharacter()] )
		{}
		
		if ( c < 0 )
			return false;
		
		qualbuffer.reset();
		while ( ((c=getNextCharacter()) >= 0) && (qualbuffer.length < pattern.patlen) )
		        if ( nospacetable[c] )
		                qualbuffer.bufferPush(c);

                if ( qualbuffer.length < pattern.patlen )
                        return false;
                        
		findNextMarker();

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
		while ( !plusscanterm[(c=getNextCharacter())] )
			if ( nospacetable[c] )
				patbuffer.bufferPush(c);

                patbuffer.assign(pattern.spattern);
		pattern.pattern = pattern.spattern.c_str();
		pattern.patlen = pattern.spattern.size();
		
		plusbuffer.reset();
		while ( !newlineterm[c=getNextCharacter()] )
			plusbuffer.bufferPush(c);
		
		if ( c < 0 )
			return false;
                plusbuffer.assign(pattern.plus);
		
		qualbuffer.reset();
		while ( ((c=getNextCharacter()) >= 0) && (qualbuffer.length < pattern.patlen) )
		        if ( nospacetable[c] )
		                qualbuffer.bufferPush(c-qualityOffset);

                if ( qualbuffer.length < pattern.patlen )
                        return false;
                        
                qualbuffer.assign ( pattern.quality );
                
                pattern.patid = nextid++;
		
		findNextMarker();

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

	static int getOffset(std::string const inputfile)
	{
	        FastQReader file ( inputfile );
	        FASTQEntry entry;
                while ( file.getNextPatternUnlocked(entry) )
                {
	                std::string const & quality = entry.quality;
	                
	                for ( unsigned int i = 0; i < quality.size(); ++i )
                                // Sanger
	                        if ( quality[i] <= 54 )
	                                return 33;
                                // Illumina
                                else if ( quality[i] >= 94 )
                                        return 64;                
                }
	        
	        return 0;
	}
};
#endif
