#if ! defined(FASTSUBDECODER_HPP)
#define FASTSUBDECODER_HPP

#include "AsynchronousBufferReader.hpp"
#include "PatternBaseBlock.hpp"

struct FastSubDecoder
{
	bool const dontcares;
	// length of single pattern
	unsigned int const patlen;
	// all patterns to be decoded
	unsigned int const numpat;
	// patterns per block
	unsigned int const patperblock;
	unsigned int const acgtcodelength;
	unsigned int const acgtncodelength;
	unsigned int const usedcodelength;
	unsigned int decoded;
	AsynchronousBufferReader buffer;
	unsigned int patid;
	
	FastSubDecoder ( 
		std::string const & filename,
		bool const rdontcares,
		unsigned int const rpatlen,
		unsigned int const rnumpat,
		unsigned int const rpatperblock,
		u_int64_t roffset,
		u_int64_t rnumbuffers = 16,
		unsigned int rpatid = 0
		)
	: dontcares(rdontcares), 
	  patlen(rpatlen), 
	  numpat(rnumpat), 
	  patperblock(rpatperblock),
	  acgtcodelength((patlen + 3)/4), acgtncodelength((patlen + 1)/2),
	  usedcodelength ( dontcares?acgtncodelength:acgtcodelength ), decoded(0),
	  buffer(filename, rnumbuffers, patperblock*usedcodelength, roffset),
	  patid(rpatid) {}
	  
	 AutoArray<char> allocateBuffer(unsigned int const blocksize)
	 {
	 	return AutoArray<char>(blocksize * patlen,false);
	 }
	 
	 void setupBlock(PatternBaseBlock & pbb, unsigned int const blocksize)
	 {
	         assert ( blocksize == patperblock );
	 	pbb.setup(allocateBuffer(blocksize), allocateBuffer(blocksize), blocksize, patlen);
	 }
	 	  
	 bool fillPatternBlock(PatternBaseBlock & pbb, unsigned int const blocksize)
	 {
	         assert ( blocksize == patperblock );
	 
	         if ( decoded == numpat )
	 	        return false;

	 	std::pair < char const *, ssize_t > data;
	 	if ( !buffer.getBuffer(data) )
	 	{
			// return false;	 
			assert ( false );
		}
		assert ( std::min(static_cast<u_int64_t>(data.second),static_cast<u_int64_t>((numpat-decoded)*usedcodelength)) % usedcodelength == 0 );

		pbb.blocksize = data.second / usedcodelength;
		pbb.blocksize = std::min ( 
			static_cast < u_int64_t >(pbb.blocksize), static_cast < u_int64_t >(numpat-decoded) );
		
		if ( dontcares )
		{
			u_int8_t const * s = reinterpret_cast<u_int8_t const *>(data.first);
			char * t = pbb.mappedBlock.get();
			char * tt = pbb.transposedBlock.get() + patlen;
			
			for ( unsigned int i = 0; i < pbb.blocksize; ++i )
			{
			        pbb.patterns[i].patlen = patlen;
				pbb.patterns[i].mapped = t;
				pbb.patterns[i].patid = patid++;
				
				for ( unsigned int j = 0; j < patlen/2; ++j )
				{
					*(t++) = ((*s)>>4)&0xf;
					*(--tt) = toollib::invertN(((*s)>>4)&0xf);
					*(t++) = ((*s)>>0)&0xf;
					*(--tt) = toollib::invertN(((*s)>>0)&0xf);
					s++;
				} 	
				switch ( patlen % 2 )
				{
					case 0:
						break;
					case 1:
						*(t++) = ((*s)>>4)&0xf;
						*(--tt) = toollib::invertN(((*s)>>4)&0xf);
						s++;
						break;
				} 
				
				pbb.patterns[i].transposed = tt;
				tt += 2*patlen;
			}		
		}
		else
		{
			u_int8_t const * s = reinterpret_cast<u_int8_t const *>(data.first);
			char * t = pbb.mappedBlock.get();
			char * tt = pbb.transposedBlock.get()+patlen;
			
			for ( unsigned int i = 0; i < pbb.blocksize; ++i )
			{
			        pbb.patterns[i].patlen = patlen;
				pbb.patterns[i].mapped = t;
				pbb.patterns[i].patid = patid++;
				
				for ( unsigned int j = 0; j < patlen/4; ++j )
				{
					*(t++) = ((*s)>>6)&0x3;
					*(t++) = ((*s)>>4)&0x3;
					*(t++) = ((*s)>>2)&0x3;
					*(t++) = ((*s)>>0)&0x3;
					*(--tt) = toollib::invertN(((*s)>>6)&0x3);
					*(--tt) = toollib::invertN(((*s)>>4)&0x3);
					*(--tt) = toollib::invertN(((*s)>>2)&0x3);
					*(--tt) = toollib::invertN(((*s)>>0)&0x3);
					s++;
				}
				switch ( patlen % 4 )
				{
					case 0:
						break;
					case 1:
						*(t++) = ((*s)>>6)&0x3;
						*(--tt) = toollib::invertN(((*s)>>6)&0x3);
						s++;
						break;
					case 2:
						*(t++) = ((*s)>>6)&0x3;
						*(t++) = ((*s)>>4)&0x3;
						*(--tt) = toollib::invertN(((*s)>>6)&0x3);
						*(--tt) = toollib::invertN(((*s)>>4)&0x3);
						s++;
						break;
					case 3:
						*(t++) = ((*s)>>6)&0x3;
						*(t++) = ((*s)>>4)&0x3;
						*(t++) = ((*s)>>2)&0x3;
						*(--tt) = toollib::invertN(((*s)>>6)&0x3);
						*(--tt) = toollib::invertN(((*s)>>4)&0x3);
						*(--tt) = toollib::invertN(((*s)>>2)&0x3);
						s++;
						break;
				}
				
				pbb.patterns[i].transposed = tt;
				tt += 2*patlen;
			}		
		}
		
		buffer.returnBuffer();
		
		decoded += pbb.blocksize;
		
		return true;
	 }
};
#endif
