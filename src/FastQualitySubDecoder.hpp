#if ! defined(FASTQUALITYSUBDECODER_HPP)
#define FASTQUALITYSUBDECODER_HPP

#include "AsynchronousBufferReader.hpp"
#include "PatternQualityBaseBlock.hpp"

struct FastQualitySubDecoder
{
	bool const dontcares;
	unsigned int const patlen;
	unsigned int const numpat;
	unsigned int const patperblock;
	unsigned int const acgtcodelength;
	unsigned int const acgtncodelength;
	unsigned int const usedcodelength;
	unsigned int decoded;
	AsynchronousBufferReader buffer;
	unsigned int patid;
	
	FastQualitySubDecoder ( 
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
	  usedcodelength ( dontcares? acgtncodelength:acgtcodelength ), decoded(0),
	  buffer(filename, rnumbuffers, patperblock*(usedcodelength+patlen), roffset),
	  patid(rpatid) {}
	  
	 AutoArray<char> allocateBuffer(unsigned int const blocksize)
	 {
	 	return AutoArray<char>(blocksize * patlen);
	 }
	 AutoArray<char> allocateQualityBuffer(unsigned int const blocksize)
	 {
	 	return AutoArray<char>(blocksize * patlen);
	 }
	 
	 void setupBlock(PatternQualityBaseBlock & pbb, unsigned int const blocksize)
	 {
	        // std::cerr << "blocksize " << blocksize << " patperblock " << patperblock << std::endl;
	        assert ( blocksize == patperblock );
	 	pbb.setup(allocateBuffer(blocksize), allocateBuffer(blocksize), allocateQualityBuffer(blocksize), patperblock, patlen);
	 }
	 	  
	 bool fillPatternBlock(PatternQualityBaseBlock & pbb, unsigned int const blocksize)
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
		// std::cerr << "Got buffer of size " << data.second << " used code length " << usedcodelength << std::endl;
		assert ( std::min(static_cast<u_int64_t>(data.second),static_cast<u_int64_t>((numpat-decoded)*(usedcodelength+patlen))) % (usedcodelength+patlen) == 0 );

		pbb.blocksize = data.second / (usedcodelength+patlen);
		pbb.blocksize = std::min ( 
			static_cast < u_int64_t >(pbb.blocksize), static_cast < u_int64_t >(numpat-decoded) );
		
		if ( dontcares )
		{
			u_int8_t const * s = reinterpret_cast<u_int8_t const *>(data.first);
			char * t = pbb.mappedBlock.get();
			char * tt = pbb.transposedBlock.get() + patlen;
			char * q = pbb.quality.get();
			
			for ( unsigned int i = 0; i < pbb.blocksize; ++i )
			{
			        pbb.patterns[i].patlen = patlen;
				pbb.patterns[i].mapped = t;
				pbb.patterns[i].patid = patid++;
				pbb.patterns[i].quality = q;
				
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

				for ( unsigned int j = 0; j < patlen; ++j )
					*(q++) = *(s++);
				
				pbb.patterns[i].transposed = tt;
				tt += 2*patlen;
			}		
		}
		else
		{
			u_int8_t const * s = reinterpret_cast<u_int8_t const *>(data.first);
			char * t = pbb.mappedBlock.get();
			char * tt = pbb.transposedBlock.get()+patlen;
			char * q = pbb.quality.get();
			
			for ( unsigned int i = 0; i < pbb.blocksize; ++i )
			{
			        pbb.patterns[i].patlen = patlen;
				pbb.patterns[i].mapped = t;
				pbb.patterns[i].patid = patid++;
				pbb.patterns[i].quality = q;
				
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

				for ( unsigned int j = 0; j < patlen; ++j )
					*(q++) = *(s++);
				
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
