#if ! defined(FASTDECODER_HPP)
#define FASTDECODER_HPP

#include <map>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include "FastDecodeInfo.hpp"
#include "FastSubDecoder.hpp"
#include "FastIDDecoder.hpp"
#include "FastFileDecoderBase.hpp"

struct FastDecoder : public FastFileDecoderBase
{
	std::string const filename;
	std::map < unsigned int, FastDecodeInfo > info;
	u_int64_t numpat;
	unsigned int blocksize;
	
	static unsigned int const numbuffers = 16;

	u_int32_t readNumber(std::istream & in)
	{
		u_int8_t buf[4];
		in.read(reinterpret_cast<char *>(&buf[0]),4);

		
		#if 0
		unsigned int red = in.gcount();
		std::cerr << "Getcount " << red << std::endl;
		
		std::cerr << "buf[0]=" << static_cast<u_int32_t>(buf[0]) << std::endl;
		std::cerr << "buf[1]=" << static_cast<u_int32_t>(buf[1]) << std::endl;
		std::cerr << "buf[2]=" << static_cast<u_int32_t>(buf[2]) << std::endl;
		std::cerr << "buf[3]=" << static_cast<u_int32_t>(buf[3]) << std::endl;
		#endif

		// assert ( red == 4 );

		return
			(static_cast<u_int32_t>(buf[0])<<24)
			|
			(static_cast<u_int32_t>(buf[1])<<16)
			|
			(static_cast<u_int32_t>(buf[2])<< 8)
			|
			(static_cast<u_int32_t>(buf[3])<< 0)
			;
	}
	
	std::vector < unsigned int > getPatternLengthVector() const
	{
		std::vector < unsigned int > V;
		for ( std::map < unsigned int, FastDecodeInfo >::const_iterator ita = info.begin(); ita != info.end(); ++ita )
			V.push_back ( ita->first );
		return V;
	}
	
	static u_int64_t countPatterns(std::string const & filename)
	{
		FastDecoder FD(filename,default_blocksize);
		return FD.numpat;
	}
		
	FastDecoder(
		std::string const & rfilename,
		unsigned int const rblocksize
	)
	: filename(rfilename), numpat(0), blocksize(rblocksize)
	{
		std::ifstream istr(filename.c_str(),std::ios::binary);
		
		while ( istr )
		{
			u_int32_t const patl = readNumber ( istr );
			
			if ( istr )
			{
				FastDecodeInfo fdi;
				
				fdi.idbase = numpat;
				fdi.patlen = patl;
				std::cerr << "Pattern length " << patl << std::endl;
				unsigned int const acgtcodelength = (patl + 3)/4;
				unsigned int const acgtncodelength = (patl + 1)/2;
				
				// std::cerr << "Handling ACGT at offset " << istr.tellg() << std::endl;
				u_int32_t const acgtbl = readNumber(istr);
				u_int32_t const acgtdata = acgtbl-4;
				u_int32_t const acgtmagic = readNumber(istr);
				fdi.patpos = istr.tellg();
				fdi.patnum = acgtdata/acgtcodelength;
				// std::cerr << "acgt block length " << acgtbl << " number of patterns " << fdi.patnum << std::endl;
				assert ( acgtmagic == 0 );
				istr.seekg ( acgtbl - 4, std::ios_base::cur );

				// std::cerr << "Handling ACGTid at offset " << istr.tellg() << std::endl;
				u_int32_t const acgtidbl = readNumber(istr);
				u_int32_t const acgtidmagic = readNumber(istr);
				assert ( acgtidmagic == 1 );
				fdi.patid = istr.tellg();
				istr.seekg ( acgtidbl - 4, std::ios_base::cur );

				// std::cerr << "Handling ACGTN at offset " << istr.tellg() << std::endl;
				fdi.idnbase = fdi.idbase + fdi.patnum;
				u_int32_t const acgtnbl = readNumber(istr);
				u_int32_t const acgtnmagic = readNumber(istr);
				u_int32_t const acgtndata = acgtnbl - 4;
				fdi.patnnum = acgtndata/acgtncodelength;
				// std::cerr << "acgtn block length " << acgtnbl << " number of patterns " << fdi.patnnum << std::endl;
				assert ( acgtnmagic == 2 );
				fdi.patnpos = istr.tellg();
				istr.seekg ( acgtnbl - 4, std::ios_base::cur );

				// std::cerr << "Handling ACGTNid at offset " << istr.tellg() << std::endl;
				u_int32_t const acgtnidbl = readNumber(istr);
				u_int32_t const acgtnidmagic = readNumber(istr);
				assert ( acgtnidmagic == 3 );
				fdi.patnid = istr.tellg();
				istr.seekg ( acgtnidbl-4, std::ios_base::cur );
				
				std::cerr << fdi << std::endl;
				
				info[patl] = fdi;
				
				numpat += fdi.getNumPatterns();
			}
		}		
	}
	
	std::auto_ptr < FastSubDecoder > getDecoder(unsigned int const patl) const
	{
		if ( info.find(patl) != info.end() )
		{
			FastDecodeInfo const & fdi = info.find(patl)->second;			
			return std::auto_ptr < FastSubDecoder > ( new FastSubDecoder (filename, false, fdi.patlen, fdi.patnum, blocksize, fdi.patpos, numbuffers, fdi.idbase) );
		}
		else
		{
			return std::auto_ptr < FastSubDecoder >();
		}
	} 
	std::auto_ptr < FastSubDecoder > getNDecoder(unsigned int const patl) const
	{
		if ( info.find(patl) != info.end() )
		{
			FastDecodeInfo const & fdi = info.find(patl)->second;			
			return std::auto_ptr < FastSubDecoder > ( new FastSubDecoder (filename, true, fdi.patlen, fdi.patnnum, blocksize, fdi.patnpos, numbuffers, fdi.idnbase) );
		}
		else
		{
			return std::auto_ptr < FastSubDecoder >();
		}
	}
	
	std::auto_ptr < FastIDDecoder > getIdDecoder(unsigned int const patl) const
	{
		if ( info.find(patl) != info.end() )
		{
			FastDecodeInfo const & fdi = info.find(patl)->second;			
			return std::auto_ptr < FastIDDecoder > ( new FastIDDecoder (filename, fdi.patnum, fdi.patid, blocksize, fdi.idbase) );
		}
		else
		{
			std::cerr << "getIdDecoder called for non-existing pattern size" << std::endl;
			return std::auto_ptr < FastIDDecoder >();
		}

	}
	std::auto_ptr < FastIDDecoder > getIdNDecoder(unsigned int const patl) const
	{
		if ( info.find(patl) != info.end() )
		{
			FastDecodeInfo const & fdi = info.find(patl)->second;			
			return std::auto_ptr < FastIDDecoder > ( new FastIDDecoder (filename, fdi.patnnum, fdi.patnid, blocksize, fdi.idnbase) );
		}
		else
		{
			std::cerr << "getIdDecoder called for non-existing pattern size" << std::endl;
			return std::auto_ptr < FastIDDecoder >();
		}

	}
};
#endif
