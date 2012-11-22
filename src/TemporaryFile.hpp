#if ! defined(TEMPORARYFILE_HPP)
#define TEMPORARYFILE_HPP

#include "types.hpp"
#include "AsynchronousWriter.hpp"
#include <cstdio>
#include <ostream>
#include <memory>

struct TemporaryFileBase
{
	void writeNumber(u_int32_t const num, std::ostream & out)
	{
		out.put ( (num >> 24) & 0xFF );
		out.put ( (num >> 16) & 0xFF );
		out.put ( (num >>  8) & 0xFF );
		out.put ( (num >>  0) & 0xFF );
	}
};

struct TemporaryName
{
	std::string filename;
	
	TemporaryName( )
	{
		char buf[] = "tempfile.XXXXXX";
		int fd = mkstemp(&buf[0]);
		if ( fd < 0 )
			throw std::runtime_error("Unable to create temporary file.");
		close(fd);
		filename = &buf[0];
	}
	TemporaryName(std::vector < std::string > & remlist )
	{
		char buf[] = "tempfile.XXXXXX";
		int fd = mkstemp(&buf[0]);
		if ( fd < 0 )
			throw std::runtime_error("Unable to create temporary file.");
		close(fd);
		filename = &buf[0];
		remlist.push_back(filename);
	}
}; 

struct TemporaryFile : public TemporaryFileBase, public TemporaryName
{
	static unsigned int const bufsize = 16*1024;

	int fd;
	// std::string filename;
	u_int64_t numwritten;
	std::auto_ptr < AsynchronousWriter > AW;
	AutoArray<char> AB;
	char * const B;
	unsigned int buffill;
	unsigned int const magic;
	
	TemporaryFile(unsigned int rmagic, std::vector < std::string > & remlist )
	: TemporaryFileBase(), TemporaryName(remlist), fd(-1), numwritten(0), AB(bufsize), B(AB.get()), buffill(0), magic(rmagic)
	{
		fd = open(filename.c_str(), O_CREAT|O_APPEND|O_RDWR, 0600);
		
		if ( fd < 0 )
			throw std::runtime_error("Unable to reopen temporary file.");
		
		AW = std::auto_ptr < AsynchronousWriter > ( new AsynchronousWriter (fd,64) );
		// std::cerr << "Writing magic " << magic << std::endl;
		writeNumber(magic);
	}
	~TemporaryFile()
	{
		AW.reset(0);		
		
		if ( fd != -1 )
		{
			close(fd);
			fd = -1;
		}

		if ( filename.size() )
			remove(filename.c_str());	
	}
	
	void flush()
	{
		#if 0
		if ( buffill == 4 )
		{
			std::cerr << "***HERE FLUSHING Temp Buffer***" <<
				static_cast<unsigned int>(reinterpret_cast<u_int8_t const *>(B)[0]) << " " <<
				static_cast<unsigned int>(reinterpret_cast<u_int8_t const *>(B)[1]) << " " <<
				static_cast<unsigned int>(reinterpret_cast<u_int8_t const *>(B)[2]) << " " <<
				static_cast<unsigned int>(reinterpret_cast<u_int8_t const *>(B)[3]) << std::endl;
		}
		#endif
	
		AW->write(B,B+buffill);
		buffill = 0;
	}
	
	template<typename iterator>
	void append(iterator sa, iterator se)
	{
		while ( sa != se )
		{
			if ( buffill == bufsize )
				flush();

			B[buffill++] = *(sa++);
		}
	}

	void writeNumber(u_int32_t const num)
	{
		u_int8_t buf[4];

		buf[0] = (num >> 24) & 0xFF;
		buf[1] = (num >> 16) & 0xFF;
		buf[2] = (num >>  8) & 0xFF;
		buf[3] = (num >>  0) & 0xFF;
		
		append( &buf[0], &buf[0]+sizeof(buf) );
	}
	void writeNumber2(u_int16_t const num)
	{
		u_int8_t buf[2];

		buf[0] = (num >>  8) & 0xFF;
		buf[1] = (num >>  0) & 0xFF;
		
		append( &buf[0], &buf[0]+sizeof(buf) );
	}

	void write(std::string const & s)
	{
		++numwritten;
		append(s.begin(),s.end());
	}
	template<typename iterator>
	void write(iterator sa, iterator se)
	{
		++numwritten;
		append(sa,se);
	}
	
	bool writeContent(std::ostream & out)
	{
		flush();
		AW->flush();
		AW.reset();
		fsync(fd);

		u_int64_t towrite = lseek(fd,0,SEEK_END);

		close(fd);
		fd = -1;
		
		fd = open(filename.c_str(), O_RDONLY);
		
		if ( fd < 0 )
		{
			std::cerr << "Failed to reopen file " << filename << std::endl;
			return false;
		}
	
		unsigned int bufsize = (1ul << 16);
		AutoArray<char> Abuffer(bufsize);
		char * buffer = Abuffer.get();
		
		// std::cerr << "towrite = " << towrite << std::endl;

		TemporaryFileBase::writeNumber(towrite,out);
	
		while ( towrite )
		{
			ssize_t red = ::read(fd,buffer,bufsize);
			
			if ( red < 0 )
				return false;
			
			out.write ( buffer, red );
			
			if ( out.bad() )
				return false;
				
			towrite -= red;
		}
		
		return true;
	}
};

template<typename pattern_type>
struct TemporaryFileBunch : public TemporaryFileBase
{
	unsigned int const patl;
	TemporaryFile TV;
	TemporaryFile TVid;
	TemporaryFile TVN;
	TemporaryFile TVNid;
	
	TemporaryFileBunch(unsigned int const rpatl, std::vector< std::string > & remlist)
	: patl(rpatl), TV(0,remlist), TVid(1,remlist), TVN(2,remlist), TVNid(3,remlist)
	{
	}
	
	void writePatternWithDontCares(pattern_type & pattern)
	{
		std::string & mapped = pattern.smapped;
		
		for ( unsigned int i = 0; i < (patl>>1); ++i )
		{
			u_int8_t const a = mapped[2*i+0];
			u_int8_t const b = mapped[2*i+1];
			
			u_int8_t const c =  (a << 4) | (b << 0);
			
			mapped[i] = c;			
		}

		if ( patl & 1 )
			mapped[patl>>1] = mapped[patl-1] << 4;
		
		unsigned int const codelength = (patl+1)>>1;
					
		TVN.write(mapped.begin(), mapped.begin()+codelength);
		TVNid.writeNumber2 ( pattern.sid.size() );
		TVNid.write ( pattern.sid );	
	}
	void writePatternDontCareFree(pattern_type & pattern)
	{
		std::string & mapped = pattern.smapped;
		
		for ( unsigned int i = 0; i < (patl>>2); ++i )
		{
			u_int8_t const c = 
				(static_cast<u_int8_t>(mapped[4*i+0]) << 6)
				|
				(static_cast<u_int8_t>(mapped[4*i+1]) << 4)
				|
				(static_cast<u_int8_t>(mapped[4*i+2]) << 2)
				|
				(static_cast<u_int8_t>(mapped[4*i+3]) << 0)
				;
			mapped[i] = c;
		}

		unsigned int const codelength = (patl+3)>>2;

		if ( ((patl>>2)<<2) != patl )
		{
			u_int8_t c = 0;

			if ( ((patl>>2)<<2)+0 < patl )
				c |= static_cast<u_int8_t>(mapped[((patl>>2)<<2)+0]) << 6;
			if ( ((patl>>2)<<2)+1 < patl )
				c |= static_cast<u_int8_t>(mapped[((patl>>2)<<2)+1]) << 4;
			if ( ((patl>>2)<<2)+2 < patl )
				c |= static_cast<u_int8_t>(mapped[((patl>>2)<<2)+2]) << 2;
				
			mapped [ (patl>>2) ] = c;
		}
					
		TV.write(mapped.begin(), mapped.begin()+codelength);
		TVid.writeNumber2 ( pattern.sid.size() );
		TVid.write ( pattern.sid );	
	}
	
	void writePattern(pattern_type & pattern)
	{
		if ( pattern.isDontCareFree() )
			writePatternDontCareFree(pattern);
		else
			writePatternWithDontCares(pattern);
	}
	
	void writeContent(std::ostream & out)
	{
		writeNumber(patl,out);
		// std::cerr << "Writing ACGT at offset " << out.tellp() << std::endl;
		bool oka = TV.writeContent(out);
		// std::cerr << "Writing ACGTid at offset " << out.tellp() << std::endl;
		bool okb = TVid.writeContent(out);
		// std::cerr << "Writing ACGTN at offset " << out.tellp() << std::endl;
		bool okc = TVN.writeContent(out);
		// std::cerr << "Writing ACGTNid at offset " << out.tellp() << std::endl;
		bool okd = TVNid.writeContent(out);
		
		assert ( oka );
		assert ( okb );
		assert ( okc );
		assert ( okd );
	}
};

template<typename pattern_type>
struct TemporaryQualityFileBunch : public TemporaryFileBase
{
	unsigned int const patl;
	TemporaryFile TV;
	TemporaryFile TVid;
	TemporaryFile TVN;
	TemporaryFile TVNid;
	
	TemporaryQualityFileBunch(unsigned int const rpatl, std::vector<std::string> & remlist)
	: patl(rpatl), TV(0,remlist), TVid(1,remlist), TVN(2,remlist), TVNid(3,remlist)
	{
	}
	
	void writePatternWithDontCares(pattern_type & pattern)
	{
		std::string & mapped = pattern.smapped;
		
		for ( unsigned int i = 0; i < (patl>>1); ++i )
		{
			u_int8_t const a = mapped[2*i+0];
			u_int8_t const b = mapped[2*i+1];
			
			u_int8_t const c =  (a << 4) | (b << 0);
			
			mapped[i] = c;			
		}

		if ( patl & 1 )
			mapped[patl>>1] = mapped[patl-1] << 4;
		
		unsigned int const codelength = (patl+1)>>1;
					
		TVN.write(mapped.begin(), mapped.begin()+codelength);
		TVN.write(pattern.quality.begin(), pattern.quality.end());
		TVNid.writeNumber2 ( pattern.sid.size() );
		TVNid.write ( pattern.sid );	
	}
	void writePatternDontCareFree(pattern_type & pattern)
	{
		std::string & mapped = pattern.smapped;
		
		for ( unsigned int i = 0; i < (patl>>2); ++i )
		{
			u_int8_t const c = 
				(static_cast<u_int8_t>(mapped[4*i+0]) << 6)
				|
				(static_cast<u_int8_t>(mapped[4*i+1]) << 4)
				|
				(static_cast<u_int8_t>(mapped[4*i+2]) << 2)
				|
				(static_cast<u_int8_t>(mapped[4*i+3]) << 0)
				;
			mapped[i] = c;
		}

		unsigned int const codelength = (patl+3)>>2;

		if ( ((patl>>2)<<2) != patl )
		{
			u_int8_t c = 0;

			if ( ((patl>>2)<<2)+0 < patl )
				c |= static_cast<u_int8_t>(mapped[((patl>>2)<<2)+0]) << 6;
			if ( ((patl>>2)<<2)+1 < patl )
				c |= static_cast<u_int8_t>(mapped[((patl>>2)<<2)+1]) << 4;
			if ( ((patl>>2)<<2)+2 < patl )
				c |= static_cast<u_int8_t>(mapped[((patl>>2)<<2)+2]) << 2;
				
			mapped [ patl>>2 ] = c;
		}
					
		TV.write(mapped.begin(), mapped.begin()+codelength);
		TV.write(pattern.quality.begin(), pattern.quality.end());
		TVid.writeNumber2 ( pattern.sid.size() );
		TVid.write ( pattern.sid );	
	}
	
	void writePattern(pattern_type & pattern)
	{
		if ( pattern.isDontCareFree() )
			writePatternDontCareFree(pattern);
		else
			writePatternWithDontCares(pattern);
	}
	
	void writeContent(std::ostream & out)
	{
		writeNumber(patl,out);
		// std::cerr << "Writing ACGT at offset " << out.tellp() << std::endl;
		bool oka = TV.writeContent(out);
		// std::cerr << "Writing ACGTid at offset " << out.tellp() << std::endl;
		bool okb = TVid.writeContent(out);
		// std::cerr << "Writing ACGTN at offset " << out.tellp() << std::endl;
		bool okc = TVN.writeContent(out);
		// std::cerr << "Writing ACGTNid at offset " << out.tellp() << std::endl;
		bool okd = TVNid.writeContent(out);
		
		assert ( oka );
		assert ( okb );
		assert ( okc );
		assert ( okd );
	}
};
#endif
