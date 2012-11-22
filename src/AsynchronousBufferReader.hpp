#if ! defined(ASYNCHRONOUSBUFFERREADER_HPP)
#define ASYNCHRONOUSBUFFERREADER_HPP

#include <iostream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <cassert>
#include "AutoArray.hpp"

#if defined(HAVE_AIO_H)
#include <aio.h>
#else
#include <fstream>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <cerrno>

struct AsynchronousBufferReader
{
	std::string filename;
	bool closefd;
	int const fd;
	unsigned int const numbuffers;
	unsigned int const bufsize;
	AutoArray < char > bufferspace;
	AutoArray < char * > buffers;
	#if defined(HAVE_AIO_H)
	AutoArray < aiocb > contexts;
	#else
	std::ifstream istr;
	std::istream * Pistr;
	#endif
	
	u_int64_t low;
	u_int64_t high;
	u_int64_t offset;
	
	#if defined(HAVE_AIO_H)
	void enqueRead()
	{
		memset ( &contexts[high%numbuffers], 0, sizeof(aiocb) );
		contexts[high%numbuffers].aio_fildes = fd;
		contexts[high%numbuffers].aio_buf = buffers[high % numbuffers];
		contexts[high%numbuffers].aio_nbytes = bufsize;
		contexts[high%numbuffers].aio_offset = offset + high*bufsize;
		contexts[high%numbuffers].aio_sigevent.sigev_notify = SIGEV_NONE;
		aio_read( & contexts[high%numbuffers] );
		high++;
	}
	void waitForNextBuffer()
	{
		aiocb *waitlist[1] = { &contexts[low%numbuffers] };
		int r = aio_suspend (waitlist,1,0);	
		
		if ( r != 0 )
		{
			std::cerr << "ERROR: " << strerror(errno) << std::endl;
		}
	}
	#endif

	void flush()
	{
		#if defined(HAVE_AIO_H)
		if ( filename != "-" )
		{
			while ( high-low )
			{
				waitForNextBuffer();
				low++;
			}
		}
		#endif
	}
	
	~AsynchronousBufferReader()
	{
		flush();
		if ( closefd )
			close(fd);
	}

	AsynchronousBufferReader ( std::string const & rfilename, u_int64_t rnumbuffers = 16, u_int64_t rbufsize = 32, u_int64_t roffset = 0 )
	: filename(rfilename),
	  #if defined(HAVE_AIO_H)
	  closefd(filename != "-"),
	  fd( (filename != "-") ? open(filename.c_str(),O_RDONLY ) : STDIN_FILENO ), 
	  #else
	  closefd(false),
	  fd(-1),
	  #endif
	  numbuffers(rnumbuffers), bufsize(rbufsize), 
	  bufferspace ( numbuffers * bufsize ),
	  buffers ( numbuffers ),
	  #if defined(HAVE_AIO_H)
	  contexts(numbuffers), 
	  #endif
	  low(0), high(0), offset(roffset)
	{
		#if defined(HAVE_AIO_H)
		if ( fd < 0 )
		{
			std::cerr << "Failed to open file " << filename << std::endl;
			throw std::runtime_error("Failed to open file.");
		}
		#endif
		
		for ( unsigned int i = 0; i < numbuffers; ++i )
			buffers[i] = bufferspace.get() + i*bufsize;

		#if defined(HAVE_AIO_H)
		if ( filename != "-" )
			while ( high < numbuffers )
				enqueRead();
		#else
		if ( filename == "-" )
		{
			Pistr = &(std::cin);
		}
		else
		{
			istr.open(filename.c_str());
			if ( ! istr.is_open() )
			{
				std::cerr << "Failed to open file " << filename << std::endl;
				throw std::runtime_error("Failed to open file.");	
			}
			Pistr = &istr;
			
			if ( offset )
			{
				istr.seekg(offset,std::ios::beg);
			}
		}
		#endif
	}
	
	
	bool getBuffer(std::pair < char const *, ssize_t > & data)
	{
		ssize_t red = -1;

		#if defined(HAVE_AIO_H)
		if ( filename != "-" )
		{
			assert ( low != high );
			waitForNextBuffer();
			red = aio_return(&contexts[low%numbuffers]);
		}
		else
		{
			std::cin.read ( buffers[low%numbuffers], bufsize );
			red = std::cin.gcount();
		}
		#else
		Pistr->read ( buffers[low%numbuffers], bufsize );
		red = Pistr->gcount();
		#endif
		
		data = std::pair < char const *, ssize_t > (buffers[low%numbuffers],red );
		
		low++;		
		return red > 0;
	}
	void returnBuffer()
	{
		#if defined(HAVE_AIO_H)
		if ( filename != "-" )
			enqueRead();	
		#endif
	}
};

#endif
