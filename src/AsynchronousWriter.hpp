#if ! defined(ASYNCHRONOUSWRITER_HPP)
#define ASYNCHRONOUSWRITER_HPP

#include <iostream>
#include <string>
#include <cstring>
#include <cerrno>
#include <stdexcept>
#include <cassert>
#include "AutoArray.hpp"
#include "Lock.hpp"

#if defined(HAVE_AIO_H)
#include <aio.h>

struct AsynchronousWriter
{
	int const fd;
	bool const releasefd;
	unsigned int const numbuffers;
	AutoArray < AutoArray<char> > buffers;
	AutoArray < aiocb > contexts;
	u_int64_t low;
	u_int64_t high;
	toollib::Lock lock;

	AsynchronousWriter ( int const rfd, u_int64_t rnumbuffers = 16 )
	: fd( rfd ), releasefd(false), numbuffers(rnumbuffers), buffers(numbuffers), contexts(numbuffers), low(0), high(0)
	{
		if ( fd < 0 )
		{
			throw std::runtime_error("Failed to open file.");
		}
	        fcntl (fd, F_SETFL, O_APPEND);
	}
	AsynchronousWriter ( std::string const & filename, u_int64_t rnumbuffers = 16 )
	: fd( open(filename.c_str(),O_WRONLY|O_CREAT|O_APPEND|O_TRUNC,0644) ), releasefd(true), numbuffers(rnumbuffers), buffers(numbuffers), contexts(numbuffers), low(0), high(0)
	{
		if ( fd < 0 )
		{
			throw std::runtime_error("Failed to open file.");
		}
	}
	~AsynchronousWriter()
	{
		flush();
		if ( releasefd )
		{
			close(fd);
		}
		else
		{
		}
	}
	
	template<typename iterator>
	void write(iterator sa, iterator se)
	{
		lock.lock();

		if ( high-low == numbuffers )
		{
			aiocb *waitlist[1] = { &contexts[low%numbuffers] };
			// std::cerr << "waiting for " << low << std::endl;
			aio_suspend (waitlist,1,0);
			low++;
		}
		
		u_int64_t const len = se-sa;
		// std::cerr << "writing " << s.size() << std::endl;

		buffers[high % numbuffers] = AutoArray<char>(len);
		std::copy ( sa, se, buffers[high%numbuffers].get() );
		memset ( &contexts[high%numbuffers], 0, sizeof(aiocb) );
		contexts[high%numbuffers].aio_fildes = fd;
		contexts[high%numbuffers].aio_buf = buffers[high % numbuffers].get();
		contexts[high%numbuffers].aio_nbytes = len;
		contexts[high%numbuffers].aio_offset = 0;
		contexts[high%numbuffers].aio_sigevent.sigev_notify = SIGEV_NONE;
		aio_write( & contexts[high%numbuffers] );		

		high++;		
		
		lock.unlock();
	}
	
	void flush()
	{
		while ( high-low )
		{
			aiocb *waitlist[1] = { &contexts[low%numbuffers] };
			aio_suspend (waitlist,1,0);
			ssize_t ret = aio_return ( &contexts[low%numbuffers] );
			assert ( ret == static_cast<ssize_t>(contexts[low%numbuffers].aio_nbytes) ) ;
			// std::cerr << "WRITTEN: " << ret << " requested " << contexts[low%numbuffers].aio_nbytes << std::endl;
			
			low++;
		}

		if ( fd != STDOUT_FILENO )
			if ( fsync(fd) )
				std::cerr << "Failure in fsync: " << strerror(errno) << std::endl;
	}
};
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
                     
struct AsynchronousWriter
{
	int fd;
	toollib::Lock lock;
	bool closefd;
	AutoArray<u_int8_t> buf;
	
	void expandBuf(u_int64_t size)
	{
		if ( buf.getN() < size )
			buf = AutoArray<u_int8_t>(size,false);
	}
	
	AsynchronousWriter ( std::string const & filename, u_int64_t )
	: fd(-1), closefd(true), buf()
	{
		fd = open(filename.c_str(),O_CREAT|O_TRUNC|O_WRONLY,0644);
		
		if ( fd < 0 )
		{
			std::cerr << "Failed to open file " << filename << std::endl;
			throw std::runtime_error("Failed to open file.");
		}
	}

	AsynchronousWriter ( int rfd, u_int64_t )
	: fd(rfd), closefd(false), buf()
	{
	
	}
	
	~AsynchronousWriter()
	{
		lock.lock();
		fsync(fd);
		if ( closefd )
			close(fd);
		lock.unlock();
	}

	template<typename iterator>
	void write(iterator rsa, iterator rse)
	{
		lock.lock();
		
		expandBuf(rse-rsa);
		std::copy ( rsa, rse, buf.get() );
		u_int8_t const * sa = buf.get();
		u_int8_t const * se = sa + (rse-rsa);
		
		while ( sa != se )
		{
			ssize_t const wr = ::write(fd,sa,se-sa);
			
			if ( wr <= 0 )
			{
				switch ( errno )
				{
					case EINTR:
						break;
					default:
						std::cerr << "write() failed with error code " << errno << ": " << strerror(errno) << " wr " << wr << std::endl;
						throw std::runtime_error("write() failed.");
						break;
				}
			}
			else
			{
				sa += wr;
			}
		}
		lock.unlock();
	}
	
	void flush()
	{
		lock.lock();
		fsync(fd);
		lock.unlock();
	}
};
#endif

#endif
