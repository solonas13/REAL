#if ! defined(CHARBUFFER_HPP)
#define CHARBUFFER_HPP

#include "AutoArray.hpp"

struct CharBuffer
{
	unsigned int buffersize;
	unsigned int length;
	AutoArray<char> abuffer;
	char * buffer;

	void bufferPush(int c)
	{
		if ( length == buffersize )
			expandBuffer();
		buffer[length++] = c;
	}
	void expandBuffer()
	{
		unsigned int newbuffersize = 2*buffersize;
		AutoArray<char> newabuffer(newbuffersize);
	
		std::copy ( abuffer.get(), abuffer.get()+buffersize, newabuffer.get() );
		
		buffersize = newbuffersize;
		abuffer = newabuffer;
		buffer = abuffer.get();
		
		std::cerr << "expanded buffer size to " << buffersize << std::endl;
	}

	CharBuffer(unsigned int const initialsize = 128)
	: buffersize(initialsize), length(0), abuffer(buffersize), buffer(abuffer.get())
	{
	
	}
	
	void assign(std::string & s)
	{
	        s.assign ( buffer, buffer + length );
	}
	void reset()
	{
	        length = 0;
	}
};
#endif
