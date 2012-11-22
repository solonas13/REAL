#if ! defined(FASTREADERBASE_HPP)
#define FASTREADERBASE_HPP

#include "AsynchronousBufferReader.hpp"

struct FastReaderBase
{
	AsynchronousBufferReader reader;
	std::pair < char const *, ssize_t > data;
	u_int8_t const * text;
	u_int64_t textlength;
	u_int64_t dataptr;
	bool firstbuf;

	FastReaderBase(std::string const & filename, unsigned int numbuffers = 16, unsigned int bufsize = 1024, u_int64_t const offset = 0)
	: reader(filename,numbuffers,bufsize,offset), data(0,0), text(0), textlength(0), dataptr(0), firstbuf(true)
	{
	}

	int getNextCharacter()
	{
		if ( dataptr == textlength )
		{
			if ( !firstbuf )
				reader.returnBuffer();
			firstbuf = false;
			bool const readok = reader.getBuffer(data);
			if ( ! readok )
				return -1;
			text = reinterpret_cast<u_int8_t const *>(data.first);
			textlength = data.second;
			dataptr = 0;
		}
		
		return text[dataptr++];
	}
	
	bool getNumber(u_int32_t & num)
	{
		int a,b,c,d;
		num = 0;

		a = getNextCharacter();
		if ( a < 0 )
			return false;
		b = getNextCharacter();
		if ( b < 0 )
			return false;
		c = getNextCharacter();
		if ( c < 0 )
			return false;
		d = getNextCharacter();
		if ( d < 0 )
			return false;
			
		num = (static_cast<u_int32_t>(a)<<24) | (static_cast<u_int32_t>(b)<<16) | (static_cast<u_int32_t>(c)<< 8) | (static_cast<u_int32_t>(d)<< 0);
		
		return true;
	}
	bool getNumber2(u_int16_t & num)
	{
		int a,b;
		num = 0;

		a = getNextCharacter();
		if ( a < 0 )
			return false;
		b = getNextCharacter();
		if ( b < 0 )
			return false;
			
		num = (static_cast<u_int32_t>(a)<<8) | (static_cast<u_int32_t>(b)<<0);
		
		return true;
	}
	bool getString(std::string & s)
	{
		u_int16_t num;
		if ( ! getNumber2(num) )
		{
			return false;
		}
		
		s.resize(num);
		
		for ( unsigned int i = 0; i < num; ++i )
		{
			int c = getNextCharacter();
			
			if ( c < 0 )
				return false;
			
			s[i] = c;
		}

		return true;
	}
};
#endif
