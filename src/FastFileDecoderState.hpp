#if ! defined(FASTFILEDECODERSTATE_HPP)
#define FASTFILEDECODERSTATE_HPP

#include <vector>

struct FastFileDecoderState
{
	unsigned int idpos;
	unsigned int subidpos;
	std::vector < unsigned int > patlen;
	
	FastFileDecoderState(std::vector < unsigned int > rpatlen) : idpos(0), subidpos(0), patlen(rpatlen)
	{
		
	}
	
	bool getNext(unsigned int & len, bool & dc)
	{
		if ( idpos >= patlen.size() )
			return false;
		
		len = patlen[idpos];
		dc = (subidpos != 0);
		
		subidpos++;
		
		if ( subidpos > 1 )
		{
			subidpos = 0;
			idpos += 1;
		}
		
		return true;
	}
};
#endif
