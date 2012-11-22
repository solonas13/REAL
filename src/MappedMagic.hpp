#if ! defined(MAPPEDMAGIC_HPP)
#define MAPPEDMAGIC_HPP

#include <string>
#include <cctype>
#include <fstream>

struct MappedMagic : std::ifstream
{
	MappedMagic ( std::string const & filename )
	: std::ifstream(filename.c_str())
	{
	
	}
	~MappedMagic()
	{
	        close();
	}
	
	int getFirstNonSpace()
	{
	        int c;
	        
	        while ( (c=get()) >= 0 )
	                if ( ! isspace(c) )
	                        return c;
                
                return c;
	}

	static int getFirstNonSpace(std::string const & filename)
	{
	        MappedMagic mm(filename);
	        return mm.getFirstNonSpace();
	}
};
#endif
