#include "real_config.hpp"

#include "FastDecodeInfo.hpp"

std::ostream & operator<< (std::ostream & out, FastDecodeInfo const & fdi)
{
	out <<
	        "idbase=" << fdi.idbase <<
	        " idnbase=" << fdi.idnbase <<
	        " patlen=" << fdi.patlen <<
		" patpos=" << fdi.patpos <<
		" patid=" << fdi.patid <<
		" patnum=" << fdi.patnum <<
		" patnpos=" << fdi.patnpos <<
		" patnid=" << fdi.patnid  <<
		" patnnum=" << fdi.patnnum << std::endl;
	return out; 
}
