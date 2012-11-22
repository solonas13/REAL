#if ! defined(REORDERFASTQ_HPP)
#define REORDERFASTQ_HPP

#include <ostream>
#include <vector>
#include "TemporaryFile.hpp"

template<typename pattern_type>
void handleFastQPattern(
	pattern_type & pattern, std::vector < TemporaryQualityFileBunch<pattern_type> * > & TV,
	std::vector<std::string> & remlist
	)
{
	pattern.computeMapped();
	unsigned int const patl = pattern.getPatternLength();
	
	while ( ! (patl < TV.size()) )
		TV.push_back(0);
	if ( ! TV[patl] )
		TV[patl] = new TemporaryQualityFileBunch<pattern_type>(patl,remlist);

	TV[patl]->writePattern(pattern);
}

template<typename reader_type>
void reorderFastQ(reader_type & FA, std::ostream & out, std::vector<std::string> & remlist)
{
	typedef typename reader_type::pattern_type pattern_type;
	std::vector < TemporaryQualityFileBunch<pattern_type> * > TV;
	pattern_type pattern;

	while ( FA.getNextPatternUnlocked(pattern) )
	{
		handleFastQPattern(pattern, TV, remlist);
		
		if ( ((pattern.getPatID()+1) & ((1ull << 20)-1)) == 0 )
		{
			std::cerr << "\r                                   \r" << (pattern.getPatID()+1) << std::flush;
		}
	}
	std::cerr << "\r                                   \r" << (pattern.getPatID()+1) << std::endl;
	
	std::cerr << "Concatenating temporary files...";
	for ( unsigned int i = 0; i < TV.size(); ++i )
		if ( TV[i] )
		{
			TV[i]->writeContent(out);
			delete TV[i];
			TV[i];
		}
	
	out.flush();
	std::cerr << "done." << std::endl;
}
#endif
