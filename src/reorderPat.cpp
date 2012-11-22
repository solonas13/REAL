#include "real_config.hpp"


#include "ReorderFastA.hpp"
#include "ReorderFastQ.hpp"

#include "RealTimeClock.hpp"
#include "RealOptions.hpp"

#include "FastAReader.hpp"
#include "FastQReader.hpp"

#include "FastFileDecoder.hpp"
#include "FastQualityFileDecoder.hpp"
#include "AsynchronousReader.hpp"

#include "FastIDFileDecoder.hpp"
#include "FastQualityIDFileDecoder.hpp"

int main(int argc, char * argv[])
{
	if ( argc < 3 )
		return EXIT_FAILURE;
	
	std::string const filename = argv[1];
	std::string const outfilename = argv[2];
	std::ofstream ostr(outfilename.c_str(), std::ios::binary);
	
	double clock_before = clock();
	toollib::RealTimeClock rtc;
	rtc.start();
	unsigned int qoffset;
	std::vector<std::string> remlist;
	if ( RealOptions::isFastQ(filename) )
	{
		typedef FastQReader reader_type;

		qoffset = FastQReader::getOffset(filename);
		reader_type reader(filename,qoffset);
		reorderFastQ(reader,ostr,remlist);
	}
	else
	{
		typedef FastAReader reader_type;
				
		reader_type reader(filename);
		reorderFastA(reader,ostr,remlist);
	}
	ostr.close();
	double clock_after = clock();
	
	if ( RealOptions::isFastQ(filename) )
	{
		typedef FastQualityFileDecoder reader_type;
		reader_type FFD(outfilename);		
		reader_type::stream_data_type AFRD(FFD,64*1024,16);
		reader_type::stream_reader_type AFR(AFRD);
		
		reader_type::block_type * pbb = 0;
		while ( (pbb = AFR.getBlock()) )
		{
			std::cerr << "blocksize: " << pbb->blocksize << std::endl;

			for ( unsigned int j = 0; j < pbb->blocksize; ++j )
			{
				PatternQualityBase const & pattern = pbb->patterns[j];

				std::string const clear = toollib::remapString (
					std::string(
						pattern.mapped,
						pattern.mapped + pbb->patlen
						)
					);
				std::cerr << pattern.patid << "\t" << clear << std::endl;
			}

			AFR.returnBlock(pbb);
		}
		
		std::cerr << "Number of patterns: " << reader_type::countPatterns(outfilename) << std::endl;

		reader_type::idfile_type IF(outfilename);
		reader_type::stream_id_type SI(IF,64*1024,16);
		reader_type::stream_idreader_type SR(SI);
		reader_type::idblock_type * ib = 0;
		
		while ( (ib = SR.getBlock()) )
		{
			for ( unsigned int i = 0; i < ib->blocksize; ++i )
				std::cerr << ib->idbase+i << "\t" << ib->ids[i] << std::endl;
		
			SR.returnBlock(ib);
		}
	}
	else
	{
		typedef FastFileDecoder reader_type;
		reader_type FFD(outfilename);		
		reader_type::stream_data_type AFRD(FFD,64*1024,16);
		reader_type::stream_reader_type AFR(AFRD);
		
		reader_type::block_type * pbb = 0;
		while ( (pbb = AFR.getBlock()) )
		{
			std::cerr << "blocksize: " << pbb->blocksize << std::endl;
			
			for ( unsigned int j = 0; j < pbb->blocksize; ++j )
			{
				PatternBase const & pattern = pbb->patterns[j];

				std::string const clear = toollib::remapString (
					std::string(
						pattern.mapped,
						pattern.mapped + pbb->patlen
						)
					);
				std::cerr << pattern.patid << "\t" << clear << std::endl;
			}
			
			AFR.returnBlock(pbb);
		}

		std::cerr << "Number of patterns: " << reader_type::countPatterns(outfilename) << std::endl;

		reader_type::idfile_type IF(outfilename);
		reader_type::stream_id_type SI(IF,64*1024,16);
		reader_type::stream_idreader_type SR(SI);
		reader_type::idblock_type * ib = 0;
		
		while ( (ib = SR.getBlock()) )
		{
			for ( unsigned int i = 0; i < ib->blocksize; ++i )
				std::cerr << ib->idbase+i << "\t" << ib->ids[i] << std::endl;
		
			SR.returnBlock(ib);
		}
	}
		
	std::cerr << "CPU time " << (clock_after-clock_before)/CLOCKS_PER_SEC << "s, real time " << rtc.getElapsedSeconds() << "s" << std::endl;
}
