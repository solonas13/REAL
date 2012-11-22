/**
    REAL: An efficient REad ALigner for next generation sequencing reads.
    Copyright (C) 2010 German Tischler, Solon Pissis

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#if defined(HAVE_CONFIG_H)
#include "real_config.hpp"
#endif

#if defined(GAPPED_MATCHING)	
#include "gappedMatching.hpp"
#endif

#include "matchAll.hpp"
#include "matchUnique.hpp"
#include "RealOptions.hpp"
#include "FastAReader.hpp"
#include "FastQReader.hpp"
#include "FastFileDecoder.hpp"
#include "FastQualityFileDecoder.hpp"
#include "types.hpp"
#include <csignal>
#include <signal.h>
#include <sys/types.h>

#include <iostream>

#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#include <set>
#include "stringFunctions.hpp"
std::string getFeatures()
{
	char featureBuffer[4096];
	size_t featureBufferLen = sizeof(featureBuffer);
	int const sysctlretname = sysctlbyname("machdep.cpu.features", &featureBuffer[0], &featureBufferLen, 0, 0);
	assert ( ! sysctlretname );
	std::string const features(&featureBuffer[0],&featureBuffer[0] + featureBufferLen);
	return features;
}
std::set<std::string> getFeatureVector()
{
	std::string const features = getFeatures();
	std::deque<std::string> featureVec = stringFunctions::tokenize(features,std::string(" "));
	std::set<std::string> featureSet;
	for ( size_t i = 0; i < featureVec.size(); ++i )
		featureSet.insert(featureVec[i]);	
	return featureSet;
}
bool hasSSE4PopCnt()
{
	std::set<std::string> featureSet = getFeatureVector();

	#if 0
	for ( std::set<std::string>::const_iterator ita = featureSet.begin(); ita != featureSet.end(); ++ita )
	{
		std::cerr << *ita << std::endl;
	}
	#endif

	return featureSet.find("POPCNT") != featureSet.end();
}
#else
void cpuid(
	u_int32_t & eax,
	u_int32_t & ebx,
	u_int32_t & ecx,
	u_int32_t & edx
)
{
	#if defined(HAVE_x86_64)
	typedef u_int64_t regsave_type;
	#else
	typedef u_int32_t regsave_type;
	#endif

	regsave_type * regsave = 0;
	
	try
	{
		regsave = new regsave_type[9];

		regsave[0] = eax;
		regsave[1] = ebx;
		regsave[2] = ecx;
		regsave[3] = edx;
		
		#if defined(HAVE_x86_64)
		asm volatile(
			"mov %%rax,(4*8)(%0)\n"
			"mov %%rbx,(5*8)(%0)\n"
			"mov %%rcx,(6*8)(%0)\n"
			"mov %%rdx,(7*8)(%0)\n"
			"mov %%rsi,(8*8)(%0)\n"
			"mov %0,%%rsi\n"

			"mov (0*8)(%%rsi),%%rax\n"
			"mov (1*8)(%%rsi),%%rbx\n"
			"mov (2*8)(%%rsi),%%rcx\n"
			"mov (3*8)(%%rsi),%%rdx\n"
			
			"cpuid\n"
			
			"mov %%rax,(0*8)(%%rsi)\n"
			"mov %%rbx,(1*8)(%%rsi)\n"
			"mov %%rcx,(2*8)(%%rsi)\n"
			"mov %%rdx,(3*8)(%%rsi)\n"

			"mov %%rsi,%0\n"
			
			"mov (4*8)(%0),%%rax\n"
			"mov (5*8)(%0),%%rbx\n"
			"mov (6*8)(%0),%%rcx\n"
			"mov (7*8)(%0),%%rdx\n"
			"mov (8*8)(%0),%%rsi\n"
			: : "a" (&regsave[0]) : "memory" ); 
		#else
		asm volatile(
			"mov %%eax,(4*4)(%0)\n"
			"mov %%ebx,(5*4)(%0)\n"
			"mov %%ecx,(6*4)(%0)\n"
			"mov %%edx,(7*4)(%0)\n"
			"mov %%esi,(8*4)(%0)\n"
			"mov %0,%%esi\n"

			"mov (0*4)(%%esi),%%eax\n"
			"mov (1*4)(%%esi),%%ebx\n"
			"mov (2*4)(%%esi),%%ecx\n"
			"mov (3*4)(%%esi),%%edx\n"
			
			"cpuid\n"
			
			"mov %%eax,(0*4)(%%esi)\n"
			"mov %%ebx,(1*4)(%%esi)\n"
			"mov %%ecx,(2*4)(%%esi)\n"
			"mov %%edx,(3*4)(%%esi)\n"

			"mov %%esi,%0\n"
			
			"mov (4*4)(%0),%%eax\n"
			"mov (5*4)(%0),%%ebx\n"
			"mov (6*4)(%0),%%ecx\n"
			"mov (7*4)(%0),%%edx\n"
			"mov (8*4)(%0),%%esi\n"
			: : "a" (&regsave[0]) : "memory" ); 
		#endif
		
		eax = regsave[0];
		ebx = regsave[1];
		ecx = regsave[2];
		edx = regsave[3];
		
		// std::cerr << "eax=" << eax << " ebx=" << ebx << " ecx=" << ecx << " edx=" << edx << std::endl;
		
		delete [] regsave;
		regsave = 0;
	}
	catch(...)
	{
		delete [] regsave;
		regsave = 0;
		throw;
	}
}

bool hasSSE4PopCnt()
{
	u_int32_t a=0,b=0,c=0,d=0;
	cpuid(a,b,c,d);

	if ( 1 <= a )
	{
		a = 1;
		cpuid(a,b,c,d);
		
		if ( c & (1ul << 23) )
			return true;
		else
			return false;
	}
	else
	{
	        std::cerr << "cpuid function 1 is not supported." << std::endl;
		return false;
	}
}
#endif

template<bool sse4, bool unique, typename reader_type, bool scores, typename word_type>
int cpuReaderScoresWordsizeMain(RealOptions const & opts)
{
        if ( unique )
                EnumerateUniqueMatches<sse4, word_type, reader_type, scores>::doMatching(opts);
        else
                EnumerateAllMatches<sse4, word_type, reader_type, scores>::doMatching(opts);
                
        return EXIT_SUCCESS;
}

template<bool sse4, bool unique, typename reader_type, bool scores>
int cpuReaderScoresMain(RealOptions const & opts)
{
        if ( opts.seedl <= 32 )
                return cpuReaderScoresWordsizeMain<sse4, unique, reader_type, scores, u_int32_t>(opts);
        else
                return cpuReaderScoresWordsizeMain<sse4, unique, reader_type, scores, u_int64_t>(opts);
}

template<bool sse4, bool unique, typename reader_type>
int cpuReaderMain(RealOptions const & opts)
{
        if ( opts.scores )
                return cpuReaderScoresMain<sse4,unique,reader_type,true>(opts);
        else
                return cpuReaderScoresMain<sse4,unique,reader_type,false>(opts);
}

#include "TemporaryFile.hpp"
#include "ReorderFastA.hpp"
#include "ReorderFastQ.hpp"
#include <fstream>

static std::vector < std::string > tempfilelist;

std::string reorderPatternFile(RealOptions const & opts)
{
        TemporaryName patcomp(tempfilelist);
        std::string patfilename = patcomp.filename;
        
        std::cerr << "Rewriting pattern file to temporary file " << patfilename << std::endl;
        
        if ( opts.fastq )
        {
        	if ( (! opts.qualityOffset) && (opts.patternfilename == "-") )
        	{
        		std::cerr << "WARNING: automatic quality offset detection not supported when" << std::endl;
        		std::cerr << "         reading patterns from standard input. Assuming input" << std::endl;
        		std::cerr << "         was produced by an Illumina  GA (i.e. -Q 64)" << std::endl;
        	}
        
                int const qualityOffset = opts.qualityOffset ? opts.qualityOffset 
                	: 
                	((opts.patternfilename == "-") ? 64 : FastQReader::getOffset(opts.patternfilename));

                if ( ! qualityOffset )
                        throw std::runtime_error("Unable to automatically detect FastQ quality format.");
                        
                std::ofstream ostr(patfilename.c_str(),std::ios::binary);
                FastQReader reader(opts.patternfilename,qualityOffset);
                reorderFastQ(reader,ostr,tempfilelist);
                ostr.flush();
                ostr.close();
        }
        else
        {
                std::ofstream ostr(patfilename.c_str(),std::ios::binary);
                FastAReader reader(opts.patternfilename);
                reorderFastA(reader,ostr,tempfilelist);
                ostr.flush();
                ostr.close();                
        }

        return patfilename;
}


void removeTemporaryFile()
{
	for ( unsigned int i = 0; i < tempfilelist.size(); ++i )
		remove ( tempfilelist[i].c_str() );
	tempfilelist.resize(0);
}

void sigIntHandler(int)
{
        removeTemporaryFile();
        _Exit(1);
}


template<bool sse4, bool unique>
int cpuMain(RealOptions & opts)
{
	if ( unique )
	{
	        if ( opts.rewritepatterns )
	        {
                        std::string const patcomp = reorderPatternFile(opts);
                        opts.patternfilename = patcomp;
                        int result;
                        
                        if ( opts.fastq )
                                result = cpuReaderMain<sse4,unique,FAST_UNIQUE_FASTQ_READER_TYPE>(opts);
                        else
                                result = cpuReaderMain<sse4,unique,FAST_UNIQUE_FASTA_READER_TYPE>(opts);

                        removeTemporaryFile();        
                        
                        return result;
                }
                else
                {
                        if ( opts.fastq )
                                return cpuReaderMain<sse4,unique,SLOW_UNIQUE_FASTQ_READER_TYPE>(opts);
                        else
                                return cpuReaderMain<sse4,unique,SLOW_UNIQUE_FASTA_READER_TYPE>(opts);                
                }
        }
        else
        {
        	if ( opts.fastq )
        	        return cpuReaderMain<sse4,unique,SLOW_ALL_FASTA_READER_TYPE>(opts);
                else
                        return cpuReaderMain<sse4,unique,SLOW_ALL_FASTA_READER_TYPE>(opts);        
        }
}

template<bool sse4>
int cpuMain(int argc, char * argv[])
{
	try
	{
		RealOptions opts(argc,argv);
		
		if ( opts.match_unique )
        		return cpuMain<sse4,true>(opts);
                else
                        return cpuMain<sse4,false>(opts);
        }
        catch(std::exception const & ex)
        {
        	std::cerr << ex.what() << std::endl;
        	return EXIT_FAILURE;
        }
        catch(...)
        {
        	std::cerr << "Caught unexpected exception, terminating." << std::endl;
        	return EXIT_FAILURE;
        }
}


int main (int argc, char **argv)
{
        std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << std::endl;
        std::cerr << "real is distributed under version 3.0 of the GNU GENERAL PUBLIC LICENSE." << std::endl;

        atexit(removeTemporaryFile);                        
        ::signal(SIGINT,sigIntHandler);
        ::signal(SIGTERM,sigIntHandler);

	bool const sse4popcnt = hasSSE4PopCnt();

	if ( sse4popcnt )
	        std::cerr << "Machine supports popcnt instruction." << std::endl;
	
	if ( sse4popcnt )
		cpuMain<true>(argc,argv);
	else
		cpuMain<false>(argc,argv);
}
