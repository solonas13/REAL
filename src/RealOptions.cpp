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

#include "real_config.hpp"

#include "RealOptions.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <memory>

#include "FastAReader.hpp"
#include "FastQReader.hpp"
#include "MappedMagic.hpp"

#include "getPhysicalMemory.hpp"
#include "stringFunctions.hpp"
#include "Scoring.hpp"
#include "UniqueMatchInfo.hpp"

#if defined(_OPENMP)
#include <omp.h>
#endif

#define ALLOW_MATCH_ALL

bool RealOptions::isFastQ(std::istream & istr)
{
	int first = istr.get();
	
	if ( first >= 0 )
	{
		istr.unget();

	        if ( first == '>' )
        	        return false;
	        else if ( first == '@' )
        	        return true;
	        else
        	        throw std::runtime_error("Unable to determine type of pattern file.");
	}
	else
	{
		throw std::runtime_error("Failed to read first character from pattern file.");
	}	
}

bool RealOptions::isFastQ(std::string const & filename)
{
	std::ifstream istr(filename.c_str());
	
	if ( ! istr.is_open() )
                throw std::runtime_error("Unable to open pattern file.");

	return isFastQ(istr);
}

#if defined(_OPENMP)
static unsigned int getDefaultNumThreads()
{
        if ( getenv("OMP_NUM_THREADS") )
                return std::max(1,atoi(getenv("OMP_NUM_THREADS")));
        else
                return omp_get_num_procs();
}
#endif

void RealOptions::printHelp()
{
        std::cerr << "Options:" << std::endl;        
        std::cerr << "-t <textfilename>" << std::endl;
        std::cerr << "-p <patternfilename>" << std::endl;
        std::cerr << "-o <outputfilename>" << std::endl;
        std::cerr << "-s <maximum number of errors in seed, default="<< default_seedkmax <<">" << std::endl;
        std::cerr << "-e <total maximum number of errors, default="<<default_totalkmax<<">" << std::endl;
        std::cerr << "-l <length of seed, default="<<default_seedl<<">" << std::endl;
#if defined(ALLOW_MATCH_ALL)
        std::cerr << "-u <search for unique match, default="<<default_match_unique<< ">" << std::endl;
#endif
        std::cerr << "-f <fraction of physical memory to use, default="<< default_fracmem <<">" << std::endl;
        std::cerr << "-q <use quality scores, default="<<default_scores<<">" << std::endl;
        std::cerr << "-Q <offset for quality scores, default=autodetect>" << std::endl;
        std::cerr << "-R <rewrite pattern file, default="<<default_rewritepatterns <<">" << std::endl;
#if defined(_OPENMP)
        std::cerr << "-T <number of matching threads, default=" << getDefaultNumThreads() << ">" << std::endl;
#endif

        std::cerr << "-similarity <sequence similarity, default=" << Scoring::DFLT_SIMILARITY << ">" << std::endl;
        /**
         * not used
         **/
        /*
        std::cerr << "-err <error rate, default=" << Scoring::DFLT_ERR << ">" << std::endl;
         */
        std::cerr << "-trans <transitions fraction of mutations, default=" << Scoring::DFLT_TRANS << ">" << std::endl;
        std::cerr << "-gc <composition bias, default=" << Scoring::DFLT_GC << ">" << std::endl;
        std::cerr << "-gcmut_bias <mutability bias of G&C, default=" << Scoring::DFLT_GCMUT_BIAS << ">" << std::endl;
        std::cerr << "-filter_level <filtering level for equal hits 0-4, default=" << default_filter_level << ">" << std::endl;
        /** implementation not completed **/
        /*
        std::cerr << "-g <search for gapped matches, default="<<default_gaps<<">" << std::endl;
         */
}


RealOptions::RealOptions(int argc, char * argv[])
: seedkmax(default_seedkmax), totalkmax(default_totalkmax), seedl(default_seedl), match_unique(default_match_unique), fracmem(default_fracmem), scores(default_scores), qualityOffset(0), usemem(0),
  rewritepatterns(default_rewritepatterns), sort_threads(default_sort_threads),
  filter_level(default_filter_level),
  similarity(Scoring::DFLT_SIMILARITY),
  err(Scoring::DFLT_ERR),
  trans(Scoring::DFLT_TRANS),
  gc(Scoring::DFLT_GC),
  gcmut_bias(Scoring::DFLT_GCMUT_BIAS),
  gaps(default_gaps)
{
        std::vector < std::string > opts;
        for ( int i = 1 ; i < argc; ++i )
                opts.push_back ( argv[i] );
                
        unsigned int i = 0;
        
        
        while ( i < opts.size() )
        {
                if ( opts[i] == "-t" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -t is missing.");
                        textfilename = opts[i+1];
                        i += 2;
                }
                else if ( opts[i] == "-p" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -p is missing.");
                        patternfilename = opts[i+1];
                        i += 2;
                }
                else if ( opts[i] == "-o" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -o is missing.");
                        outputfilename = opts[i+1];
                        i += 2;
                }
                else if ( opts[i] == "-s" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -s is missing.");
                        seedkmax = atoi(opts[i+1].c_str());
                        i += 2;
                }
                else if ( opts[i] == "-e" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -e is missing.");
                        totalkmax = atoi(opts[i+1].c_str());
                        
                        if ( totalkmax > UniqueMatchInfoBase::getMaxErrors() )
                        {
                                totalkmax = UniqueMatchInfoBase::getMaxErrors();
                                std::cerr << "Warning: reducing maximum amount of errors to " << totalkmax << std::endl;
                        }
                        
                        i += 2;
                }
                else if ( opts[i] == "-l" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -l is missing.");
                        seedl = atoi(opts[i+1].c_str());
                        i += 2;
                }
#if defined(ALLOW_MATCH_ALL)
                else if ( opts[i] == "-u" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -u is missing.");
                        match_unique = atoi(opts[i+1].c_str());
                        i += 2;
                }
#endif
                else if ( opts[i] == "-g" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -g is missing.");
                        gaps = atoi(opts[i+1].c_str());
                        i += 2;
                }
                else if ( opts[i] == "-R" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -R is missing.");
                        rewritepatterns = atoi(opts[i+1].c_str());
                        i += 2;
                }
                else if ( opts[i] == "-m" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -m is missing.");
                        fracmem = atof(opts[i+1].c_str());
                        i += 2;
                }
                else if ( opts[i] == "-q" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -q is missing.");
                        scores = atoi(opts[i+1].c_str());
                        i += 2;
                }
                else if ( opts[i] == "-Q" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -Q is missing.");
                        qualityOffset = atoi(opts[i+1].c_str());
                        i += 2;
                }
                else if ( opts[i] == "-f" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -f is missing.");
                        fracmem = atof(opts[i+1].c_str());
                        i += 2;
                }
#if defined(_OPENMP)
                else if ( opts[i] == "-T" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -T is missing.");
                        int const numthreads = atoi(opts[i+1].c_str());
                        
                        if ( numthreads < 1 )
                                throw std::runtime_error("Argument for -T parameter is invalid (<1)");
                        
                        omp_set_num_threads(numthreads);
                        
                        i += 2;
                }
#endif
                else if ( opts[i] == "-similarity" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -similarity is missing.");
                        similarity = atof(opts[i+1].c_str());
                        
                        if ( similarity < 0 )
                        {
                                std::cerr << "Warning: setting similarity up to 0." << std::endl;
                                similarity = 0;
                        }
                        if ( similarity > 1 )
                        {
                                std::cerr << "Warning: setting similarity down to 1." << std::endl;
                                similarity = 1;
                        }
                        if ( similarity == 0 )
                        {
                                std::cerr << "Warning: similarity value is " << similarity << std::endl;
                        }
                        
                        i += 2;
                }
                else if ( opts[i] == "-err" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -err is missing.");
                        err = atof(opts[i+1].c_str());

                        if ( err < 0 )
                        {
                                std::cerr << "Warning: setting err up to 0." << std::endl;
                                err = 0;
                        }
                        if ( err > 1 )
                        {
                                std::cerr << "Warning: setting err down to 1." << std::endl;
                                err = 1;
                        }
                        if ( err == 1 )
                        {
                                std::cerr << "Warning: err value is " << err << std::endl;
                        }

                        i += 2;
                }
                else if ( opts[i] == "-trans" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -trans is missing.");
                        trans = atof(opts[i+1].c_str());

                        if ( trans < 0 )
                        {
                                std::cerr << "Warning: setting trans up to 0." << std::endl;
                                trans = 0;
                        }
                        if ( trans > 1 )
                        {
                                std::cerr << "Warning: setting trans down to 1." << std::endl;
                                trans = 1;
                        }
                        if ( trans == 0 || trans == 1 )
                        {
                                std::cerr << "Warning: trans value is " << trans << std::endl;
                        }

                        i += 2;
                }
                else if ( opts[i] == "-gc" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -gc is missing.");
                        gc = atof(opts[i+1].c_str());

                        if ( gc < 0 )
                        {
                                std::cerr << "Warning: setting gc up to 0." << std::endl;
                                gc = 0;
                        }
                        if ( gc > 1 )
                        {
                                std::cerr << "Warning: setting gc down to 1." << std::endl;
                                gc = 1;
                        }
                        if ( gc == 0 || gc == 1 )
                        {
                                std::cerr << "Warning: gc value is " << gc << std::endl;
                        }

                        i += 2;
                }
                else if ( opts[i] == "-gcmut_bias" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -gcmut_bias is missing.");
                        gcmut_bias = atof(opts[i+1].c_str());

                        if ( gcmut_bias < 0 )
                        {
                                std::cerr << "Warning: setting gcmut_bias up to 0." << std::endl;
                                gcmut_bias = 0;
                        }
                        if ( gcmut_bias == 0 )
                        {
                                std::cerr << "Warning: gcmut_bias value is " << gcmut_bias << std::endl;
                        }

                        i += 2;
                }
                else if ( opts[i] == "-filter_level" )
                {
                        if ( i+1 >= opts.size() )
                                throw std::runtime_error("Parameter for argument -filter_level is missing.");
                        filter_level = atoi(opts[i+1].c_str());

                        if ( filter_level < 0 )
                        {
                                std::cerr << "Warning: setting filter_level up to 0." << std::endl;
                                filter_level = 0;
                        }
                        if ( filter_level > 4 )
                        {
                                std::cerr << "Warning: setting filter_level down to 4." << std::endl;
                                filter_level = 4;
                        }

                        i += 2;
                }
                else if ( opts[i] == "-h" )
                {
                        printHelp();
                        throw std::runtime_error("Help requested.");
                }
                else
                {
                        std::cerr << "Ignoring argument " << opts[i] << std::endl;
                        i += 1;
                }
        }
        
        bool const optionscomplete = textfilename.size() && patternfilename.size() && outputfilename.size();
        
        if ( ! optionscomplete )
                printHelp();        
        
        if ( ! textfilename.size() )
                throw std::runtime_error("Mandatory argument -t (text file name) is not given.");
        if ( ! patternfilename.size() )
                throw std::runtime_error("Mandatory argument -p (pattern file name) is not given.");
        if ( ! outputfilename.size() )
                throw std::runtime_error("Mandatory argument -o (output file name) is not given.");

        fracmem = std::min ( 1.0, fracmem );

        std::cerr << "fraction of memory used will be " << fracmem << std::endl;

        size_t physmem = getPhysicalMemory();
        usemem = static_cast<size_t>(physmem * fracmem);
        std::cerr << "detected " << physmem/(1024*1024) << " mega bytes of physical memory, of which we will use " << usemem/(1024*1024) << "(" << usemem << ")" << std::endl;

        if ( patternfilename == "-" )
        {
        	fastq = isFastQ(std::cin);
        	if ( ! rewritepatterns )
        	{
        		std::cerr << "Reading patterns from stdin, switching on pattern rewriting." << std::endl;
	        	rewritepatterns = true;
		}
	}
        else
        {
	        fastq = isFastQ(patternfilename);
	}

        std::cerr << "pattern file is " << (fastq ? "FASTQ" : "FASTA") << " rewrite is " << (rewritepatterns?"on":"off") << std::endl;

        if ( seedl > 64 )
        {
                seedl = 64;
                std::cerr << "reduced seed size to " << seedl << " to not exceed 64." << std::endl;
        }
        if ( seedl % 4 )
        {
                seedl -= (seedl%4);
                std::cerr << "reduced seed size to " << seedl << " to have a multiple of 4." << std::endl;
        }
        if ( seedl < static_cast<int>(nu) )
        {
                throw std::runtime_error("cannot handle seed length < 4");
        }

        if ( seedkmax > 2 )
        {
                seedkmax = 2;
                std::cerr << "reduced number of mismatches in seed to " << seedkmax << " as we cannot handle more." << std::endl;
        }
        
        switch ( filter_level )
        {
        	case 1: filter_mult = 0.5 * totalkmax; break;
        	case 2: filter_mult = 1 * totalkmax; break;
        	case 3: filter_mult = 2 * totalkmax; break;
        	case 4: filter_mult = 3 * totalkmax; break;
        	case 0: default: filter_mult = 0 * totalkmax; break;
        }
        filter_mult /= 70.0;
        
        std::cerr << "filter_mult=" << filter_mult << std::endl;
}        

double const RealOptions::default_fracmem = 0.75;
