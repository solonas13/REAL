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

#if ! defined(REALOPTIONS_HPP)
#define REALOPTIONS_HPP

#include <string>
#include "types.hpp"

struct RealOptions
{
        static unsigned int const default_seedkmax = 2;
        static unsigned int const default_totalkmax = 5;
        static unsigned int const default_seedl = 32;
        static bool const default_match_unique = true;
        static double const default_fracmem;
        static bool const default_scores = true;
        static bool const default_rewritepatterns = true;
        static unsigned int const default_sort_threads = 2;
        static unsigned int const default_filter_level = 2;
        static bool const default_gaps = 0;

        static unsigned int const nu = 4; // number of segments to split signature in
        std::string textfilename;
        std::string patternfilename;
        std::string outputfilename;
        unsigned int seedkmax;
        unsigned int totalkmax;
        int seedl;
        bool match_unique;
        double fracmem;
        bool scores;
        unsigned int qualityOffset;
        u_int64_t usemem;
        bool rewritepatterns;
        unsigned int sort_threads;
        int filter_level;
        double filter_mult;

        /* scoring parameters */
        double similarity;
        double err;
        double trans;
        double gc;
        double gcmut_bias;

        bool gaps;
        
        bool fastq;
                
        public:

        void printHelp();
        static bool isFastQ(std::string const &);
        static bool isFastQ(std::istream &);

        RealOptions(int argc, char * argv[]);
        
        inline double getFilterValue(unsigned int const patl) const
        {
        	return filter_mult * patl;
        }
};
#endif
