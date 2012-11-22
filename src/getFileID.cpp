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

#include "getFileID.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

/**
 * obtain file id
 **/
std::string getFileID(std::string const filename)
{
        std::ifstream infile(filename.c_str());
        
        if ( ! infile.is_open() )
        {
        	std::cerr << "Failed to open file " << filename << std::endl;
        	throw std::runtime_error("Failed to open file.");
        }

        std::string line;
        std::getline(infile,line);
        
        if ( (! line.size()) || line[0] != '>' )
        {
                std::cerr << "File is malformed, leading comment missing." << std::endl;                
        	throw std::runtime_error("File is malformed.");
        }
        
        line = line.substr(1);
        
        return line;
}
