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

#include "getPhysicalMemory.hpp"

#if defined(HAVE_WINDOWS_H)
#include <windows.h>
#endif

#if defined(HAVE_WINDOWS_H)
static size_t getMemoryWindowsApi()
{
	MEMORYSTATUS memstat;
	GlobalMemoryStatus(&memstat);
	return memstat.dwTotalPhys;
}
#endif

#include <fstream>
#include <sstream>
#include <iostream>

static size_t getMemorySizeProc()
{
	std::ifstream meminfo("/proc/meminfo");
	
	if ( ! meminfo.is_open() )
		return 0;
	
	while ( meminfo )
	{
		std::string line; std::getline(meminfo,line);
		
		if ( meminfo )
		{
			std::string const tag = "MemTotal:";
		
			if ( 
				line.size() >= tag.size()
				&&
				line.substr(0,tag.size()) == tag 
			)
			{
				line = line.substr(tag.size());
				while ( line.size() && isspace(line[0]) )
					line = line.substr(1);

				if ( line.size() >= 2 && line.substr(line.size()-2) == "kB" )
				{
					line = line.substr(0,line.size()-2);
					while ( line.size() && isspace(line[line.size()-1]) )
						line = line.substr(0,line.size()-1);
						
					bool allnum = true;
					for ( unsigned int i = 0; i < line.size(); ++i )
						if ( ! isdigit(line[i]) )
							allnum = false;

					if ( allnum && line.size() )
					{
						size_t mem;
						std::istringstream istr(line);
						istr >> mem;
						
						mem *= 1024;
						
						return mem;
					}
				}
			}
		}
	}	
	
	return 0;
}

size_t getPhysicalMemory()
{
	size_t physmem = 0;
	
	#if defined(HAVE_WINDOWS_H)
	physmem = getMemoryWindowsApi();
	#endif
	
	if ( ! physmem )
		physmem = getMemorySizeProc();
	
	
	if ( ! physmem )
	{
		std::cerr << "Failed to detect amount of physical memory, using 1GB." << std::endl;
		physmem = 1024*1024*1024ull;
	}

	return physmem;
}
