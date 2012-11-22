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

#include <iostream>
#include <cstdlib>
#include <ctime>

int main(int argc, char * argv[])
{
	if ( argc < 2 )
	{
		std::cerr << "usage: " << argv[0] << " <length>" << std::endl;
		return EXIT_FAILURE;
	}

	int len = atoi(argv[1]);
	srand(time(0));

	std::cout << "> random_" << len << "\n";	
	for ( unsigned int i = 0; i < len; ++i )
	{
		if ( i && (i % 60 == 0) )
			std::cout.put('\n');
			
		switch ( rand() % 4 )
		{
			case 0: std::cout.put('A'); break;
			case 1: std::cout.put('C'); break;
			case 2: std::cout.put('G'); break;
			case 3: std::cout.put('T'); break;
		}
	}
	std::cout << std::endl;
}