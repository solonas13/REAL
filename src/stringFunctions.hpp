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

#if ! defined(STRINGFUNCTIONS_HPP)
#define STRINGFUNCTIONS_HPP

#include <deque>

struct stringFunctions {
	template<typename string_type>
	static std::deque<string_type> tokenize(string_type input, string_type br) 
	{
		typename string_type::iterator a,b,c;
		std::deque<string_type> tokens;

		while ( input.length() )
		{
			a = input.begin();
			c = input.end();
			typename string_type::size_type e = input.find(br);

			if ( e != string_type::npos )
			{
				b = a+e;
				tokens.push_back(string_type(a,b));
				input = string_type(b+br.length(),c);
			}
			else
			{
				tokens.push_back(input);
				input = string_type();
			}
		}

		return tokens;
	}
	template<typename string_type>
	static bool endsOn(string_type const & name, string_type const & suffix)
	{
		if ( name.size() < suffix.size() )
			return false;
		for ( unsigned int i = 0; i < suffix.size(); ++i )
			if ( name[ name.size()-suffix.size()+i ] != suffix[i] )
				return false;
		return true;
	}

};
#endif
