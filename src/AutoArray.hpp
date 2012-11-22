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

#if ! defined(AUTOARRAY_HPP)
#define AUTOARRAY_HPP

#include <cstdlib>
#include <vector>

template<typename N>
struct AutoArray
{
	protected:
	mutable N * array;
	mutable size_t n;

	public:
	unsigned int size() const
	{
		return  sizeof(size_t) + sizeof(N *) + n * sizeof(N);
	}
	unsigned int getN() const
	{
		return n;
	}
	
	AutoArray() : array(0), n(0) {}
	AutoArray(AutoArray const & o) : array(0) {
		array = o.array;
		n = o.n;
		o.array = 0;
		o.n = 0;
	}
	AutoArray(size_t rn, bool erase = true) : array(0), n(rn) {
		array = new N[n];

		if ( erase )
			for ( size_t i = 0; i < n; ++i )
				array[i] = N();
	}
	AutoArray(size_t rn, N const * c) : array(0), n(rn) 
	{
		array = new N[n];
		for ( unsigned int i = 0; i < rn; ++i )
			array[i] = c[i];
	}
	virtual ~AutoArray() { 
		delete [] array; 
	}
	N * get() { return array; }
	N const * get() const { return array; }
	
	N       & operator[](size_t i)       { return array[i]; }
	N const & operator[](size_t i) const { return array[i]; }
	
	AutoArray<N> & operator=(AutoArray<N> const & o)
	{
		if ( this != &o )
		{
			delete [] array;
			this->array = 0;
			this->n = 0;
			this->array = o.array;
			this->n = o.n;
			o.array = 0;
			o.n = 0;
		}
		return *this;
	}
	
	void release()
	{
		delete [] array;
		array = 0;
		n = 0;
	}

	AutoArray<N> clone() const
	{
		AutoArray c(n,false);
		N * d = c.get();
		for ( size_t i = 0; i < n; ++i )
			d[i] = array[i];
		return c;
	}
	std::vector<N> toVector() const
	{
		std::vector<N> V;
		for ( size_t i = 0; i < n; ++i ) V.push_back(array[i]);
		return V;
	}
};
#endif
