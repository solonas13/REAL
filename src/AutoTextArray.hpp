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

#if ! defined(AUTOTEXTARRAY_HPP)
#define AUTOTEXTARRAY_HPP

#include "AutoArray.hpp"
#include "ERank222B.hpp"
#include "acgtnMap.hpp"
#include <iostream>

template<typename iterator>
AutoArray<u_int64_t> getTextArray(iterator utext, u_int64_t const length)
{
	u_int64_t const numbits = length*2;
	u_int64_t const numbitsr = ((numbits+63)/64)*64;
	u_int64_t const numwords = numbitsr/64;

	typedef Rank::FastWriteBitWriter8 writer_type;
	AutoArray<writer_type::data_type> atext(numwords);
	writer_type writer(atext.get());
	
	for ( unsigned int i = 0; i < length; ++i )
		writer.write(utext[i] & 0x3, 2);
	writer.flush();

	return atext;
}

template<typename iterator>
AutoArray<u_int64_t> getWildcardArray(iterator utext, u_int64_t const length)
{
	u_int64_t const numbits = length;
	u_int64_t const numbitsr = ((numbits+63)/64)*64;
	u_int64_t const numwords = numbitsr/64;

	typedef Rank::FastWriteBitWriter8 writer_type;
	AutoArray<writer_type::data_type> atext(numwords);
	writer_type writer(atext.get());
	
	for ( unsigned int i = 0; i < length; ++i )
		writer.writeBit(utext[i] > 3);
	writer.flush();

	return atext;
}

template<bool sse4>
struct AutoTextArray
{
        private:
	typedef Rank::ERank222B<sse4> eclass;
	
        u_int64_t const n;
	AutoArray<u_int64_t> const atext;
	AutoArray<u_int64_t> const awild;
	u_int64_t const * const text;
	u_int64_t const * const wild;
	eclass wildrank;
	
	AutoTextArray & operator=(AutoTextArray const &)
	{
	        return *this;
	}
	
	public:
	u_int64_t getTextWord(unsigned int i) const
	{
	        assert ( i < atext.getN() );
	        return text[i];
	}
	
	u_int64_t getN() const
	{
	        return n;
	}
	
	u_int64_t size() const
	{
	        return atext.size() + awild.size() + 2*sizeof(u_int64_t *) + wildrank.size();
	}
	
	template<typename iterator>
	AutoTextArray(iterator text, u_int64_t const rn)
	: n(rn), atext(getTextArray(text,n)), awild(getWildcardArray(text,n)), text(atext.get()), wild(awild.get()),
	  wildrank(awild.get(), ((n+63)/64)*64 )
	{
	        #if 1
	        std::cerr << "Checking AutoTextArray...";
  	        for ( unsigned int i = 0; i < n; ++i )
  		        assert ( text[i] == (*this)[i] );
		std::cerr << "done." << std::endl;
                #endif
	}
	
	~AutoTextArray()
	{
	        // std::cerr << "Deallocating AutoTextArray." << std::endl;
	}
	
	u_int64_t operator[](u_int64_t i) const
	{
	        assert ( i < n );
		return (static_cast<unsigned int>(Rank::getBit(wild,i)) << 2) | Rank::getBits64(text,i<<1,2);
	}
	
	u_int64_t getTextWord(u_int64_t const i, u_int64_t const l) const
	{
                return Rank::getBits64(text,i<<1,l<<1);
	}
	u_int64_t getGappedTextWord(u_int64_t const i, u_int64_t const p, u_int64_t const g, u_int64_t const q) const
	{
		return ( getTextWord(i,p) << (2*q) ) | getTextWord(i+p+g,q);
	}
	u_int64_t countDontCareFreeGappedTextWords(u_int64_t const p, u_int64_t const g, u_int64_t const q) const
	{
		u_int64_t l = p+g+q;

		if ( (n+1) < l )
			return 0;
		
		u_int64_t c = 0;
		for ( u_int64_t i = 0; i < (n+1)-l; ++i )
			if ( isDontCareFree(i,p) && isDontCareFree(i+p+q,q) )
				c++;

		return c;
	}
	u_int64_t fillDontCareFreeGappedTextWordPositions(
		u_int32_t * const A,
		u_int64_t const s,
		u_int64_t & ii,
		u_int64_t const p,
		u_int64_t const g,
		u_int64_t const q
		) const
	{
		u_int64_t l = p+g+q;

		if ( (n+1) < l )
			return 0;
		
		u_int64_t c = 0;
		u_int64_t i;
		for ( i = ii ; (i < (n+1)-l) && (c < s); ++i )
			if ( isDontCareFree(i,p) && isDontCareFree(i+p+q,q) )
				A[c++] = i;

		ii = i;
		return c;
	}
	bool isDontCareFree(u_int64_t const i, u_int64_t const l) const
	{
                u_int64_t const bef = i ? wildrank.rank1(i-1) : 0;
                u_int64_t const aft = wildrank.rank1(i+l-1);
                return (aft-bef) == 0;
	}
	static std::string wordToString(u_int64_t const word, unsigned int const l)
	{
	        unsigned int shift = ((l-1)<<1);
	        u_int64_t mask = static_cast<u_int64_t>(0x3ull) << shift;
	        
	        std::string s(l,' ');
	        for ( unsigned int i = 0; i < l; ++i )
	        {
                        s[i] = (word & mask) >> shift;
	                mask >>= 2;
	                shift -= 2;
	        }
	        
	        return s;
	}
};

template<bool sse4>
struct AutoTextArrayGappedComparator
{
	AutoTextArray<sse4> const & ATA;
	u_int64_t const p;
	u_int64_t const g;
	u_int64_t const q;
	
	AutoTextArrayGappedComparator(AutoTextArray<sse4> const & rATA, u_int64_t const rp, u_int64_t const rg, u_int64_t const rq)
	: ATA(rATA), p(rp), g(rg), q(rq)
	{
		
	}
	
	bool operator()(u_int32_t const i, u_int32_t const j) const
	{
	        // std::cerr << "Comparing " << i << " and " << j << std::endl;
		return ATA.getGappedTextWord(i,p,g,q) < ATA.getGappedTextWord(j,p,g,q);
	}
	
	u_int64_t operator()(u_int32_t const i) const
	{
	        return ATA.getGappedTextWord(i,p,g,q);
	}
	std::string getString(u_int32_t const i) const
	{
	        // std::cerr << "p=" << p << " q=" << q << std::endl;
	
	        std::string s = ATA.wordToString ( (*this)(i) , p + q );
	        
	        for ( unsigned int j = 0; j < s.size(); ++j )
	                s[j] = toollib::remapChar(s[j]);
	                
                // std::cerr << "word=" << (*this)(i) << std::endl;
                
                return s;
	}
};
#endif
