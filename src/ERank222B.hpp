/**
    REAL: An efficient REad ALigner for next generation sequencing reads.
    Copyright (C) 2010 German Tischler

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

#if ! defined(ERANK222B_HPP)
#define ERANK222B_HPP

#include "AutoArray.hpp"
#include "PopCountTable.hpp"
#include "types.hpp"
#include <cassert>
#include <stdexcept>
#include <limits>

namespace Rank
{
	template<unsigned int k>
	struct MetaLog2
	{
		static unsigned int const log = 1 + MetaLog2<k/2>::log;
	};
	template<>
	struct MetaLog2<1>
	{
		static unsigned int const log = 0;
	};

	template<typename N>
	inline unsigned int numbits(N n)
	{
		unsigned int c = 0;
		while ( n )
		{
			c++;
			n >>= 1;
		}
		return c;
	}

	// get bits from an array
	template<typename iterator>
	inline u_int64_t getBits64(iterator A, u_int64_t offset, unsigned int numbits)
	{
		if ( ! numbits )
			return 0;

		typedef typename std::iterator_traits<iterator>::value_type value_type;
		unsigned int const bcnt = 8 * sizeof(value_type);
		unsigned int const bshf = MetaLog2<bcnt>::log;
		unsigned int const bmsk = (1ul<<bshf)-1ul;

		unsigned int byteSkip = (offset >> bshf);
		unsigned int const bitSkip = (offset & bmsk);
		unsigned int const restBits = bcnt-bitSkip;
			
		u_int64_t b = static_cast<u_int64_t>(A[byteSkip]);
		// skip bits by masking them
		b &= std::numeric_limits<value_type>::max() >> (bcnt - restBits);
		
		if ( restBits == numbits )
			return b;
		else if ( numbits < restBits )
			return b >> (restBits - numbits);
			
		numbits -= restBits;
		
		if ( numbits )
			b = (b<<numbits) | (static_cast<u_int64_t>( A[++byteSkip] ) >> (bcnt-numbits));
		
		return b;
	}


	/**
	 * read one bit from a bit vector
	 * @param A bit vector
	 * @param offset
	 * @return bit
	 **/
	template<typename iterator>
	inline bool getBit1(iterator A, u_int64_t const offset)
	{
		static const u_int8_t maskone[] = { (1u<<7), (1u<<6), (1u<<5), (1u<<4), (1u<<3), (1u<<2), (1u<<1), (1u<<0) };
		return A[(offset >> 3)] & maskone[offset & 0x7u];
	}
	/**
	 * read one bit from a bit vector
	 * @param A bit vector
	 * @param offset
	 * @return bit
	 **/
	template<typename iterator>
	inline bool getBit2(iterator A, u_int64_t const offset)
	{
		static const 
			u_int16_t maskone[] = { 
				(1u<<15), (1u<<14), (1u<<13), (1u<<12), (1u<<11), (1u<<10), (1u<<9), (1u<<8),
				(1u<<7), (1u<<6), (1u<<5), (1u<<4), (1u<<3), (1u<<2), (1u<<1), (1u<<0) 
			};

		// std::cerr << "getBit2::Offset: " << (offset >> 4) << " num " << A[offset>>4] << std::endl;
		
		return A[(offset >> 4)] & maskone[offset & 0xfu];
	}
	/**
	 * read one bit from a bit vector
	 * @param A bit vector
	 * @param offset
	 * @return bit
	 **/
	template<typename iterator>
	inline bool getBit4(iterator A, u_int64_t const offset)
	{
		static const 
			u_int32_t maskone[] = { 
				(1u<<31), (1u<<30), (1u<<29), (1u<<28), (1u<<27), (1u<<26), (1u<<25), (1u<<24),
				(1u<<23), (1u<<22), (1u<<21), (1u<<20), (1u<<19), (1u<<18), (1u<<17), (1u<<16),
				(1u<<15), (1u<<14), (1u<<13), (1u<<12), (1u<<11), (1u<<10), (1u<<9), (1u<<8),
				(1u<<7), (1u<<6), (1u<<5), (1u<<4), (1u<<3), (1u<<2), (1u<<1), (1u<<0) 
			};
		return A[(offset >> 5)] & maskone[offset & 0x1fu];
	}
	/**
	 * read one bit from a bit vector
	 * @param A bit vector
	 * @param offset
	 * @return bit
	 **/
	template<typename iterator>
	inline bool getBit8(iterator A, u_int64_t const offset)
	{
		static const 
			u_int64_t maskone[] = { 
				(1ULL<<63), (1ULL<<62), (1ULL<<61), (1ULL<<60), (1ULL<<59), (1ULL<<58), (1ULL<<57), (1ULL<<56),
				(1ULL<<55), (1ULL<<54), (1ULL<<53), (1ULL<<52), (1ULL<<51), (1ULL<<50), (1ULL<<49), (1ULL<<48),
				(1ULL<<47), (1ULL<<46), (1ULL<<45), (1ULL<<44), (1ULL<<43), (1ULL<<42), (1ULL<<41), (1ULL<<40),
				(1ULL<<39), (1ULL<<38), (1ULL<<37), (1ULL<<36), (1ULL<<35), (1ULL<<34), (1ULL<<33), (1ULL<<32),
				(1ULL<<31), (1ULL<<30), (1ULL<<29), (1ULL<<28), (1ULL<<27), (1ULL<<26), (1ULL<<25), (1ULL<<24),
				(1ULL<<23), (1ULL<<22), (1ULL<<21), (1ULL<<20), (1ULL<<19), (1ULL<<18), (1ULL<<17), (1ULL<<16),
				(1ULL<<15), (1ULL<<14), (1ULL<<13), (1ULL<<12), (1ULL<<11), (1ULL<<10), (1ULL<<9), (1ULL<<8),
				(1ULL<<7), (1ULL<<6), (1ULL<<5), (1ULL<<4), (1ULL<<3), (1ULL<<2), (1ULL<<1), (1ULL<<0) 
			};
		return A[(offset >> 6)] & maskone[offset & 0x3Fu];
	}
	/**
	 * read one bit from a bit vector
	 * @param A bit vector
	 * @param offset
	 * @return bit
	 **/
	template<typename iterator>
	inline bool getBit(iterator A, u_int64_t const offset)
	{
		typedef typename std::iterator_traits<iterator>::value_type value_type;
		
		switch( sizeof(value_type) )
		{
			case 8:          return getBit8(A,offset);
			case 4:          return getBit4(A,offset);
			case 2:          return getBit2(A,offset);
			case 1: default: return getBit1(A,offset);
		}
	}


	/**
	 * bit stream writer class
	 **/
	template<typename _data_type, typename _data_iterator, _data_type basemask>
	struct BitWriterTemplate
	{
		public:
		/**
		 * data type used by this bit writer class
		 **/
		typedef _data_type data_type;
		typedef _data_iterator data_iterator;
	
		private:
		data_iterator U;
		data_type mask;
		data_type cur;
		
		public:
		/**
		 * initialize writer with pointer to array
		 **/
		BitWriterTemplate(data_iterator rU) : U(rU), mask(basemask), cur(0) {}
		
		void writeUnary(u_int64_t k)
		{
			for ( u_int64_t i = 0; i < k; ++i )
				writeBit(0);
			writeBit(1);
		}
		
		/**
		 *
		 **/
		template<typename N>
		void writeElias2(N n)
		{
			// number of bits to store n
			unsigned int log_1 = numbits(n);
			// number of bits to store log_1
			unsigned int log_2 = numbits(log_1);
			
			// write log_2 in unary form
			writeUnary(log_2);
			
			// write log_1 using log_2 bits
			write(log_1,log_2);
			
			// write n using log_1 bits
			write(n,log_1);
		}
		
		/**
		 * write a b bit number n
		 * @param n number to be written
		 * @param b number of bits to write
		 **/
		template<typename N>
		void write(N n, unsigned int b)
		{
			if ( b )
			{
				N m = static_cast<N>(1) << (b-1);
					
				// write number, msb to lsb
				for ( unsigned int i = 0; i < b; ++i, m >>= 1 )
				{
					writeBit( (n & m) != 0 );
#if defined(RANK_WRITE_DEBUG)
					n &= ~m;
#endif
				}
#if defined(RANK_WRITE_DEBUG)
				assert ( n == 0 );
#endif
			}
		}
		/**
		 * write one bit to stream
		 * @param bit
		 **/
		void writeBit(bool const bit)
		{
			if ( bit ) 
			{
				cur |= mask;
			}
			
			mask >>= 1;
			
			if ( ! mask )
			{
				*(U++) = cur;
				mask = basemask;
				cur = 0;
			}
		}
		/**
		 * flush output (align to byte boundary) by writing zero bits
		 **/
		void flush()
		{
			while ( mask != basemask )
				writeBit(0);
		}
	};
	
	typedef BitWriterTemplate<u_int64_t, u_int64_t *, 0x8000000000000000ULL> BitWriter8;

	/**
	 * bit stream writer class
	 **/
	template<typename _data_type, typename _data_iterator, _data_type basemask, _data_type fullmask>
	struct FastWriteBitWriterTemplate
	{
		public:
		/**
		 * data type used by this bit writer class
		 **/
		typedef _data_type data_type;
		typedef _data_iterator data_iterator;
	
		private:
		data_iterator U;
		data_type mask;
		data_type cur;
		unsigned int bitsleft;
		
		public:
		/**
		 * initialize writer with pointer to array
		 **/
		FastWriteBitWriterTemplate(data_iterator rU) : U(rU), mask(basemask), cur(0), bitsleft(8 * sizeof(data_type)) {}
		
		void writeUnary(u_int64_t k)
		{
			for ( u_int64_t i = 0; i < k; ++i )
				writeBit(0);
			writeBit(1);
		}
		
		/**
		 *
		 **/
		template<typename N>
		void writeElias2(N n)
		{
			// number of bits to store n
			unsigned int log_1 = numbits(n);
			// number of bits to store log_1
			unsigned int log_2 = numbits(log_1);
			
			// write log_2 in unary form
			writeUnary(log_2);
			
			// write log_1 using log_2 bits
			write(log_1,log_2);
			
			// write n using log_1 bits
			write(n,log_1);
		}
		
		/**
		 * write a b bit number n
		 * @param n number to be written
		 * @param b number of bits to write
		 **/
		template<typename N>
		void write(N n, unsigned int b)
		{
			if ( b < bitsleft )
			{
				cur |= static_cast<data_type>(n) << (bitsleft - b);
				bitsleft -= b;
				mask >>= b;
			}
			else
			{
				b -= bitsleft;
				*(U++) = cur | (n >> b);

				if ( sizeof(N) > sizeof(data_type) )
					while ( b >= 8 * sizeof(data_type) )
					{
						b -= 8*sizeof(data_type);
						*(U++) = (n >> b) & fullmask;
					}
				
				n &= ((static_cast<N>(1)<<b)-1);

				cur = static_cast<data_type>(n) << (8*sizeof(data_type) - b);
				mask = basemask >> b;
				bitsleft = 8 * sizeof(data_type) - b;
			}
		}
		
		/**
		 * write one bit to stream
		 * @param bit
		 **/
		void writeBit(bool const bit)
		{
			if ( bit ) 
			{
				cur |= mask;
			}
			
			mask >>= 1;
			bitsleft -= 1;
			
			if ( ! mask )
			{
				*(U++) = cur;
				mask = basemask;
				cur = 0;
				bitsleft = 8 * sizeof(data_type);
			}
		}
		/**
		 * flush output (align to byte boundary) by writing zero bits
		 **/
		void flush()
		{
			while ( mask != basemask )
				writeBit(0);
		}
	};

	typedef FastWriteBitWriterTemplate<u_int64_t, u_int64_t *, 0x8000000000000000ULL, 0xFFFFFFFFFFFFFFFFULL> FastWriteBitWriter8;

	/**
	 * two level rank dictionary using approximately n/4 bits for index.
	 * the superblock size is 2^16 bits, the miniblock size is 64 bits.
	 * rank queries are answered using the dictionary and a 64 bit
	 * population count function. if the machine instruction set
	 * does not provide a 64 bit popcount function, these calls
	 * are simulated by using a precomputed 16 bit lookup table.
	 **/
	template<bool sse4>
	struct ERank222B
	{
		public:
		typedef Rank::BitWriter8 writer_type;
	
		private:
		// super block size 2^16 bits
		static unsigned int const sbbitwidth = 16;
		// mini block size 2^6 = 64 bits
		static unsigned int const mbbitwidth = 6;
		
		// derived numbers
		
		// actual block sizes
		static u_int64_t const sbsize = 1u << sbbitwidth;
		static u_int64_t const mbsize = 1u << mbbitwidth;
		
		// miniblocks per superblock
		static u_int64_t const mps = (sbsize + mbsize-1)/mbsize;
		// superblock mask
		static u_int64_t const sbmask = (1u<<sbbitwidth)-1;
		// miniblock mask
		static u_int64_t const mbmask = (1u<<mbbitwidth)-1;

		u_int64_t const * const UUUUUUUU;
		u_int64_t const n;
		u_int64_t const numsuper;
		u_int64_t const nummini;
		
		AutoArray<u_int64_t> S; // n / 2^16 * 64 bits = n / 2^10 = n/1024 bits
		AutoArray<unsigned short> M; // n / 2^16 * 2^16 / 64 * 16 = n/4 bits

		static inline u_int64_t divUp(u_int64_t a, u_int64_t b)
		{
			return (a+b-1)/b;
		}

		/**
		 * return superblock containing i th 1 bit,
		 * i.e. largest j such that S[j] < ii
		 **/
		u_int64_t selectSuper(u_int64_t const ii) const
		{
			// search largest superblock index s such that ii < S[s]
			u_int64_t left = 0, right = numsuper;

			while ( right-left > 1 )
			{
				u_int64_t const d = right-left;
				u_int64_t const d2 = d>>1;
				u_int64_t const mid = left + d2;

				// number of 1s is too large
				if ( S[mid] < ii )
					left = mid;
				else
					right = mid;
			}

			return left;
		}
		/**
		 * return miniblock containing i th 1 bit,
		 * i.e. largest j such that M[j] < ii
		 **/
		u_int64_t selectMini(u_int64_t const s, u_int64_t const iii) const
		{
			u_int64_t const ii = iii - S[s];
			u_int64_t left = (s << sbbitwidth) >>  mbbitwidth;
			u_int64_t right = std::min( nummini, ((s+1) << sbbitwidth) >>  mbbitwidth);
		
			while ( right-left > 1 )
			{
				u_int64_t const d = right-left;
				u_int64_t const d2 = d>>1;
				u_int64_t const mid = left + d2;

				// number of 1s is too large
				if ( M[mid] < ii )
					left = mid;
				else
					right = mid;
			}

			return left;
		}

		
		public:		
		/**
		 * @param rUUUUUUUU bit vector
		 * @param rn number of bits in vector (has to be a multiple of 64)
		 **/
		ERank222B(u_int64_t const * const rUUUUUUUU, u_int64_t const rn) 
		: UUUUUUUU(rUUUUUUUU), n(rn),
		  numsuper((n + (sbsize-1)) >> sbbitwidth), nummini((n + (mbsize-1)) >> mbbitwidth),
		  S( divUp(n,sbsize) , false ), M( divUp(n,mbsize), false)
		{
			if ( n & mbmask )
				throw std::runtime_error("Rank::ERank222B: n is not multiple of miniblock size 64.");
		
			u_int64_t c = 0;

			// superblock counter
			int64_t s = -1;
			// miniblock counter
			int64_t m = -1;

			u_int64_t i = 0;
			for ( u_int64_t mi = 0; mi < nummini; ++mi, i += mbsize )
			{
				u_int64_t b = UUUUUUUU[mi];

				if ( (i & sbmask) == 0 )
				{
					S[ ++s ] = c;
					assert ( S[s] == c);
				}

				M[ ++m ] = c - S[s];
				assert( S[s] + M[m] == c );
			
				c += toollib::PopCount<sse4>::popcount8(b);
			}
		}
		
		/**
		 * @return estimated space in bytes
		 **/		
		u_int64_t size() const
		{
			return 
				sizeof(u_int64_t *) + 
				3*sizeof(u_int64_t) +
				S.size() + 
				M.size();
		}
		
		/**
		 * return number of 1 bits up to (and including) index i
		 * @param i
		 * @return population count
		 **/
		u_int64_t rank1(u_int64_t i) const
		{
			u_int64_t const mi = i >> mbbitwidth;
			return S[i >> sbbitwidth] + M[mi] + toollib::PopCount<sse4>::popcount8( (UUUUUUUU[mi]), i - (mi<<mbbitwidth));
		}
		/**
		 * return number of 0 bits up to (and including) index i
		 * @param i
		 * @return inverse population count
		 **/
		u_int64_t rank0(u_int64_t i) const
		{
			return (i+1) - rank1(i);
		}
		/**
		 * Return the position of the ii'th 1 bit. This function is implemented using a 
		 * binary search on the rank1 function.
		 **/
		u_int64_t select1(u_int64_t const ii) const
		{
			u_int64_t i = ii+1;

			u_int64_t const s = selectSuper(i);
			u_int64_t const m = selectMini(s,i);
			i -= S[s]; i -= M[m];
			u_int64_t const v = (UUUUUUUU[m]);
			
			u_int64_t left = 0, right = 1u<<mbbitwidth;
			while ( right-left )
			{
				u_int64_t const d = right-left;
				u_int64_t const d2 = d>>1;
				u_int64_t const mid = left + d2;
				u_int64_t const rmid = toollib::PopCount<sse4>::popcount8( v, mid );

				// number of ones is too small
				if ( rmid < i )
					left = mid+1;
				// number of ones is too large
				else if ( rmid > i )
					right = mid;
				// if this is the leftmost occurence in the interval, return it
				else if ( (!mid) || ( toollib::PopCount<sse4>::popcount8( v, mid-1 ) != i) )
				{
					return (m<<mbbitwidth)+mid;
				}
				// otherwise, go on and search to the left
				else
					right = mid;
			}
			
			return n;
		}
		/**
		 * Return the position of the ii'th 0 bit. This function is implemented using a 
		 * binary search on the rank1 function.
		 **/
		u_int64_t select0(u_int64_t const ii) const
		{
			u_int64_t const i = ii+1;

			u_int64_t left = 0, right = n;
			
			while ( (right-left) )
			{
				u_int64_t const d = right-left;
				u_int64_t const d2 = d>>1;
				u_int64_t const mid = left + d2;

				// number of ones is too small
				if ( rank0(mid) < i )
					left = mid+1;
				// number of ones is too large
				else if ( rank0(mid) > i )
					right = mid;
				// if this is the leftmost occurence in the interval, return it
				else if ( (!mid) || (rank0(mid-1) != i) )
					return mid;
				// otherwise, go on and search to the left
				else
					right = mid;
			}
			
			return n;		
		}
	};
}
#endif
