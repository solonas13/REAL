/*
Copyright (c) 2008, 2009 Andrew I. Schein
Copyright (c) 2010, German Tischler, Solon Pissis

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#if defined(_OPENMP)
#include <omp.h>
#endif
#include "Lock.hpp"

#include <iostream>
#include <cassert>

template<typename mask_type>
inline typename mask_type::signature_type u4_0(mask_type const & v) { return ((v.sign)         & 0x7FF); }
template<typename mask_type>
inline typename mask_type::signature_type u4_1(mask_type const & v) { return (((v.sign) >> 11) & 0x7FF); }
template<typename mask_type>
inline typename mask_type::signature_type u4_2(mask_type const & v) { return (((v.sign) >> 22) & 0x7FF); }

template<typename mask_type>
inline void u4_sort( AutoArray < mask_type > & Areader, size_t const n) 
{
	static unsigned int const HIST_SIZE = 2048;

	mask_type * reader = Areader.get();
	if (n < HIST_SIZE) { std::sort(reader,reader+n); return; }
	
	/* allocate 3 lists of HIST_SIZE elements */
	AutoArray< u_int64_t > Ab0(3 * HIST_SIZE, false);
	u_int64_t * const b0   = Ab0.get();
	u_int64_t * const b1   = b0 + HIST_SIZE;
	u_int64_t * const b2   = b1 + HIST_SIZE;

	/* erase lists */
	std::fill(b0, b0+3*HIST_SIZE, 0);

#if defined(_OPENMP)
	toollib::Lock lock;
#pragma omp parallel
	{
		size_t const packsize = (n + omp_get_num_threads()-1)/omp_get_num_threads();
		size_t const packlow = omp_get_thread_num()*packsize;
		size_t const packhigh = std::min((omp_get_thread_num()+1)*packsize,n);
		
		AutoArray< u_int64_t > localAb0(3 * HIST_SIZE, false);
		u_int64_t * const lb0   = localAb0.get();
		u_int64_t * const lb1   = lb0 + HIST_SIZE;
		u_int64_t * const lb2   = lb1 + HIST_SIZE;

		/* erase lists */
		for ( size_t i = 0; i < 3*HIST_SIZE; ++i )
			lb0[i] = 0;
		
		/* create histograms */
		for (size_t i=packlow; i < packhigh; i++)
		{
			lb0[u4_0(reader[i])]++;
			lb1[u4_1(reader[i])]++;
			lb2[u4_2(reader[i])]++;
		}
		
		lock.lock();
		for ( size_t i = 0; i < 3*HIST_SIZE; ++i )
			b0[i] += lb0[i];
		lock.unlock();
	}
#else	
	/* create histograms */
	for (size_t i=0; i < n; i++)
	{
		b0[u4_0(reader[i])]++;
		b1[u4_1(reader[i])]++;
		b2[u4_2(reader[i])]++;
	}
#endif

#if 0
	for ( size_t i = 0; i < HIST_SIZE; ++i )
		if ( b2[i] )
			std::cerr << i << " -> " << b2[i] << std::endl;
#endif

	/* compute prefix sums */
	u_int64_t sum0=0,sum1=0,sum2=0,tsum=0;
	for (unsigned int j = 0; j < HIST_SIZE; j++ )
	{
		tsum  = b0[j] + sum0; b0[j] = sum0 - 1; sum0  = tsum;
		tsum  = b1[j] + sum1; b1[j] = sum1 - 1; sum1  = tsum;
		tsum  = b2[j] + sum2; b2[j] = sum2 - 1; sum2  = tsum;
	}
	
	/* sort */
	AutoArray< mask_type > Awriter(n,false);
	mask_type *writer = Awriter.get();

	// correct reader -> writer
	for (size_t i=0; i < n; i++) writer[++ b0[u4_0(reader[i])] ] = reader[i]; std::swap(reader,writer);
	// inverse writer -> reader
	for (size_t i=0; i < n; i++) writer[++ b1[u4_1(reader[i])] ] = reader[i]; std::swap(reader,writer);
	// correct reader -> writer
	for (size_t i=0; i < n; i++) writer[++ b2[u4_2(reader[i])] ] = reader[i]; std::swap(reader,writer);

	// assign
	Areader = Awriter;
}
