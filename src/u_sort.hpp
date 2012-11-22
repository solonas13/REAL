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

#if ! defined(U_SORT_HPP)
#define U_SORT_HPP

#include <cstring> 
#include <cstdlib>
#include <algorithm>

#include "AutoArray.hpp"
#include "Mask.hpp"

template<typename mask1_type, typename mask2_type>
inline void updatePrev(
	mask1_type const * a, 
	mask2_type * b, 
	size_t const n) 
{
#if defined(_OPENMP)
#pragma omp parallel for
#endif
	for ( size_t i = 0; i < n; ++i )
		b [ a[i].ptr ] . ptr = i;
}

template<typename mask_type>
struct MaskProjector
{
	inline typename mask_type::signature_type operator()(mask_type const & mask) const
	{
		return mask.sign;
	}
};

#if defined(OLD_SORT)
#include "u4_sort.hpp"
#include "u8_sort.hpp"
#else
#include "ParallelRadixSort.hpp"
#endif
#include <iostream>

template<typename signature_type, typename ptr_type>
struct MaskSort
{
	static void sort(
		AutoArray< Mask<signature_type,ptr_type> > & a, 
		AutoArray< BaseMask<signature_type,ptr_type> > & b, 
		size_t const n,
		int inumthreads = -1)
	{
		std::sort(a.get(),a.get()+n);
		updatePrev(a.get(),b.get(),n);
                
		std::sort(b.get(),b.get()+n);
		updatePrev(b.get(),a.get(),n);
	}
};

template<typename ptr_type>
struct MaskSort<u_int64_t, ptr_type>
{
	static void sort(
		AutoArray < Mask<u_int64_t,ptr_type> > & a,
		AutoArray < BaseMask<u_int64_t,ptr_type> > & b,
		size_t const n,
		int inumthreads = -1)
	{
#if defined(OLD_SORT)
		u8_sort(a,n);
#else
		typedef MaskProjector< Mask<u_int64_t,ptr_type> > a_projector_type;
		a_projector_type a_projector;
		RadixSort64< Mask<u_int64_t,ptr_type>, a_projector_type , false >::radixSort(a,n,a_projector,inumthreads);
#endif
		updatePrev(a.get(),b.get(),n);
                
#if defined(OLD_SORT)
		u8_sort(b,n);
#else
		typedef MaskProjector< BaseMask<u_int64_t,ptr_type> > b_projector_type;
		b_projector_type b_projector;
		RadixSort64< BaseMask<u_int64_t,ptr_type> , b_projector_type, false >::radixSort(b,n,b_projector,inumthreads);
#endif
		updatePrev(b.get(),a.get(),n);
	}
};

template<typename ptr_type>
struct MaskSort<u_int32_t, ptr_type>
{
	static void sort(
		AutoArray< Mask<u_int32_t,ptr_type> > & a, 
		AutoArray< BaseMask<u_int32_t,ptr_type> > & b, 
		size_t const n,
		int inumthreads = -1)
	{
#if defined(OLD_SORT)
		u4_sort(a,n);
#else
		typedef MaskProjector< Mask<u_int32_t,ptr_type> > a_projector_type;
		a_projector_type a_projector;
		RadixSort32< Mask<u_int32_t,ptr_type>, a_projector_type , false >::radixSort(a,static_cast<u_int64_t>(n),a_projector,inumthreads);
#endif
		updatePrev(a.get(),b.get(),n);
                
#if defined(OLD_SORT)
		u4_sort(b,n);
#else
		typedef MaskProjector< BaseMask<u_int32_t,ptr_type> > b_projector_type;
		b_projector_type b_projector;
		RadixSort32< BaseMask<u_int32_t,ptr_type>, b_projector_type , false >::radixSort(b,static_cast<u_int64_t>(n),b_projector,inumthreads);
#endif
		updatePrev(b.get(),a.get(),n);
	}
};
#endif
