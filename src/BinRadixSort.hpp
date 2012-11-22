#if ! defined(BINRADIXSORT_HPP)
#define BINRADIXSORT_HPP

#include <sys/types.h>
#include <iterator>
#include <vector>
#include <algorithm>

#if defined(_OPENMP)
#include <omp.h>
#endif

template<typename iterator, typename projector_type>
struct BinRadixSort
{
	static u_int64_t binsort(
		iterator a, 
		iterator e, 
		typename std::iterator_traits<iterator>::value_type const mask,
		projector_type const & projector
	)
	{
		iterator const ta = a;
		iterator const te = e;
		
		if ( a != e )
		{
			e--;
		
			while ( a < e )
			{
				while ( (a < e) && (!(projector(*a) & mask)) )
					++a;
				while ( (a < e) && ((projector(*e) & mask)) )
					--e;
				if ( a < e )
					std::swap(*(a++),*(e--));
			}
		}
		
		while ( (a != te) && (!(projector(*a) & mask)) )
			a++;
		
		return a-ta;
	}

	static void binradixsort(
		iterator a, 
		iterator e, 
		typename std::iterator_traits<iterator>::value_type const startmask, 
		projector_type const & projector,
		u_int64_t parts = 0)
	{
#if defined(_OPENMP)
		if ( ! parts )
			parts = omp_get_max_threads()*4;
#endif
		std::vector < u_int64_t > sortvec;
		
		u_int64_t numparts = 1;
		sortvec.push_back(0);
		sortvec.push_back(e-a);
		
		unsigned int loop = 0;
		
		while ( (numparts < parts) && ((startmask >> loop)!=0) )
		{
#if 0
			for ( unsigned int i = 0; i < sortvec.size(); ++i )
				std::cerr << sortvec[i] << ";";
			std::cerr << std::endl;
#endif
		
			std::vector < u_int64_t > newsortvec(2*sortvec.size()-1);
			newsortvec.back() = sortvec.back();

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,1)
#endif		
			for ( unsigned int j = 0; j < numparts; ++j )
			{
				u_int64_t const left = sortvec[j];
				u_int64_t const right = sortvec[j+1];
				u_int64_t const c0 = binsort(a + left, a + right,  startmask >> loop, projector);
				
				newsortvec[2*j] = left;
				newsortvec[2*j+1] = left + c0;
			}

			loop++;
			numparts *= 2;	
			sortvec.swap(newsortvec);
		}

		if ( (startmask >> loop) != 0 )
		{
#if 0
			for ( unsigned int i = 0; i < sortvec.size(); ++i )
				std::cerr << sortvec[i] << ";";
			std::cerr << std::endl;
#endif

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,1)
#endif	
			for ( u_int64_t j = 0; j < numparts; ++j )
				std::sort( a + sortvec[j], a + sortvec[j+1], projector );
		}
		else
		{
#if 0
			std::cerr << "Skipping final stage, already sorted." << std::endl;
#endif
		}
	}
};
#endif
