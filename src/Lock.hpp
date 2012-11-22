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

#if ! defined(LOCK_HPP)
#define LOCK_HPP

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace toollib
{
	struct Lock
	{
#if defined(_OPENMP)
		omp_lock_t lock_obj;
#endif

		Lock()
		{
#if defined(_OPENMP)
			omp_init_lock(&lock_obj);
#endif
		}
		~Lock()
		{
#if defined(_OPENMP)
			omp_destroy_lock(&lock_obj);
#endif
		}
		void lock()
		{
#if defined(_OPENMP)
			omp_set_lock(&lock_obj);
#endif		
		}
		void unlock()
		{
#if defined(_OPENMP)
			omp_unset_lock(&lock_obj);
#endif	
		}
	};
}
#endif
