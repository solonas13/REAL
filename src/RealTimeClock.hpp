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

#if ! defined(REALTIMECLOCK_HPP)
#define REALTIMECLOCK_HPP

#include "types.hpp"

#if defined(HAVE_WINDOWS_H)
#include <windows.h>
#endif

#if defined(HAVE_SYS_TIME_H)
#include <sys/time.h>
#endif

namespace toollib {
#if defined(__linux) || defined(__FreeBSD__) || defined(__APPLE__)
  typedef u_int64_t rtc_u_int64_t;

  struct RealTimeClock {
    private:
    struct timeval started;
    mutable struct timezone tz;
    
    public:
    RealTimeClock() {
      tz.tz_minuteswest = 0;
      tz.tz_dsttime = 0;
    }
    ~RealTimeClock() throw() {}
    
    bool start() throw() {
      return gettimeofday(&started,&tz) == 0;
    }
    
    //! elapsed time in u-secs
    rtc_u_int64_t getElapsed() const {
      struct timeval stopped;
      gettimeofday(&stopped,&tz);
      struct timeval dif;
      timersub(&stopped,&started,&dif);
      return 
        static_cast<rtc_u_int64_t>(dif.tv_usec) + (
          static_cast<rtc_u_int64_t>(dif.tv_sec) * 
          static_cast<rtc_u_int64_t>(1000000ul)
      );
    }
    double getElapsedSeconds() const {
    	rtc_u_int64_t const t = getElapsed();
    	rtc_u_int64_t const s = t / 1000000;
    	rtc_u_int64_t const u = t % 1000000;
    	double ds = static_cast<double>(s);
    	double du = static_cast<double>(u)/1000000.0;
    	return (ds+du);
    }
  };
#else
  typedef unsigned __int64 rtc_u_int64_t;

  struct RealTimeClock {
    private:
    LARGE_INTEGER freq;
    LARGE_INTEGER started;  

    static __int64 ltd(LARGE_INTEGER const & l) {
	__int64 result = l.HighPart;
	result <<= 32;
	result |= l.LowPart;
	return result;
    } 
    
    public:
    RealTimeClock() {
      QueryPerformanceFrequency(&freq);
    }
    ~RealTimeClock() throw() {}
    
    bool start() throw() {
      if ( ltd(freq) == 0 )
        return false;
      else {
        return QueryPerformanceCounter(&started) > 0;
      }
    }
    
    rtc_u_int64_t getElapsed() const {
      if ( ! ltd(freq) )
        return 0;
      else {
        LARGE_INTEGER now;
        QueryPerformanceCounter(&now);
        return 
          ( static_cast<rtc_u_int64_t>(ltd(now)-ltd(started)) *
          static_cast<rtc_u_int64_t>(1000000ul) )
          / static_cast<rtc_u_int64_t>(ltd(freq));
      }
    }

    double getElapsedSeconds() const {
    	rtc_u_int64_t const t = getElapsed();
    	rtc_u_int64_t const s = t / 1000000;
    	rtc_u_int64_t const u = t % 1000000;
    	double ds = static_cast<double>(s);
    	double du = static_cast<double>(u)/1000000.0;
    	return (ds+du);
    }
  };
#endif
}

#endif
