#if ! defined(POSIXSEMAPHORE_HPP)
#define POSIXSEMAPHORE_HPP

#if defined(HAVE_PTHREADS) || defined(__APPLE__)

#include <ctime>
#include <cerrno>
#include <semaphore.h>
#include <sys/time.h>
#include <string>
#include <sstream>

struct PosixSemaphore
{
	#if defined(__APPLE__)
	std::string semname;
	#endif
        sem_t semaphore;
	sem_t * psemaphore;

        PosixSemaphore()
        {
		#if defined(__APPLE__)
		std::ostringstream semnamestr;
		semnamestr << "real_" << getpid() << "_" << reinterpret_cast<u_int64_t>(this);
		semname = semnamestr.str();
		#endif

                memset ( & semaphore, 0, sizeof(sem_t) );        

		#if defined(__APPLE__)
		psemaphore = sem_open(semname.c_str(), O_CREAT | O_EXCL, 0777, 0);
		if ( psemaphore == SEM_FAILED )
		{
			std::cerr << "sem_open failed: " << strerror(errno) << std::endl;
			throw std::runtime_error("Failed to sam_open.");
		}
		#else
                if ( sem_init ( & semaphore, 0, 0 ) )
		{
			std::cerr << "sem_init failed: " << strerror(errno) << std::endl;
                        throw std::runtime_error("Failed to initialise semaphore.");
		}
		psemaphore = &semaphore;
		#endif
        }
        ~PosixSemaphore()
        {
		#if defined(__APPLE__)
		sem_close(psemaphore);
		sem_unlink(semname.c_str());
		#else
                sem_destroy(psemaphore);        
		#endif
        }

        void post()
        {
                sem_post (psemaphore);
        }

        void wait()
        {
                sem_wait (psemaphore);
        }
        
        bool timedWait()
        {
		#if ! defined(__APPLE__)
                struct timeval tv;
                struct timezone tz;
                struct timespec waittime;

                gettimeofday(&tv,&tz);
                waittime.tv_sec = tv.tv_sec + 1;
                waittime.tv_nsec = static_cast<u_int64_t>(tv.tv_usec)*1000;

                int const v = sem_timedwait(psemaphore,&waittime);
                
                if ( v < 0 )
                {
                        if ( errno == ETIMEDOUT )
                                return false;
                        else
                                throw std::runtime_error("timedWait failed on semaphore.");
                }
                else
                {
                        return true;
                }
		#else
		int const v = sem_trywait(psemaphore);

		if ( v < 0 )
		{
			if ( errno == EAGAIN )
				return false;
			else
			{
				std::cerr << "sem_trywait failed with error " << strerror(errno) << " on semaphore " << semname << std::endl;
				throw std::runtime_error("sem_trywait failed on semaphore.");
			}
		}   
		else
		{
			return true;
		}
		#endif
        }
};
#endif
#endif
