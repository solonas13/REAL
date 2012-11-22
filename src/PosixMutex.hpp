#if ! defined(POSIXMUTEX_HPP)
#define POSIXMUTEX_HPP

#if defined(HAVE_PTHREADS) || defined(__APPLE__)
#include <pthread.h>

struct PosixMutex
{
        pthread_mutex_t mutex;
        
        PosixMutex()
        {
                pthread_mutex_init(&mutex,0);
        }
        ~PosixMutex()
        {
                pthread_mutex_destroy(&mutex);
        }
        
        void lock()
        {
                pthread_mutex_lock ( &mutex );
        }
        void unlock()
        {
                pthread_mutex_unlock ( &mutex );
        }
};
#endif
#endif
