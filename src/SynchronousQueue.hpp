#if ! defined(SYNCHRONOUSQUEUE_HPP)
#define SYNCHRONOUSQUEUE_HPP

#if defined(HAVE_PTHREADS) || defined(__APPLE__)
#include "PosixMutex.hpp"
#include "PosixSemaphore.hpp"
#include <deque>

template<typename value_type>
struct SynchronousQueue
{
        std::deque < value_type > Q;
        PosixMutex lock;
        PosixSemaphore semaphore;
        
        unsigned int getFillState()
        {
                lock.lock();
                unsigned int const fill = Q.size();
                lock.unlock();
                return fill;
        }
        
        void enque(value_type const q)
        {
                lock.lock();
                Q.push_back(q);
                lock.unlock();
                semaphore.post();
        }
        value_type deque()
        {
                semaphore.wait();
                lock.lock();
                value_type const v = Q.front();
                Q.pop_front();
                lock.unlock();
                return v;
        }
};
#endif
#endif
