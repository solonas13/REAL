#if ! defined(TERMINATABLESYNCHRONOUSQUEUE_HPP)
#define TERMINATABLESYNCHRONOUSQUEUE_HPP

#include "SynchronousQueue.hpp"

#if defined(HAVE_PTHREADS) || defined(__APPLE__)
template<typename value_type>
struct TerminatableSynchronousQueue : public SynchronousQueue<value_type>
{
        typedef SynchronousQueue<value_type> parent_type;

        PosixMutex terminatelock;
        volatile bool terminated;
        
        TerminatableSynchronousQueue()
        {
                terminated = false;
        }
        
        void terminate()
        {
                terminatelock.lock();
                terminated = true;
                terminatelock.unlock();
        }
        bool isTerminated()
        {
                terminatelock.lock();
                bool const isTerm = terminated;
                terminatelock.unlock();
                return isTerm;
        }
        
        value_type deque()
        {
                while ( !parent_type::semaphore.timedWait() )
                {
                        terminatelock.lock();
                        bool const isterminated = terminated;
                        terminatelock.unlock();
                        
                        if ( isterminated )
                                throw std::runtime_error("Queue is terminated");
                }
                
                parent_type::lock.lock();
                value_type const v = parent_type::Q.front();
                parent_type::Q.pop_front();
                parent_type::lock.unlock();
                return v;
        }
};
#endif
#endif
