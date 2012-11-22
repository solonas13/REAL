#if !defined(ASYNCHRONOUSREADER_HPP)
#define ASYNCHRONOUSREADER_HPP

#include "AutoArray.hpp"
#include "TerminatableSynchronousQueue.hpp"

#if defined(HAVE_PTHREADS) || defined(__APPLE__)
template<typename reader_type>
struct AsynchronousStreamReaderData
{
        typedef typename reader_type::pattern_type pattern_type;
        typedef typename reader_type::block_type data_type;

        reader_type & file;
        
        unsigned int const blocksize;
        unsigned int const numblocks;

        AutoArray < pattern_type > blocks;
        AutoArray < data_type > blockp;

        AsynchronousStreamReaderData(
                reader_type & rfile , 
                unsigned int const rblocksize, 
                unsigned int const rnumblocks 
        )
        : file(rfile), blocksize(rblocksize), numblocks(rnumblocks), blocks(blocksize*numblocks), blockp(numblocks)
        {
                std::cerr << "Using cache for " << blocksize*numblocks << " patterns." << std::endl;
        
                for ( unsigned int i = 0; i < numblocks; ++i )
                {
                        blockp[i].patterns = blocks.get() + i * blocksize;
                        blockp[i].blocksize = 0;
                        blockp[i].blockid = i;
                }
        }
        
        unsigned int getNumEntities()
        {
                return numblocks;
        }
        
        data_type * getData(unsigned int const id)
        {
                return &blockp[id];
        }
        unsigned int dataToId(data_type const * data)
        {
               return data->blockid; 
        }
        bool generateData(unsigned int const blockid)
        {
		// std::cerr << "Filling pattern block " << blockid << std::endl;
                blockp[blockid].blocksize = file.fillPatternBlock(blockp[blockid].patterns,blocksize);        
		// std::cerr << "Block size " << blockp[blockid].blocksize << std::endl;
                return blockp[blockid].blocksize != 0;
        }
};

template<typename reader_type>
struct AsynchronousFastReaderData
{
        typedef typename reader_type::pattern_type pattern_type;
        typedef typename reader_type::block_type data_type;

        reader_type & file;
        
        unsigned int const numblocks;
        unsigned int const blocksize;

        AutoArray < data_type > blockp;

        AsynchronousFastReaderData(reader_type & rfile , unsigned int const rblocksize, unsigned int const rnumblocks)
        : file(rfile), numblocks(rnumblocks), blocksize(rblocksize), blockp(numblocks)
        {
                for ( unsigned int i = 0; i < numblocks; ++i )
                        blockp[i].blockid = i;
        }
        
        unsigned int getNumEntities()
        {
                return numblocks;
        }
        
        data_type * getData(unsigned int const id)
        {
                return &blockp[id];
        }
        unsigned int dataToId(data_type const * data)
        {
               return data->blockid; 
        }
        bool generateData(unsigned int const blockid)
        {
                blockp[blockid].blocksize = file.fillPatternBlock(blockp[blockid],blocksize);
                return blockp[blockid].blocksize != 0;
        }
};

template<typename reader_type>
struct AsynchronousFastIdData
{
        typedef typename reader_type::idfile_type file_type;
        typedef typename reader_type::idblock_type data_type;

        file_type & file;
        
        unsigned int const numblocks;
        unsigned int const blocksize;

        AutoArray < data_type > blockp;

        AsynchronousFastIdData(file_type & rfile , unsigned int const rblocksize, unsigned int const rnumblocks)
        : file(rfile), numblocks(rnumblocks), blocksize(rblocksize), blockp(numblocks)
        {
                for ( unsigned int i = 0; i < numblocks; ++i )
                        blockp[i].blockid = i;
        }
        
        unsigned int getNumEntities()
        {
                return numblocks;
        }
        
        data_type * getData(unsigned int const id)
        {
                return &blockp[id];
        }
        unsigned int dataToId(data_type const * data)
        {
               return data->blockid; 
        }
        bool generateData(unsigned int const blockid)
        {
                blockp[blockid].blocksize = file.fillPatternBlock(blockp[blockid],blocksize);
                return blockp[blockid].blocksize != 0;
        }
};

template<typename reader_type>
struct AsynchronousIdData
{
        typedef reader_type file_type;
        typedef typename reader_type::idblock_type data_type;

        file_type & file;
        
        unsigned int const numblocks;
        unsigned int const blocksize;

        AutoArray < data_type > blockp;

        AsynchronousIdData(file_type & rfile , unsigned int const rblocksize, unsigned int const rnumblocks)
        : file(rfile), numblocks(rnumblocks), blocksize(rblocksize), blockp(numblocks)
        {
                for ( unsigned int i = 0; i < numblocks; ++i )
                        blockp[i].blockid = i;
        }
        
        unsigned int getNumEntities()
        {
                return numblocks;
        }
        
        data_type * getData(unsigned int const id)
        {
                return &blockp[id];
        }
        unsigned int dataToId(data_type const * data)
        {
               return data->blockid; 
        }
        bool generateData(unsigned int const blockid)
        {
                blockp[blockid].blocksize = file.fillIdBlock(blockp[blockid],blocksize);
                return blockp[blockid].blocksize != 0;
        }
};

template<typename data_type>
struct AsynchronousStreamReader
{
        data_type & data;        
        
        SynchronousQueue<unsigned int> unfilled;
        TerminatableSynchronousQueue<unsigned int> filled;
        pthread_t thread;
        
        typename data_type::data_type * getBlock()
        {
                unsigned int blockid;
                try
                {
                        blockid = filled.deque();
			// std::cerr << "dequed block id " << blockid << std::endl;
                }
                catch(...)
                {
			//std::cerr << ex.what() << std::endl;
                        return 0;
                }

                return data.getData(blockid);
        }
        
        void returnBlock(typename data_type::data_type * block)
        {
                unfilled.enque(data.dataToId(block));
        }
        
        unsigned int getFillState()
        {
                return filled.getFillState();
        }

        void * readerThread()
        {
                while ( !filled.isTerminated() )
                {
                        unsigned int const blockid = unfilled.deque();

                        bool const generated = data.generateData(blockid);
                        
                        if ( ! generated )
                                filled.terminate();
                        else
			{
				// std::cerr << "Queued block " << blockid << std::endl;
                                filled.enque(blockid);
			}
                }

                return 0;
        }
        
        static void * dispatcher ( void * object )
        {
                AsynchronousStreamReader * reader = reinterpret_cast < AsynchronousStreamReader * > (object);
                pthread_exit ( reader->readerThread() );
                return 0;
        }

        AsynchronousStreamReader ( data_type & rdata ) : data(rdata)
        {
                for ( unsigned int i = 0; i < data.getNumEntities(); ++i )
                        unfilled.enque(i);
                
                if ( pthread_create ( & thread, 0, dispatcher, this ) )
                {
                        throw std::runtime_error("Cannot create reader thread.");
                }
        }
        ~AsynchronousStreamReader()
        {
                void * p = 0;
                pthread_join ( thread, & p );
        }
};
#endif

#endif
