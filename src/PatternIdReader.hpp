#if ! defined(PATTERNIDREADER_HPP)
#define PATTERNIDREADER_HPP

#include <memory>
#include <string>
#include "FastFileDecoderBase.hpp"
#include "Lock.hpp"

template<typename reader_type>
struct PatternIdReader
{
        /**
         * pattern file
         **/
        std::auto_ptr < reader_type > patfile;
        std::auto_ptr < typename reader_type::stream_data_type > ad;
        std::auto_ptr < typename reader_type::stream_reader_type > ar;

        /**
         * id file
         **/
        std::auto_ptr < typename reader_type::idfile_type > idfile;
        std::auto_ptr < typename reader_type::stream_id_type > id;
        std::auto_ptr < typename reader_type::stream_idreader_type > idr;
        
        /**
         * lock
         **/
        toollib::Lock lock;

        PatternIdReader(std::string const & patternfilename, int qualityOffset, unsigned int numthreads)
        : patfile ( new reader_type (patternfilename,qualityOffset) ),
          ad ( new typename reader_type::stream_data_type ( *(patfile.get()), FastFileDecoderBase::default_blocksize,std::max(4u,3*numthreads) ) ),
          ar ( new typename reader_type::stream_reader_type ( *(ad.get()) ) ),
          idfile ( new typename reader_type::idfile_type(patternfilename)  ),
          id( new typename reader_type::stream_id_type( *(idfile.get()), FastFileDecoderBase::default_blocksize,std::max(4u,3*numthreads) ) ),
          idr( new typename reader_type::stream_idreader_type ( *(id.get()) ) )
        {
        
        }
        
        bool getBlock(std::pair < typename reader_type::block_type *, typename reader_type::idblock_type * > & block)
        {
                lock.lock();
                block.first = ar->getBlock();
                block.second = idr->getBlock();
                lock.unlock();
                
                return block.first;
        }
        
        void returnBlock(std::pair < typename reader_type::block_type *, typename reader_type::idblock_type * > block)
        {
                lock.lock();
                ar->returnBlock(block.first);
                idr->returnBlock(block.second);
                lock.unlock();
        }
};
#endif
