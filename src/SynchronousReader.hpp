#if ! defined(SYNCHRONOUSREADER_HPP)
#define SYNCHRONOUSREADER_HPP

#include "AutoArray.hpp"
#include "PatternBlock.hpp"

template<typename reader_type>
struct SynchronousReader
{
        typedef typename reader_type::pattern_type pattern_type;
        
        reader_type & file;
        
        unsigned int const blocksize;
        AutoArray < pattern_type > block;
        PatternBlock< pattern_type > pblock;
        
        SynchronousReader (reader_type & rfile , unsigned int const rblocksize)
        : file(rfile), blocksize(rblocksize), block(blocksize)
        {
                pblock.patterns = block.get();
                pblock.blockid = 0;
                pblock.blocksize = 0;
        }

        PatternBlock<pattern_type> * getBlock()
        {
                pblock.blocksize = file.fillPatternBlock(pblock.patterns,blocksize);
                
                if ( pblock.blocksize )
                        return & pblock;
                else
                        return 0;
        }
        void returnBlock(PatternBlock<pattern_type> *)
        {
        }
        
        unsigned int getFillState() const
        {
                return 1;
        }
};
#endif
