#if ! defined(PATTERNBLOCK_HPP)
#define PATTERNBLOCK_HPP

template<typename pattern_type>
struct PatternBlock
{
        pattern_type * patterns;
        unsigned int blockid;
        unsigned int blocksize;
        
        pattern_type const & getPattern(u_int64_t i) const
        {
                patterns[i].computeMapped();
                return patterns[i];
        }
};
#endif
