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

#if ! defined(UNIQUEMATCHINFO_HPP)
#define UNIQUEMATCHINFO_HPP

#include <limits>

struct UniqueMatchInfoBase
{
        typedef u_int64_t data_type;
        static unsigned int const databits = sizeof(data_type)*8;
        
        static unsigned int const statebits = 3;
        static unsigned int const errorbits = 4;
        static unsigned int const fileidbits = 6;
        static unsigned int const fragmentbits = 16;
        static unsigned int const posbits = databits - statebits - errorbits - fileidbits - fragmentbits; // 64 - (2+4+6+16) = 35

        static unsigned int const posshift = 0;
        static unsigned int const fileidshift = posshift + posbits;
        static unsigned int const errorshift = fileidshift + fileidbits;
        static unsigned int const fragmentshift = errorshift + errorbits;
        static unsigned int const stateshift = fragmentshift + fragmentbits;
        
        static data_type const fragmentmask  = ((static_cast<data_type>(1) << fragmentbits) - 1) << fragmentshift;
        static data_type const statemask  = ((static_cast<data_type>(1) << statebits) - 1) << stateshift;
        static data_type const errormask  = ((static_cast<data_type>(1) << errorbits) - 1) << errorshift;
        static data_type const fileidmask = ((static_cast<data_type>(1) << fileidbits) - 1) << fileidshift;
        static data_type const posmask = ((static_cast<data_type>(1) << posbits) - 1) << posshift;
        
        static data_type const nonstatemask =  ~statemask;
        static data_type const nonposmask = ~posmask;
        static data_type const nonerrormask = ~errormask;
        static data_type const nonfileidmask = ~fileidmask;
        static data_type const nonfragmentmask = ~fragmentmask;

        static data_type const straightmask = static_cast<u_int64_t>(1) << stateshift;
        static data_type const reversemask = static_cast<u_int64_t>(2) << stateshift;
        static data_type const gappedmask = static_cast<u_int64_t>(3) << stateshift;
        static data_type const nonuniquemask = static_cast<u_int64_t>(4) << stateshift;
        
        static unsigned int getMaxErrors()
        {
                return (1ul << errorbits) - 1;
        }
        
        static unsigned int getMaxFragmentsPerFile()
        {
                return (1ul << fragmentbits);
        }
        
        // data
        u_int64_t data;
                
        enum MatchState
        {
                NoMatch = 0,
                Straight = 1,
                Reverse = 2,
                Gapped = 3,
                NonUnique = 4
        };
        
        MatchState getState() const
        {
                switch ( data >> stateshift )
                {
                        case 0:
                                return NoMatch;
                        case 1:
                                return Straight;
                        case 2:
                                return Reverse;
                        case 3:
                                return Gapped;
                        case 4:
                        default:
                                return NonUnique;
                }
        }
        
        void setState(MatchState state)
        {
                // mask out state
                data &= nonstatemask;

                switch ( state )
                {
                        case NoMatch:
                                break;
                        case Straight:
                                data |= straightmask;
                                break;
                        case Reverse:
                                data |= reversemask;
                                break;
                        case Gapped:
                                data |= gappedmask;
                                break;
                        case NonUnique:
                                data |= nonuniquemask;
                                break;
                }
        }
        
        void setPosition(u_int64_t pos)
        {
                data &= nonposmask;
                data |= pos;
        }
        
        u_int64_t getPosition() const
        {
                return data & posmask;
        }

        unsigned int getErrors() const
        {
                return (data & errormask) >> errorshift;
        }
        
        void setErrors(unsigned int errors)
        {
                data &= nonerrormask;
                data |= static_cast<u_int64_t>(errors) << errorshift;
        }

        unsigned int getFileid() const
        {
                return (data & fileidmask) >> fileidshift;
        }
        
        void setFileid(unsigned int fileid)
        {
                data &= nonfileidmask;
                data |= static_cast<u_int64_t>(fileid) << fileidshift;
        }

        unsigned int getFragment() const
        {
                return (data & fragmentmask) >> fragmentshift;
        }
        
        void setFragment(unsigned int fragment)
        {
                assert ( (fragment >> fragmentbits) == 0 );
                data &= nonfragmentmask;
                data |= static_cast<u_int64_t>(fragment) << fragmentshift;
        }

        UniqueMatchInfoBase() : data(0) {}
};

template<bool scores>
struct UniqueMatchInfo {};

template<>
struct UniqueMatchInfo<false> : public UniqueMatchInfoBase
{
        UniqueMatchInfo() : UniqueMatchInfoBase() {}

        void setScore(float const)
        {
                
        }
        float getScore() const
        {
                return 0.0f;
        }
};

template<>
struct UniqueMatchInfo<true> : public UniqueMatchInfoBase
{
        UniqueMatchInfo() : UniqueMatchInfoBase(), score( - std::numeric_limits<float>::max() ) {}

        float score;

        void setScore(float const rscore)
        {
                score = rscore;
        }
        float getScore() const
        {
                return score;
        }
};
#endif
