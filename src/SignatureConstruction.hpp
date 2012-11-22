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

#if ! defined(SIGNATURECONSTRUCTION_HPP)
#define SIGNATURECONSTRUCTION_HPP

#include "types.hpp"

#include <string>
#include <cassert>

template<typename signature_type>
struct SignatureConstruction
{
        unsigned int const l;
        unsigned int const nu;
        
        unsigned int const syms_m0;
        unsigned int const syms_m1;
        unsigned int const syms_m2;
        unsigned int const syms_m3;
        unsigned int const bits_m0;
        unsigned int const bits_m1;
        unsigned int const bits_m2;
        unsigned int const bits_m3;
        
        u_int32_t const mask_m0;
        u_int32_t const mask_m1;
        u_int32_t const mask_m2;
        u_int32_t const mask_m3;
        
        SignatureConstruction(unsigned int const rl, unsigned int const rnu)
        : l(rl), nu(rnu), syms_m0(l/nu), syms_m1(syms_m0), syms_m2(syms_m0), syms_m3(l - syms_m0 - syms_m1 - syms_m2),
          bits_m0 ( syms_m0 << 1), bits_m1 ( syms_m1 << 1), bits_m2 ( syms_m2 << 1), bits_m3 ( syms_m3 << 1),
          mask_m0 ( (bits_m0 < 8*sizeof(mask_m0)) ? ((static_cast<u_int32_t>(1) << bits_m0)-1 ) : (~static_cast<u_int32_t>(0))),
          mask_m1 ( (bits_m1 < 8*sizeof(mask_m1)) ? ((static_cast<u_int32_t>(1) << bits_m1)-1 ) : (~static_cast<u_int32_t>(0))),
          mask_m2 ( (bits_m2 < 8*sizeof(mask_m2)) ? ((static_cast<u_int32_t>(1) << bits_m2)-1 ) : (~static_cast<u_int32_t>(0))),
          mask_m3 ( (bits_m3 < 8*sizeof(mask_m3)) ? ((static_cast<u_int32_t>(1) << bits_m3)-1 ) : (~static_cast<u_int32_t>(0)))
        {}
        
        signature_type s0m0(signature_type const s0) const { return s0 >> bits_m1; }
        signature_type s0m1(signature_type const s0) const { return s0 & mask_m1; }

        signature_type s5m2(signature_type const s5) const { return s5 >> bits_m3; }
        signature_type s5m3(signature_type const s5) const { return s5 & mask_m3; }

        signature_type s0(u_int64_t const m0, u_int64_t const m1) const { return (m0 << bits_m1) | m1; }
        signature_type s1(u_int64_t const m0, u_int64_t const m2) const { return (m0 << bits_m2) | m2; }
        signature_type s2(u_int64_t const m0, u_int64_t const m3) const { return (m0 << bits_m3) | m3; }
        signature_type s3(u_int64_t const m1, u_int64_t const m2) const { return (m1 << bits_m2) | m2; }
        signature_type s4(u_int64_t const m1, u_int64_t const m3) const { return (m1 << bits_m3) | m3; }
        signature_type s5(u_int64_t const m2, u_int64_t const m3) const { return (m2 << bits_m3) | m3; }        
        
        unsigned int s0bits() const { return bits_m0 + bits_m1; }
        unsigned int s1bits() const { return bits_m0 + bits_m2; }
        unsigned int s2bits() const { return bits_m0 + bits_m3; }
        unsigned int s3bits() const { return bits_m1 + bits_m2; }
        unsigned int s4bits() const { return bits_m1 + bits_m3; }
        unsigned int s5bits() const { return bits_m2 + bits_m3; }
        
        unsigned int sampleshift(unsigned int const bits, unsigned int const samplebits) const
        {
                return (bits >= samplebits) ? (bits - samplebits) : 0;
        }
        
        unsigned int s0shift(unsigned int const samplebits) const { return sampleshift( s0bits() , samplebits ); }
        unsigned int s1shift(unsigned int const samplebits) const { return sampleshift( s1bits() , samplebits ); }
        unsigned int s2shift(unsigned int const samplebits) const { return sampleshift( s2bits() , samplebits ); }
        unsigned int s3shift(unsigned int const samplebits) const { return sampleshift( s3bits() , samplebits ); }
        unsigned int s4shift(unsigned int const samplebits) const { return sampleshift( s4bits() , samplebits ); }
        unsigned int s5shift(unsigned int const samplebits) const { return sampleshift( s5bits() , samplebits ); }

        void shift(u_int32_t & lm0, u_int32_t & lm1, u_int32_t & lm2, u_int32_t & lm3) const
        {
                lm0 <<= 2; lm0 &= mask_m0; lm0 |= (lm1 >> (bits_m1-2));
                lm1 <<= 2; lm1 &= mask_m1; lm1 |= (lm2 >> (bits_m2-2));
                lm2 <<= 2; lm2 &= mask_m2; lm2 |= (lm3 >> (bits_m3-2));
                lm3 <<= 2; lm3 &= mask_m3;
        }

        static unsigned int invert(unsigned int v)
        {
                switch ( v )
                {
                        case 0: return 3;
                        case 3: return 0;
                        case 1: return 2;
                        case 2: return 1;
                        default: return 0;
                }
        }
        static char invertChar(char v)
        {
                switch ( v )
                {
                        case 'A': return 'T';
                        case 'T': return 'A';
                        case 'C': return 'G';
                        case 'G': return 'C';
                        default: return v;
                }
        }

        void reverse(
        	u_int32_t & lm0, 
        	u_int32_t & lm1, 
        	u_int32_t & lm2, 
        	u_int32_t & lm3
	) const
        {
                // copy
                u_int32_t tm0 = lm0, tm1 = lm1, tm2 = lm2, tm3 = lm3; lm0 = lm1 = lm2 = lm3 = 0;

                for ( unsigned int i = 0; i < syms_m3; ++i ) { shift(lm0,lm1,lm2,lm3); lm3 |= invert(tm3 & 0x3); tm3 >>= 2;  }
                for ( unsigned int i = 0; i < syms_m2; ++i ) { shift(lm0,lm1,lm2,lm3); lm3 |= invert(tm2 & 0x3); tm2 >>= 2;  }
                for ( unsigned int i = 0; i < syms_m1; ++i ) { shift(lm0,lm1,lm2,lm3); lm3 |= invert(tm1 & 0x3); tm1 >>= 2;  }
                for ( unsigned int i = 0; i < syms_m0; ++i ) { shift(lm0,lm1,lm2,lm3); lm3 |= invert(tm0 & 0x3); tm0 >>= 2;  }
        }
        
        static void reverseString(std::string & s)
        {
                unsigned int const n = s.size();
                unsigned int const n2 = n>>1;

                unsigned int i = 0;
                unsigned int j = n-1;
                
                for ( ; i != n2; ++i, --j )
                {
                        char const a = invertChar(s[i]);
                        s[i] = invertChar(s[j]);
                        s[j] = a;
                }
                
                if ( n & 1 )
                        s[n2] = invertChar(s[n2]);
        }

        template<typename iterator>
        inline bool signature(
                iterator w, 
                u_int32_t * vm
                ) const
        {
                assert ( l <= 64 );

                u_int32_t m0 = 0;
                
                for ( unsigned int i = 0; i < syms_m0; ++i )
                        switch( *(w++) )
                        {
                                case 'A': m0 = (m0 << 2) | 0; break;
                                case 'C': m0 = (m0 << 2) | 1; break;
                                case 'G': m0 = (m0 << 2) | 2; break;
                                case 'T': m0 = (m0 << 2) | 3; break;
                                default: return false;
                        }

                u_int32_t m1 = 0;

                for ( unsigned int i = 0; i < syms_m1; ++i )
                        switch( *(w++) )
                        {
                                case 'A': m1 = (m1 << 2) | 0; break;
                                case 'C': m1 = (m1 << 2) | 1; break;
                                case 'G': m1 = (m1 << 2) | 2; break;
                                case 'T': m1 = (m1 << 2) | 3; break;
                                default: return false;
                        }

                u_int32_t m2 = 0;

                for ( unsigned int i = 0; i < syms_m2; ++i )
                        switch( *(w++) )
                        {
                                case 'A': m2 = (m2 << 2) | 0; break;
                                case 'C': m2 = (m2 << 2) | 1; break;
                                case 'G': m2 = (m2 << 2) | 2; break;
                                case 'T': m2 = (m2 << 2) | 3; break;
                                default: return false;
                        }

                u_int32_t m3 = 0;

                for ( unsigned int i = 0; i < syms_m3; ++i )
                        switch( *(w++) )
                        {
                                case 'A': m3 = (m3 << 2) | 0; break;
                                case 'C': m3 = (m3 << 2) | 1; break;
                                case 'G': m3 = (m3 << 2) | 2; break;
                                case 'T': m3 = (m3 << 2) | 3; break;
                                default: return false;
                        }

                vm[0] = m0;
                vm[1] = m1;
                vm[2] = m2;
                vm[3] = m3;
                
                return true;
        }

        template<typename iterator>
        inline bool signatureMapped(
                iterator w, 
                u_int32_t * vm
                ) const
        {
                assert ( l <= 64 );

                u_int32_t m0 = 0;
                
                for ( unsigned int i = 0; i < syms_m0; ++i )
                        switch( *(w++) )
                        {
                                case 0: m0 = (m0 << 2) | 0; break;
                                case 1: m0 = (m0 << 2) | 1; break;
                                case 2: m0 = (m0 << 2) | 2; break;
                                case 3: m0 = (m0 << 2) | 3; break;
                                default: return false;
                        }

                u_int32_t m1 = 0;

                for ( unsigned int i = 0; i < syms_m1; ++i )
                        switch( *(w++) )
                        {
                                case 0: m1 = (m1 << 2) | 0; break;
                                case 1: m1 = (m1 << 2) | 1; break;
                                case 2: m1 = (m1 << 2) | 2; break;
                                case 3: m1 = (m1 << 2) | 3; break;
                                default: return false;
                        }

                u_int32_t m2 = 0;

                for ( unsigned int i = 0; i < syms_m2; ++i )
                        switch( *(w++) )
                        {
                                case 0: m2 = (m2 << 2) | 0; break;
                                case 1: m2 = (m2 << 2) | 1; break;
                                case 2: m2 = (m2 << 2) | 2; break;
                                case 3: m2 = (m2 << 2) | 3; break;
                                default: return false;
                        }

                u_int32_t m3 = 0;

                for ( unsigned int i = 0; i < syms_m3; ++i )
                        switch( *(w++) )
                        {
                                case 0: m3 = (m3 << 2) | 0; break;
                                case 1: m3 = (m3 << 2) | 1; break;
                                case 2: m3 = (m3 << 2) | 2; break;
                                case 3: m3 = (m3 << 2) | 3; break;
                                default: return false;
                        }

                vm[0] = m0;
                vm[1] = m1;
                vm[2] = m2;
                vm[3] = m3;
                
                return true;
        }

        template<typename iterator>
        inline bool reverseSignature(
                iterator w,
                u_int32_t * vm
                ) const
        {
                assert ( l <= 64 );

                u_int32_t m;

                w += syms_m3;
                m = 0;
                for ( unsigned int i = 0; i < syms_m3; ++i )
                        switch( *(--w) )
                        {
                                case 'T': m = (m << 2) | 0; break;
                                case 'G': m = (m << 2) | 1; break;
                                case 'C': m = (m << 2) | 2; break;
                                case 'A': m = (m << 2) | 3; break;
                                default: return false;
                        }
		vm[3] = m;		

                w += (syms_m3+syms_m2);
                m = 0;
                for ( unsigned int i = 0; i < syms_m2; ++i )
                        switch( *(--w) )
                        {
                                case 'T': m = (m << 2) | 0; break;
                                case 'G': m = (m << 2) | 1; break;
                                case 'C': m = (m << 2) | 2; break;
                                case 'A': m = (m << 2) | 3; break;
                                default: return false;
                        }
		vm[2] = m;		

                w += (syms_m2+syms_m1);
                m = 0;
                for ( unsigned int i = 0; i < syms_m1; ++i )
                        switch( *(--w) )
                        {
                                case 'T': m = (m << 2) | 0; break;
                                case 'G': m = (m << 2) | 1; break;
                                case 'C': m = (m << 2) | 2; break;
                                case 'A': m = (m << 2) | 3; break;
                                default: return false;
                        }
		vm[1] = m;		

                w += (syms_m1+syms_m0);
                m = 0;
                for ( unsigned int i = 0; i < syms_m0; ++i )
                        switch( *(--w) )
                        {
                                case 'T': m = (m << 2) | 0; break;
                                case 'G': m = (m << 2) | 1; break;
                                case 'C': m = (m << 2) | 2; break;
                                case 'A': m = (m << 2) | 3; break;
                                default: return false;
                        }
		vm[0] = m;
                
                return true;
        }

        template<typename iterator>
        inline bool reverseMappedSignature(
                iterator w,
                u_int32_t * vm
                ) const
        {
                assert ( l <= 64 );

                u_int32_t m;

                w += syms_m3;
                m = 0;
                for ( unsigned int i = 0; i < syms_m3; ++i )
                        switch( *(--w) )
                        {
                                case 3: m = (m << 2) | 0; break;
                                case 2: m = (m << 2) | 1; break;
                                case 1: m = (m << 2) | 2; break;
                                case 0: m = (m << 2) | 3; break;
                                default: return false;
                        }
		vm[3] = m;		

                w += (syms_m3+syms_m2);
                m = 0;
                for ( unsigned int i = 0; i < syms_m2; ++i )
                        switch( *(--w) )
                        {
                                case 3: m = (m << 2) | 0; break;
                                case 2: m = (m << 2) | 1; break;
                                case 1: m = (m << 2) | 2; break;
                                case 0: m = (m << 2) | 3; break;
                                default: return false;
                        }
		vm[2] = m;		

                w += (syms_m2+syms_m1);
                m = 0;
                for ( unsigned int i = 0; i < syms_m1; ++i )
                        switch( *(--w) )
                        {
                                case 3: m = (m << 2) | 0; break;
                                case 2: m = (m << 2) | 1; break;
                                case 1: m = (m << 2) | 2; break;
                                case 0: m = (m << 2) | 3; break;
                                default: return false;
                        }
		vm[1] = m;		

                w += (syms_m1+syms_m0);
                m = 0;
                for ( unsigned int i = 0; i < syms_m0; ++i )
                        switch( *(--w) )
                        {
                                case 3: m = (m << 2) | 0; break;
                                case 2: m = (m << 2) | 1; break;
                                case 1: m = (m << 2) | 2; break;
                                case 0: m = (m << 2) | 3; break;
                                default: return false;
                        }
		vm[0] = m;
                
                return true;
        }

        std::string maskToString(u_int32_t m0, u_int32_t m1, u_int32_t m2, u_int32_t m3) const
        {
                std::string s(l,' ');
                unsigned int j = 0;
                
                for ( unsigned int i = 0; i < syms_m0; ++i )
                        switch ( 
                                ( m0 >> ( 2*(syms_m0-i-1) ) ) & 0x3
                        )
                        {
                                case 0: s[j++] = 'A'; break; case 1: s[j++] = 'C'; break; case 2: s[j++] = 'G'; break; case 3: s[j++] = 'T'; break;
                        }
                for ( unsigned int i = 0; i < syms_m1; ++i )
                        switch ( 
                                ( m1 >> ( 2*(syms_m1-i-1) ) ) & 0x3
                        )
                        {
                                case 0: s[j++] = 'A'; break; case 1: s[j++] = 'C'; break; case 2: s[j++] = 'G'; break; case 3: s[j++] = 'T'; break;
                        }
                for ( unsigned int i = 0; i < syms_m2; ++i )
                        switch ( 
                                ( m2 >> ( 2*(syms_m2-i-1) ) ) & 0x3
                        )
                        {
                                case 0: s[j++] = 'A'; break; case 1: s[j++] = 'C'; break; case 2: s[j++] = 'G'; break; case 3: s[j++] = 'T'; break;
                        }
                for ( unsigned int i = 0; i < syms_m3; ++i )
                        switch ( 
                                ( m3 >> ( 2*(syms_m3-i-1) ) ) & 0x3
                        )
                        {
                                case 0: s[j++] = 'A'; break; case 1: s[j++] = 'C'; break; case 2: s[j++] = 'G'; break; case 3: s[j++] = 'T'; break;
                        }

                return s;
        }
};
#endif
