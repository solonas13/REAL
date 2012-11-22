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

#if ! defined(MAPTEXTFILE_HPP)
#define MAPTEXTFILE_HPP

#include <string>
#include <fstream>

#include "types.hpp"

#include <iostream>
#include <stdexcept>
#include <deque>
#include <sstream>

#include "Mask.hpp"
#include "stringFunctions.hpp"
#include "SignatureConstruction.hpp"
#include "AutoTextArray.hpp"
#include "ERank222B.hpp"

template<typename signature_type, bool sse4>
struct MapTextFile : public SignatureConstruction<signature_type>
{
        typedef SignatureConstruction<signature_type> base_type;

        AutoTextArray<sse4> const & ATA;

        signature_type m0;
        signature_type m1;
        signature_type m2;
        signature_type m3;
        
        bool readFirst;
        size_t syms_read;
        size_t m_pos;
        size_t r_pos;
        
        bool peeked;

        MapTextFile(
                AutoTextArray<sse4> const & rATA, 
                unsigned int const rl, 
                unsigned int const rnu
        )
        : SignatureConstruction<signature_type>(rl,rnu), ATA(rATA),
          m0(0), m1(0), m2(0), m3(0),
          readFirst(false),
          syms_read(0),
          m_pos(0), r_pos(0), peeked(false) {}
        
        void shift()
        {
                m0 <<= 2; m0 &= base_type::mask_m0; m0 |= (m1 >> (base_type::bits_m1-2));
                m1 <<= 2; m1 &= base_type::mask_m1; m1 |= (m2 >> (base_type::bits_m2-2));
                m2 <<= 2; m2 &= base_type::mask_m2; m2 |= (m3 >> (base_type::bits_m3-2));
                m3 <<= 2; m3 &= base_type::mask_m3;
        }
             
        inline int readSymbol()
        {
                int const sym = (r_pos >= ATA.getN()) ? -1 : ATA[r_pos++];
                
                switch ( sym )
                {
                        case 0:
                        case 1:
                        case 2:
                        case 3:
                        case 4:
                                syms_read += 1;
                                break;
                }
                
                return sym;
        }
        
        inline int readNextSymbol()
        {
                int sym;
                
                while ( (sym = readSymbol()) >= 0 )
                {
                        switch ( sym )
                        {
                                case 0:
                                case 1:
                                case 2:
                                case 3:
                                case 4:
                                        return sym;
                                        break;
                        }
                }
                
                return sym;
        }
        
        bool readNextSignature()
        {
                if ( ! readFirst )
                {
                        return readFullSignature();
                }
                else
                {
                        int sym = readSymbol();
                
                        while ( true )
                        {
                                // end of stream
                                if ( sym < 0 )
                                {
                                        return false;
                                }
                                else
                                {
                                        switch ( sym )
                                        {
                                                case 0: case 1: case 2: case 3: shift(); m3 |= sym; m_pos = syms_read - base_type::l; return true;
                                                default: return readFullSignature();
                                        }
                                }
                        }
                }
        }
        
        bool readFullSignature()
        {
                bool read = false;

                do
                {
                        read = true;
                
                        unsigned int toread = base_type::l;
                        
                        for ( unsigned int i = 0; read && (i < toread); ++i )
                        {
                                int const sym = readSymbol();
                                
                                // end of stream
                                if ( sym < 0 )
                                {
                                        return false;
                                }
                                else
                                {
                                        switch ( sym )
                                        {
                                                case 0: case 1: case 2: case 3: shift(); m3 |= sym; break;
                                                default: read = false; break;
                                        }
                                }
                        }
                        
                } while ( ! read );
                
                readFirst = true;
                m_pos = syms_read - base_type::l;
                
                return true;
        }
        
        template<typename ptr_type>
        unsigned int readLists(
                Mask<signature_type, ptr_type> * const list_0,
                Mask<signature_type, ptr_type> * const list_1,
                Mask<signature_type, ptr_type> * const list_2,
                BaseMask<signature_type, ptr_type> * const list_3,
                BaseMask<signature_type, ptr_type> * const list_4,
                BaseMask<signature_type, ptr_type> * const list_5,
                unsigned int const n,
                bool & havenext
        )
        {
                unsigned int i = 0;
                
                if ( i < n && peeked )
                {
                        list_0[ i ].sign = (m0 << base_type::bits_m1) | (m1); list_0[ i ].ptr = i; list_0[ i ].setPos(m_pos);
                        list_1[ i ].sign = (m0 << base_type::bits_m2) | (m2); list_1[ i ].ptr = i; list_1[ i ].setPos(m_pos);
                        list_2[ i ].sign = (m0 << base_type::bits_m3) | (m3); list_2[ i ].ptr = i; list_2[ i ].setPos(m_pos);
                        list_3[ i ].sign = (m1 << base_type::bits_m2) | (m2); list_3[ i ].ptr = i;
                        list_4[ i ].sign = (m1 << base_type::bits_m3) | (m3); list_4[ i ].ptr = i;
                        list_5[ i ].sign = (m2 << base_type::bits_m3) | (m3); list_5[ i ].ptr = i;

                        i++;

                        peeked = false;
                }

                while ( i < n && readNextSignature() )
                {
                        list_0[ i ].sign = (m0 << base_type::bits_m1) | (m1); list_0[ i ].ptr = i; list_0[ i ].setPos(m_pos);
                        list_1[ i ].sign = (m0 << base_type::bits_m2) | (m2); list_1[ i ].ptr = i; list_1[ i ].setPos(m_pos);
                        list_2[ i ].sign = (m0 << base_type::bits_m3) | (m3); list_2[ i ].ptr = i; list_2[ i ].setPos(m_pos);
                        list_3[ i ].sign = (m1 << base_type::bits_m2) | (m2); list_3[ i ].ptr = i;
                        list_4[ i ].sign = (m1 << base_type::bits_m3) | (m3); list_4[ i ].ptr = i;
                        list_5[ i ].sign = (m2 << base_type::bits_m3) | (m3); list_5[ i ].ptr = i;

                        i++;
                }
                                
                havenext = false;
                                
                if ( readNextSignature() )
                {
                        peeked = true;
                        havenext = true;
                }
                        
                return i;
        }
};
#endif
