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

#if ! defined(MASK_HPP)
#define MASK_HPP

template<typename _signature_type, typename _ptr_type>
struct BaseMask
{
        typedef _signature_type signature_type;
        typedef _ptr_type ptr_type;

        signature_type sign;
        ptr_type ptr;
        
        BaseMask() {}
        BaseMask(signature_type rsign) : sign(rsign) {}
        BaseMask(signature_type rsign, ptr_type rptr)
        : sign(rsign), ptr(rptr) {}

        template<typename list_ptr_type>
        inline u_int32_t getPos(list_ptr_type complementlist) const
        {
        	return complementlist[ptr].getPos(0);
        }
};

template<typename signature_type, typename ptr_type>
struct Mask : public BaseMask<signature_type,ptr_type>
{
	private:
        u_int32_t pos;
        
        public:
        Mask() {}
        Mask(signature_type rsign) : BaseMask<signature_type,ptr_type>(rsign) {}
        Mask(signature_type rsign, ptr_type rptr, size_t rpos)
        : BaseMask<signature_type,ptr_type>(rsign,rptr), pos(rpos) {}
        
        template<typename list_ptr_type>
        inline u_int32_t getPos(list_ptr_type) const
        {
        	return pos;
        }
        inline void setPos(u_int32_t const rpos)
        {
        	pos = rpos;
        }
};

template<typename signature_type, typename ptr_type>
inline bool operator<(BaseMask<signature_type, ptr_type> const & A, BaseMask<signature_type, ptr_type> const & B)
{
        return A.sign < B.sign;
}
#endif
