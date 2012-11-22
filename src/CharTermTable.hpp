#if ! defined(CHARTERMTABLE_HPP)
#define CHARTERMTABLE_HPP

#include "AutoArray.hpp"

struct CharTermTable
{
	private:
	AutoArray<bool> atable;
	bool * table;
	
	public:
	bool operator[](int i) const
	{
		return table[i];
	}
	
	CharTermTable(u_int8_t c)
	: atable(257), table(atable.get()+1)
	{
		for ( unsigned int i = 0; i < 256; ++i )
			table[i] = false;
		table[-1] = true;
		table[c] = true;
	}
};
#endif
