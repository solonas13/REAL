/*
Copyright (c) 2008, 2009 Andrew I. Schein
Copyright (c) 2010, German Tischler, Solon Pissis

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

template<typename mask_type>
inline u_int64_t u8_0(mask_type const & v) { return ((v.sign)         & 0x7FFULL); }
template<typename mask_type>
inline u_int64_t u8_1(mask_type const & v) { return (((v.sign) >> 11) & 0x7FFULL); }
template<typename mask_type>
inline u_int64_t u8_2(mask_type const & v) { return (((v.sign) >> 22) & 0x7FFULL); }
template<typename mask_type>
inline u_int64_t u8_3(mask_type const & v) { return (((v.sign) >> 33) & 0x7FFULL); }
template<typename mask_type>
inline u_int64_t u8_4(mask_type const & v) { return (((v.sign) >> 44) & 0x7FFULL); }
template<typename mask_type>
inline u_int64_t u8_5(mask_type const & v) { return (((v.sign) >> 55) & 0x7FFULL); }

template<typename mask_type>
inline void u8_sort(AutoArray < mask_type > & Areader, size_t const n) 
{
	static unsigned int const HIST_SIZE = 2048;

	mask_type * reader = Areader.get();
	if (n < HIST_SIZE) { std::sort(reader,reader+n); return; }
	
	/* allocate 6 lists of HIST_SIZE elements */
	AutoArray< u_int64_t > Ab0(6 * HIST_SIZE, false);
	u_int64_t * const b0   = Ab0.get();
	u_int64_t * const b1   = b0 + HIST_SIZE;
	u_int64_t * const b2   = b1 + HIST_SIZE;
	u_int64_t * const b3   = b2 + HIST_SIZE;
	u_int64_t * const b4   = b3 + HIST_SIZE;
	u_int64_t * const b5   = b4 + HIST_SIZE;

	/* erase lists */
	std::fill(b0, b0+6*HIST_SIZE, 0);

	/* create histograms */
	for (size_t i=0; i < n; i++)
	{
		b0[u8_0(reader[i])]++;
		b1[u8_1(reader[i])]++;
		b2[u8_2(reader[i])]++;
		b3[u8_3(reader[i])]++;
		b4[u8_4(reader[i])]++;
		b5[u8_5(reader[i])]++;
	}

	/* compute prefix sums */
	u_int64_t sum0=0,sum1=0,sum2=0,sum3=0,sum4=0,sum5=0,tsum=0;
	for (unsigned int j = 0; j < HIST_SIZE; j++)
	{
		tsum  = b0[j] + sum0; b0[j] = sum0 - 1; sum0  = tsum;
		tsum  = b1[j] + sum1; b1[j] = sum1 - 1; sum1  = tsum;
		tsum  = b2[j] + sum2; b2[j] = sum2 - 1; sum2  = tsum;
		tsum  = b3[j] + sum3; b3[j] = sum3 - 1; sum3  = tsum;
		tsum  = b4[j] + sum4; b4[j] = sum4 - 1; sum4  = tsum;
		tsum  = b5[j] + sum5; b5[j] = sum5 - 1; sum5  = tsum;
	}

	AutoArray< mask_type > Awriter(n,false);
	mask_type *writer = Awriter.get();

	/* sort */
	for (size_t i=0; i < n; i++) writer[++b0[u8_0(reader[i])]] = reader[i]; std::swap(reader,writer);
	for (size_t i=0; i < n; i++) writer[++b1[u8_1(reader[i])]] = reader[i]; std::swap(reader,writer);
	for (size_t i=0; i < n; i++) writer[++b2[u8_2(reader[i])]] = reader[i]; std::swap(reader,writer);
	for (size_t i=0; i < n; i++) writer[++b3[u8_3(reader[i])]] = reader[i]; std::swap(reader,writer);
	for (size_t i=0; i < n; i++) writer[++b4[u8_4(reader[i])]] = reader[i]; std::swap(reader,writer);
	for (size_t i=0; i < n; i++) writer[++b5[u8_5(reader[i])]] = reader[i];
}
