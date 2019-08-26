/* An implementation of Cantor's algorithm for multiplication of
   polynomials over GF(2).
   
   Copyright 2008,2009,2010,2012 Emmanuel Thomé.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at
   your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with CADO-NFS; see the file COPYING.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/

#ifndef FFT_ADAPTER_HPP_
#define FFT_ADAPTER_HPP_


#define CAT(X,Y) X ## _ ## Y

/* _setup and _dft, _ift integer arguments correspond to lengths in
 * number of bits of polynomials.
 *
 * _alloc(), _free(), _get(), take number of complete polynomials.
 *
 * the size() result is the number of type-t objects which represent
 * _one_ polynomial. alloc(1) is thus the same as malloc(size())
 */

#define DEFINE_FFT_ADAPTER(visible)					\
struct visible {							\
	CAT(visible,info_t) o;						\
	typedef CAT(visible,ptr) ptr;						\
	typedef CAT(visible,srcptr) srcptr;					\
									\
	visible() {}							\
	visible(size_t n1, size_t n2) { CAT(visible,init)(o, n1, n2); }	\
	~visible() { CAT(visible,clear)(o); }                 		\
									\
private: /* Make sure we forbid copies */                               \
        visible(visible const& other GF2X_MAYBE_UNUSED) {}                   \
public:                                                                 \
        /* The extra <acc> argument gives the number of times the */    \
        /* fft is going to be reused. It might be replaced by some */   \
        /* other mechanism at some point. */                            \
									\
        /* n1 and n2 give the length of the two operands. The */        \
        /* product is taken as having length n1+n2-1 */                 \
									\
	visible(size_t n1, size_t n2, int acc GF2X_MAYBE_UNUSED)		\
	{								\
		CAT(visible,init)(o, n1, n2, acc);				\
	}								\
        visible(size_t n1, size_t n2, visible const& other)             \
        {                                                               \
            CAT(visible,init_similar)(o, n1, n2, other.o);                    \
        }                                                               \
	inline ptr alloc(size_t n) const { return CAT(visible,alloc)(o, n); }	\
	inline void free(ptr x, size_t n) const { CAT(visible,free)(o,x,n); }	\
	inline void zero(ptr x, size_t n) const { CAT(visible,zero)(o,x,n); }	\
	inline ptr get(ptr x, size_t n) const { return CAT(visible,get)(o,x,n); }\
	inline srcptr get(srcptr x, size_t n) const { return CAT(visible,get_const)(o,x,n); }\
	inline void dft(ptr x, const unsigned long * F, size_t n) const {	\
		return CAT(visible,dft)(o,x,F,n);				\
	}								\
        /* this destroys the input ! */                                 \
	inline void ift(unsigned long * F, size_t n, ptr x) const {	\
		return CAT(visible,ift)(o,F,n,x);				\
	}								\
	inline void addcompose(ptr y, srcptr x1, srcptr x2) const {	\
		return CAT(visible,addcompose)(o,y,x1,x2);		\
	}								\
	inline void addcompose_n(ptr y, srcptr * x1, srcptr * x2, size_t n) const {	\
		return CAT(visible,addcompose_n)(o,y,x1,x2,n);		\
	}								\
	inline void compose(ptr y, srcptr x1, srcptr x2) const {	\
		return CAT(visible,compose)(o,y,x1,x2);			\
	}								\
	inline void add(ptr y, srcptr x1, srcptr x2) const {		\
		return CAT(visible,add)(o,y,x1,x2);				\
	}								\
	inline void cpy(ptr y, srcptr x) const {		        \
		return CAT(visible,cpy)(o,y,x);				\
	}								\
	inline size_t size() const {	        	                \
		return CAT(visible,size)(o);	        			\
	}								\
        static inline const char * name() { return #visible; }          \
};


#endif	/* FFT_ADAPTER_HPP_ */
