/*
 * GLAM - GLSL Linear Algebra Math Library
 * 
 * Copyright (c) 2014, Gregor Riepl <onitake@gmail.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 *     Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *     Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef GLAM_VECTOR_SSE_H
#define GLAM_VECTOR_SSE_H

#include <glam/config.h>

#ifdef GLAM_VECTOR_SSE

#include <x86intrin.h>
#define GLAM_VECTOR_TEMPLATE_DEFER
#include <glam/vector.h>
#undef GLAM_VECTOR_TEMPLATE_DEFER

namespace glam {

// Macros

// Tell the template implementation that we override the specialized templates
#undef GLAM_VECTOR_TEMPLATE_OVERLOADED

// Permutation generator macros
#define GLAM_PERM_MAKE2_2(type, length, name, member0, member1) \
	Vector<type, length, Permutatator<type, length, member0, member1> > name;
#define GLAM_PERM_MAKE3_3(type, length, name, member0, member1, member2) \
	Vector<type, length, Permutatator<type, length, member0, member1, member2> > name;
#define GLAM_PERM_MAKE4_4(type, length, name, member0, member1, member2, member3) \
	Vector<type, length, Permutatator<type, length, member0, member1, member2, member3> > name;
#define GLAM_PERM_MAKE4_3(type, name0, name1, name2, name3, length, prefix, member0, member1, member2) \
	GLAM_PERM_MAKE4_4(type, length, prefix##name0, member0, member1, member2, 0) \
	GLAM_PERM_MAKE4_4(type, length, prefix##name1, member0, member1, member2, 1) \
	GLAM_PERM_MAKE4_4(type, length, prefix##name2, member0, member1, member2, 2) \
	GLAM_PERM_MAKE4_4(type, length, prefix##name3, member0, member1, member2, 3) \
	GLAM_PERM_MAKE3_3(type, length, prefix, member0, member1, member2)
#define GLAM_PERM_MAKE4_2(type, name0, name1, name2, name3, length, prefix, member0, member1) \
	GLAM_PERM_MAKE4_3(type, name0, name1, name2, name3, length, prefix##name0, member0, member1, 0) \
	GLAM_PERM_MAKE4_3(type, name0, name1, name2, name3, length, prefix##name1, member0, member1, 1) \
	GLAM_PERM_MAKE4_3(type, name0, name1, name2, name3, length, prefix##name2, member0, member1, 2) \
	GLAM_PERM_MAKE4_3(type, name0, name1, name2, name3, length, prefix##name3, member0, member1, 3) \
	GLAM_PERM_MAKE2_2(type, length, prefix, member0, member1)
#define GLAM_PERM_MAKE4_1(type, name0, name1, name2, name3, length, prefix, member0) \
	GLAM_PERM_MAKE4_2(type, name0, name1, name2, name3, length, prefix##name0, member0, 0) \
	GLAM_PERM_MAKE4_2(type, name0, name1, name2, name3, length, prefix##name1, member0, 1) \
	GLAM_PERM_MAKE4_2(type, name0, name1, name2, name3, length, prefix##name2, member0, 2) \
	GLAM_PERM_MAKE4_2(type, name0, name1, name2, name3, length, prefix##name3, member0, 3)
#define GLAM_PERM_MAKE4_0(type, name0, name1, name2, name3) \
	GLAM_PERM_MAKE4_1(type, name0, name1, name2, name3, 4, name0, 0) \
	GLAM_PERM_MAKE4_1(type, name0, name1, name2, name3, 4, name1, 1) \
	GLAM_PERM_MAKE4_1(type, name0, name1, name2, name3, 4, name2, 2) \
	GLAM_PERM_MAKE4_1(type, name0, name1, name2, name3, 4, name3, 3) \
	struct { type name0; type name1; type name2; type name3; }
#define GLAM_PERM_MAKE3_2(type, name0, name1, name2, length, prefix, member0, member1) \
	GLAM_PERM_MAKE3_3(type, length, prefix##name0, member0, member1, 0) \
	GLAM_PERM_MAKE3_3(type, length, prefix##name1, member0, member1, 1) \
	GLAM_PERM_MAKE3_3(type, length, prefix##name2, member0, member1, 2) \
	GLAM_PERM_MAKE2_2(type, length, prefix, member0, member1)
#define GLAM_PERM_MAKE3_1(type, name0, name1, name2, length, prefix, member0) \
	GLAM_PERM_MAKE3_2(type, name0, name1, name2, length, prefix##name0, member0, 0) \
	GLAM_PERM_MAKE3_2(type, name0, name1, name2, length, prefix##name1, member0, 1) \
	GLAM_PERM_MAKE3_2(type, name0, name1, name2, length, prefix##name2, member0, 2)
#define GLAM_PERM_MAKE3_0(type, name0, name1, name2) \
	GLAM_PERM_MAKE3_1(type, name0, name1, name2, 3, name0, 0) \
	GLAM_PERM_MAKE3_1(type, name0, name1, name2, 3, name1, 1) \
	GLAM_PERM_MAKE3_1(type, name0, name1, name2, 3, name2, 2) \
	struct { type name0; type name1; type name2; }
#define GLAM_PERM_MAKE2_1(type, name0, name1, length, prefix, member0) \
	GLAM_PERM_MAKE2_2(type, length, prefix##name0, member0, 0) \
	GLAM_PERM_MAKE2_2(type, length, prefix##name1, member0, 1)
#define GLAM_PERM_MAKE2_0(type, name0, name1) \
	GLAM_PERM_MAKE2_1(type, name0, name1, 2, name0, 0) \
	GLAM_PERM_MAKE2_1(type, name0, name1, 2, name1, 1) \
	struct { type name0; type name1; }


// Specializations

template <>
class Vector<float, 2, Pass<float, 2> > {
public:
	union {
		float _v[2];
		__m128 _s;
		// Produces permutations: x, y, xx, xy, yx, yy
		GLAM_PERM_MAKE2_0(float, x, y);
		// Produces permutations: r, g, rr, rg, gr, gg
		GLAM_PERM_MAKE2_0(float, r, g);
		// Produces permutations: s, t, ss, st, ts tt
		GLAM_PERM_MAKE2_0(float, s, t);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<float, 2, Pass<float, 2> > &operator =(const Vector<float, SizeV, PermutationV> &v);
	
	inline float &operator [](size_t index);
	inline const float &operator [](size_t index) const;
	
	// Direct access to internal pointer
	inline const float *internal() const;
};

template <>
class Vector<float, 3, Pass<float, 3> > {
public:
	union {
		float _v[3];
		__m128 _s;
		// Produces permutations: x, y, z, xx, xy, xz, yx, yy, yz, zx, zy, zz, xxx, xxy, xxz, xyx, ...
		GLAM_PERM_MAKE3_0(float, x, y, z);
		// Produces permutations: r, g, b, rr, rg, rb, gr, gg, gb, br, bg, bb, rrr, rrg, rrb, rgr, ...
		GLAM_PERM_MAKE3_0(float, r, g, b);
		// Produces permutations: s, t, p, ss, st, sp, ts, tt, tp, ps, pt, pp, sss, sst, ssp, sts, ...
		GLAM_PERM_MAKE3_0(float, s, t, p);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<float, 3, Pass<float, 3> > &operator =(const Vector<float, SizeV, PermutationV> &v);
	
	inline float &operator [](size_t index);
	inline const float &operator [](size_t index) const;
	
	// Direct access to internal pointer
	inline const float *internal() const;
};

template <>
class Vector<float, 4, Pass<float, 4> > {
public:
	union {
		float _v[4];
		__m128 _s;
		// Produces permutations: x, y, z, w, xx, xy, xz, xw, yx, yy, yz, yw, zx, zy, zz, zw, wx, wy, wz, ww, xxx, xxy, xxz, xyx, ..., xxxx, xxxy, ...
		GLAM_PERM_MAKE4_0(float, x, y, z, w);
		// Produces permutations: r, g, b, a, rr, rg, rb, ra, gr, gg, gb, ga, br, bg, bb, ba, ar, ag, ab, aa, rrr, rrg, rrb, rgr, ..., rrrr, rrrg, ...
		GLAM_PERM_MAKE4_0(float, r, g, b, a);
		// Produces permutations: s, t, p, q, ss, st, sp, sq, ts, tt, tp, tq, ps, pt, pp, pq, qs, qt, qp, qq, sss, sst, ssp, sts, ..., ssss, ssst, ...
		GLAM_PERM_MAKE4_0(float, s, t, p, q);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<float, 4, Pass<float, 4> > &operator =(const Vector<float, SizeV, PermutationV> &v);
	
	inline float &operator [](size_t index);
	inline const float &operator [](size_t index) const;
	
	// Direct access to internal pointer
	inline const float *internal() const;
};

template <>
class Vector<int, 2, Pass<int, 2> > {
public:
	union {
		int _v[2];
		__m64 _s;
		// Produces permutations: x, y, xx, xy, yx, yy
		GLAM_PERM_MAKE2_0(int, x, y);
		// Produces permutations: r, g, rr, rg, gr, gg
		GLAM_PERM_MAKE2_0(int, r, g);
		// Produces permutations: s, t, ss, st, ts tt
		GLAM_PERM_MAKE2_0(int, s, t);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<int, 2, Pass<int, 2> > &operator =(const Vector<int, SizeV, PermutationV> &v);
	
	inline int &operator [](size_t index);
	inline const int &operator [](size_t index) const;
	
	// Direct access to internal pointer
	inline const int *internal() const;
};

template <>
class Vector<int, 3, Pass<int, 3> > {
public:
	union {
		int _v[3];
		__m128i _s;
		// Produces permutations: x, y, z, xx, xy, xz, yx, yy, yz, zx, zy, zz, xxx, xxy, xxz, xyx, ...
		GLAM_PERM_MAKE3_0(int, x, y, z);
		// Produces permutations: r, g, b, rr, rg, rb, gr, gg, gb, br, bg, bb, rrr, rrg, rrb, rgr, ...
		GLAM_PERM_MAKE3_0(int, r, g, b);
		// Produces permutations: s, t, p, ss, st, sp, ts, tt, tp, ps, pt, pp, sss, sst, ssp, sts, ...
		GLAM_PERM_MAKE3_0(int, s, t, p);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<int, 3, Pass<int, 3> > &operator =(const Vector<int, SizeV, PermutationV> &v);
	
	inline int &operator [](size_t index);
	inline const int &operator [](size_t index) const;
	
	// Direct access to internal pointer
	inline const int *internal() const;
};

template <>
class Vector<int, 4, Pass<int, 4> > {
public:
	union {
		int _v[4];
		__m128i _s;
		// Produces permutations: x, y, z, w, xx, xy, xz, xw, yx, yy, yz, yw, zx, zy, zz, zw, wx, wy, wz, ww, xxx, xxy, xxz, xyx, ..., xxxx, xxxy, ...
		GLAM_PERM_MAKE4_0(int, x, y, z, w);
		// Produces permutations: r, g, b, a, rr, rg, rb, ra, gr, gg, gb, ga, br, bg, bb, ba, ar, ag, ab, aa, rrr, rrg, rrb, rgr, ..., rrrr, rrrg, ...
		GLAM_PERM_MAKE4_0(int, r, g, b, a);
		// Produces permutations: s, t, p, q, ss, st, sp, sq, ts, tt, tp, tq, ps, pt, pp, pq, qs, qt, qp, qq, sss, sst, ssp, sts, ..., ssss, ssst, ...
		GLAM_PERM_MAKE4_0(int, s, t, p, q);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<int, 4, Pass<int, 4> > &operator =(const Vector<int, SizeV, PermutationV> &v);
	
	inline int &operator [](size_t index);
	inline const int &operator [](size_t index) const;
	
	// Direct access to internal pointer
	inline const int *internal() const;
};

template <>
class Vector<unsigned int, 2, Pass<unsigned int, 2> > {
public:
	union {
		unsigned int _v[2];
		__m64 _s;
		// Produces permutations: x, y, xx, xy, yx, yy
		GLAM_PERM_MAKE2_0(unsigned int, x, y);
		// Produces permutations: r, g, rr, rg, gr, gg
		GLAM_PERM_MAKE2_0(unsigned int, r, g);
		// Produces permutations: s, t, ss, st, ts tt
		GLAM_PERM_MAKE2_0(unsigned int, s, t);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<unsigned int, 2, Pass<unsigned int, 2> > &operator =(const Vector<unsigned int, SizeV, PermutationV> &v);
	
	inline unsigned int &operator [](size_t index);
	inline const unsigned int &operator [](size_t index) const;
	
	// Direct access to internal pointer
	inline const unsigned int *internal() const;
};

template <>
class Vector<unsigned int, 3, Pass<unsigned int, 3> > {
public:
	union {
		unsigned int _v[3];
		__m128i _s;
		// Produces permutations: x, y, z, xx, xy, xz, yx, yy, yz, zx, zy, zz, xxx, xxy, xxz, xyx, ...
		GLAM_PERM_MAKE3_0(unsigned int, x, y, z);
		// Produces permutations: r, g, b, rr, rg, rb, gr, gg, gb, br, bg, bb, rrr, rrg, rrb, rgr, ...
		GLAM_PERM_MAKE3_0(unsigned int, r, g, b);
		// Produces permutations: s, t, p, ss, st, sp, ts, tt, tp, ps, pt, pp, sss, sst, ssp, sts, ...
		GLAM_PERM_MAKE3_0(unsigned int, s, t, p);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<unsigned int, 3, Pass<unsigned int, 3> > &operator =(const Vector<unsigned int, SizeV, PermutationV> &v);
	
	inline unsigned int &operator [](size_t index);
	inline const unsigned int &operator [](size_t index) const;
	
	// Direct access to internal pointer
	inline const unsigned int *internal() const;
};

template <>
class Vector<unsigned int, 4, Pass<unsigned int, 4> > {
public:
	union {
		unsigned int _v[4];
		__m128i _s;
		// Produces permutations: x, y, z, w, xx, xy, xz, xw, yx, yy, yz, yw, zx, zy, zz, zw, wx, wy, wz, ww, xxx, xxy, xxz, xyx, ..., xxxx, xxxy, ...
		GLAM_PERM_MAKE4_0(unsigned int, x, y, z, w);
		// Produces permutations: r, g, b, a, rr, rg, rb, ra, gr, gg, gb, ga, br, bg, bb, ba, ar, ag, ab, aa, rrr, rrg, rrb, rgr, ..., rrrr, rrrg, ...
		GLAM_PERM_MAKE4_0(unsigned int, r, g, b, a);
		// Produces permutations: s, t, p, q, ss, st, sp, sq, ts, tt, tp, tq, ps, pt, pp, pq, qs, qt, qp, qq, sss, sst, ssp, sts, ..., ssss, ssst, ...
		GLAM_PERM_MAKE4_0(unsigned int, s, t, p, q);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<unsigned int, 4, Pass<unsigned int, 4> > &operator =(const Vector<unsigned int, SizeV, PermutationV> &v);
	
	inline unsigned int &operator [](size_t index);
	inline const unsigned int &operator [](size_t index) const;
	
	// Direct access to internal pointer
	inline const unsigned int *internal() const;
};

// Undefine generator macros
#undef GLAM_PERM_MAKE2_2
#undef GLAM_PERM_MAKE3_3
#undef GLAM_PERM_MAKE4_4
#undef GLAM_PERM_MAKE4_3
#undef GLAM_PERM_MAKE4_2
#undef GLAM_PERM_MAKE4_1
#undef GLAM_PERM_MAKE4_0
#undef GLAM_PERM_MAKE3_2
#undef GLAM_PERM_MAKE3_1
#undef GLAM_PERM_MAKE3_0
#undef GLAM_PERM_MAKE2_1
#undef GLAM_PERM_MAKE2_0


// Implementation

template <typename... Args>
inline Vector<float, 2, Pass<float, 2> >::Vector(Args... args) {
	initialize<0, float, 2, Pass<float, 2> >(*this, args...);
}

template <>
inline Vector<float, 2, Pass<float, 2> > &Vector<float, 2, Pass<float, 2> >::operator =(const Vector<float, 2, Pass<float, 2> > &v) {
	this->_s = v._s;
	return *this;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 2, Pass<float, 2> > &Vector<float, 2, Pass<float, 2> >::operator =(const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	for (size_t c = 0; c < 2; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

inline float &Vector<float, 2, Pass<float, 2> >::operator [](size_t index) {
	Pass<float, 2> perm;
	return perm(_v, index);
}

inline const float &Vector<float, 2, Pass<float, 2> >::operator [](size_t index) const {
	Pass<float, 2> perm;
	return perm(_v, index);
}

inline const float *Vector<float, 2, Pass<float, 2> >::internal() const {
	return _v;
}

inline Vector<float, 2> &operator ++(Vector<float, 2> &self) {
	self._s = _mm_add_ps(self._s, _mm_set1_ps(1.0f));
	return self;
}

inline Vector<float, 2> &operator --(Vector<float, 2> &self) {
	self._s = _mm_sub_ps(self._s, _mm_set1_ps(1));
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 2> &operator +=(Vector<float, 2> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<float, 2> packed(v);
	self._s = _mm_add_ps(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 2> &operator -=(Vector<float, 2> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<float, 2> packed(v);
	self._s = _mm_sub_ps(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 2> &operator *=(Vector<float, 2> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<float, 2> packed(v);
	self._s = _mm_mul_ps(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 2> &operator /=(Vector<float, 2> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<float, 2> packed(v);
	self._s = _mm_div_ps(self._s, packed._s);
	return self;
}

inline Vector<float, 2> operator +(const Vector<float, 2> &u) {
	return u;
}

inline Vector<float, 2> operator -(const Vector<float, 2> &u) {
	Vector<float, 2> v;
	v._s = _mm_sub_ps(_mm_set1_ps(0.0f), u._s);
	return v;
}

template <typename... Args>
inline Vector<float, 3, Pass<float, 3> >::Vector(Args... args) {
	initialize<0, float, 3, Pass<float, 3> >(*this, args...);
}

template <>
inline Vector<float, 3, Pass<float, 3> > &Vector<float, 3, Pass<float, 3> >::operator =(const Vector<float, 3, Pass<float, 3> > &v) {
	this->_s = v._s;
	return *this;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 3, Pass<float, 3> > &Vector<float, 3, Pass<float, 3> >::operator =(const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	for (size_t c = 0; c < 3; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

inline float &Vector<float, 3, Pass<float, 3> >::operator [](size_t index) {
	Pass<float, 3> perm;
	return perm(_v, index);
}

inline const float &Vector<float, 3, Pass<float, 3> >::operator [](size_t index) const {
	Pass<float, 3> perm;
	return perm(_v, index);
}

inline const float *Vector<float, 3, Pass<float, 3> >::internal() const {
	return _v;
}

inline Vector<float, 3> &operator ++(Vector<float, 3> &self) {
	self._s = _mm_add_ps(self._s, _mm_set1_ps(1.0f));
	return self;
}

inline Vector<float, 3> &operator --(Vector<float, 3> &self) {
	self._s = _mm_sub_ps(self._s, _mm_set1_ps(1));
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 3> &operator +=(Vector<float, 3> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<float, 3> packed(v);
	self._s = _mm_add_ps(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 3> &operator -=(Vector<float, 3> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<float, 3> packed(v);
	self._s = _mm_sub_ps(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 3> &operator *=(Vector<float, 3> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<float, 3> packed(v);
	self._s = _mm_mul_ps(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 3> &operator /=(Vector<float, 3> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<float, 3> packed(v);
	self._s = _mm_div_ps(self._s, packed._s);
	return self;
}

inline Vector<float, 3> operator +(const Vector<float, 3> &u) {
	return u;
}

inline Vector<float, 3> operator -(const Vector<float, 3> &u) {
	Vector<float, 3> v;
	v._s = _mm_sub_ps(_mm_set1_ps(0.0f), u._s);
	return v;
}

template <typename... Args>
inline Vector<float, 4, Pass<float, 4> >::Vector(Args... args) {
	initialize<0, float, 4, Pass<float, 4> >(*this, args...);
}

template <>
inline Vector<float, 4, Pass<float, 4> > &Vector<float, 4, Pass<float, 4> >::operator =(const Vector<float, 4, Pass<float, 4> > &v) {
	this->_s = v._s;
	return *this;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 4, Pass<float, 4> > &Vector<float, 4, Pass<float, 4> >::operator =(const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	for (size_t c = 0; c < 4; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

inline float &Vector<float, 4, Pass<float, 4> >::operator [](size_t index) {
	Pass<float, 4> perm;
	return perm(_v, index);
}

inline const float &Vector<float, 4, Pass<float, 4> >::operator [](size_t index) const {
	Pass<float, 4> perm;
	return perm(_v, index);
}

inline const float *Vector<float, 4, Pass<float, 4> >::internal() const {
	return _v;
}

inline Vector<float, 4> &operator ++(Vector<float, 4> &self) {
	self._s = _mm_add_ps(self._s, _mm_set1_ps(1.0f));
	return self;
}

inline Vector<float, 4> &operator --(Vector<float, 4> &self) {
	self._s = _mm_sub_ps(self._s, _mm_set1_ps(1));
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 4> &operator +=(Vector<float, 4> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<float, 4> packed(v);
	self._s = _mm_add_ps(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 4> &operator -=(Vector<float, 4> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<float, 4> packed(v);
	self._s = _mm_sub_ps(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 4> &operator *=(Vector<float, 4> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<float, 4> packed(v);
	self._s = _mm_mul_ps(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<float, 4> &operator /=(Vector<float, 4> &self, const Vector<float, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<float, 4> packed(v);
	self._s = _mm_div_ps(self._s, packed._s);
	return self;
}

inline Vector<float, 4> operator +(const Vector<float, 4> &u) {
	return u;
}

inline Vector<float, 4> operator -(const Vector<float, 4> &u) {
	Vector<float, 4> v;
	v._s = _mm_sub_ps(_mm_set1_ps(0.0f), u._s);
	return v;
}

template <typename... Args>
inline Vector<int, 2, Pass<int, 2> >::Vector(Args... args) {
	initialize<0, int, 2, Pass<int, 2> >(*this, args...);
}

template <>
inline Vector<int, 2, Pass<int, 2> > &Vector<int, 2, Pass<int, 2> >::operator =(const Vector<int, 2, Pass<int, 2> > &v) {
	this->_s = v._s;
	return *this;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 2, Pass<int, 2> > &Vector<int, 2, Pass<int, 2> >::operator =(const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	for (size_t c = 0; c < 2; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

inline int &Vector<int, 2, Pass<int, 2> >::operator [](size_t index) {
	Pass<int, 2> perm;
	return perm(_v, index);
}

inline const int &Vector<int, 2, Pass<int, 2> >::operator [](size_t index) const {
	Pass<int, 2> perm;
	return perm(_v, index);
}

inline const int *Vector<int, 2, Pass<int, 2> >::internal() const {
	return _v;
}

inline Vector<int, 2> &operator ++(Vector<int, 2> &self) {
	self._s = _m_paddd(self._s, _mm_set1_pi32(1));
	return self;
}

inline Vector<int, 2> &operator --(Vector<int, 2> &self) {
	self._s = _m_psubd(self._s, _mm_set1_pi32(1));
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 2> &operator +=(Vector<int, 2> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<int, 2> packed(v);
	self._s = _m_paddd(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 2> &operator -=(Vector<int, 2> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<int, 2> packed(v);
	self._s = _m_psubd(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 2> &operator &=(Vector<int, 2> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<int, 2> packed(v);
	self._s = _m_pand(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 2> &operator |=(Vector<int, 2> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<int, 2> packed(v);
	self._s = _m_por(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 2> &operator ^=(Vector<int, 2> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<int, 2> packed(v);
	self._s = _m_pxor(self._s, packed._s);
	return self;
}

inline Vector<int, 2> &operator <<=(Vector<int, 2> &self, size_t s) {
	self._s = _m_pslldi(self._s, (int) s);
	return self;
}

inline Vector<int, 2> &operator >>=(Vector<int, 2> &self, size_t s) {
	self._s = _m_psradi(self._s, (int) s);
	return self;
}

inline Vector<int, 2> operator +(const Vector<int, 2> &u) {
	return u;
}

inline Vector<int, 2> operator -(const Vector<int, 2> &u) {
	Vector<int, 2> v;
	v._s = _m_psubd(_mm_set1_pi8(0x00), u._s);
	return v;
}

inline Vector<int, 2> operator ~(const Vector<int, 2> &u) {
	Vector<int, 2> v;
	v._s = _m_pxor(u._s, _mm_set1_pi8(0xff));
	return v;
}

template <typename... Args>
inline Vector<int, 3, Pass<int, 3> >::Vector(Args... args) {
	initialize<0, int, 3, Pass<int, 3> >(*this, args...);
}

template <>
inline Vector<int, 3, Pass<int, 3> > &Vector<int, 3, Pass<int, 3> >::operator =(const Vector<int, 3, Pass<int, 3> > &v) {
	this->_s = v._s;
	return *this;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 3, Pass<int, 3> > &Vector<int, 3, Pass<int, 3> >::operator =(const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	for (size_t c = 0; c < 3; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

inline int &Vector<int, 3, Pass<int, 3> >::operator [](size_t index) {
	Pass<int, 3> perm;
	return perm(_v, index);
}

inline const int &Vector<int, 3, Pass<int, 3> >::operator [](size_t index) const {
	Pass<int, 3> perm;
	return perm(_v, index);
}

inline const int *Vector<int, 3, Pass<int, 3> >::internal() const {
	return _v;
}

inline Vector<int, 3> &operator ++(Vector<int, 3> &self) {
	self._s = _mm_add_epi32(self._s, _mm_set1_epi32(1));
	return self;
}

inline Vector<int, 3> &operator --(Vector<int, 3> &self) {
	self._s = _mm_sub_epi32(self._s, _mm_set1_epi32(1));
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 3> &operator +=(Vector<int, 3> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<int, 3> packed(v);
	self._s = _mm_add_epi32(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 3> &operator -=(Vector<int, 3> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<int, 3> packed(v);
	self._s = _mm_sub_epi32(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 3> &operator &=(Vector<int, 3> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<int, 3> packed(v);
	self._s = _mm_and_si128(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 3> &operator |=(Vector<int, 3> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<int, 3> packed(v);
	self._s = _mm_or_si128(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 3> &operator ^=(Vector<int, 3> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<int, 3> packed(v);
	self._s = _mm_xor_si128(self._s, packed._s);
	return self;
}

inline Vector<int, 3> &operator <<=(Vector<int, 3> &self, size_t s) {
	self._s = _mm_slli_epi32(self._s, (int) s);
	return self;
}

inline Vector<int, 3> &operator >>=(Vector<int, 3> &self, size_t s) {
	self._s = _mm_srai_epi32(self._s, (int) s);
	return self;
}

inline Vector<int, 3> operator +(const Vector<int, 3> &u) {
	return u;
}

inline Vector<int, 3> operator -(const Vector<int, 3> &u) {
	Vector<int, 3> v;
	v._s = _mm_sub_epi32(_mm_set1_epi8(0x00), u._s);
	return v;
}

inline Vector<int, 3> operator ~(const Vector<int, 3> &u) {
	Vector<int, 3> v;
	v._s = _mm_xor_si128(u._s, _mm_set1_epi8(0xff));
	return v;
}

template <typename... Args>
inline Vector<int, 4, Pass<int, 4> >::Vector(Args... args) {
	initialize<0, int, 4, Pass<int, 4> >(*this, args...);
}

template <>
inline Vector<int, 4, Pass<int, 4> > &Vector<int, 4, Pass<int, 4> >::operator =(const Vector<int, 4, Pass<int, 4> > &v) {
	this->_s = v._s;
	return *this;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 4, Pass<int, 4> > &Vector<int, 4, Pass<int, 4> >::operator =(const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	for (size_t c = 0; c < 4; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

inline int &Vector<int, 4, Pass<int, 4> >::operator [](size_t index) {
	Pass<int, 4> perm;
	return perm(_v, index);
}

inline const int &Vector<int, 4, Pass<int, 4> >::operator [](size_t index) const {
	Pass<int, 4> perm;
	return perm(_v, index);
}

inline const int *Vector<int, 4, Pass<int, 4> >::internal() const {
	return _v;
}

inline Vector<int, 4> &operator ++(Vector<int, 4> &self) {
	self._s = _mm_add_epi32(self._s, _mm_set1_epi32(1));
	return self;
}

inline Vector<int, 4> &operator --(Vector<int, 4> &self) {
	self._s = _mm_sub_epi32(self._s, _mm_set1_epi32(1));
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 4> &operator +=(Vector<int, 4> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<int, 4> packed(v);
	self._s = _mm_add_epi32(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 4> &operator -=(Vector<int, 4> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<int, 4> packed(v);
	self._s = _mm_sub_epi32(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 4> &operator &=(Vector<int, 4> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<int, 4> packed(v);
	self._s = _mm_and_si128(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 4> &operator |=(Vector<int, 4> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<int, 4> packed(v);
	self._s = _mm_or_si128(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<int, 4> &operator ^=(Vector<int, 4> &self, const Vector<int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<int, 4> packed(v);
	self._s = _mm_xor_si128(self._s, packed._s);
	return self;
}

inline Vector<int, 4> &operator <<=(Vector<int, 4> &self, size_t s) {
	self._s = _mm_slli_epi32(self._s, (int) s);
	return self;
}

inline Vector<int, 4> &operator >>=(Vector<int, 4> &self, size_t s) {
	self._s = _mm_srai_epi32(self._s, (int) s);
	return self;
}

inline Vector<int, 4> operator +(const Vector<int, 4> &u) {
	return u;
}

inline Vector<int, 4> operator -(const Vector<int, 4> &u) {
	Vector<int, 4> v;
	v._s = _mm_sub_epi32(_mm_set1_epi8(0x00), u._s);
	return v;
}

inline Vector<int, 4> operator ~(const Vector<int, 4> &u) {
	Vector<int, 4> v;
	v._s = _mm_xor_si128(u._s, _mm_set1_epi8(0xff));
	return v;
}

template <typename... Args>
inline Vector<unsigned int, 2, Pass<unsigned int, 2> >::Vector(Args... args) {
	initialize<0, unsigned int, 2, Pass<unsigned int, 2> >(*this, args...);
}

template <>
inline Vector<unsigned int, 2, Pass<unsigned int, 2> > &Vector<unsigned int, 2, Pass<unsigned int, 2> >::operator =(const Vector<unsigned int, 2, Pass<unsigned int, 2> > &v) {
	this->_s = v._s;
	return *this;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 2, Pass<unsigned int, 2> > &Vector<unsigned int, 2, Pass<unsigned int, 2> >::operator =(const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	for (size_t c = 0; c < 2; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

inline unsigned int &Vector<unsigned int, 2, Pass<unsigned int, 2> >::operator [](size_t index) {
	Pass<unsigned int, 2> perm;
	return perm(_v, index);
}

inline const unsigned int &Vector<unsigned int, 2, Pass<unsigned int, 2> >::operator [](size_t index) const {
	Pass<unsigned int, 2> perm;
	return perm(_v, index);
}

inline const unsigned int *Vector<unsigned int, 2, Pass<unsigned int, 2> >::internal() const {
	return _v;
}

inline Vector<unsigned int, 2> &operator ++(Vector<unsigned int, 2> &self) {
	self._s = _m_paddd(self._s, _mm_set1_pi32(1));
	return self;
}

inline Vector<unsigned int, 2> &operator --(Vector<unsigned int, 2> &self) {
	self._s = _m_psubd(self._s, _mm_set1_pi32(1));
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 2> &operator +=(Vector<unsigned int, 2> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<unsigned int, 2> packed(v);
	self._s = _m_paddd(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 2> &operator -=(Vector<unsigned int, 2> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<unsigned int, 2> packed(v);
	self._s = _m_psubd(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 2> &operator &=(Vector<unsigned int, 2> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<unsigned int, 2> packed(v);
	self._s = _m_pand(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 2> &operator |=(Vector<unsigned int, 2> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<unsigned int, 2> packed(v);
	self._s = _m_por(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 2> &operator ^=(Vector<unsigned int, 2> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	Vector<unsigned int, 2> packed(v);
	self._s = _m_pxor(self._s, packed._s);
	return self;
}

inline Vector<unsigned int, 2> &operator <<=(Vector<unsigned int, 2> &self, size_t s) {
	self._s = _m_pslldi(self._s, (unsigned int) s);
	return self;
}

inline Vector<unsigned int, 2> &operator >>=(Vector<unsigned int, 2> &self, size_t s) {
	self._s = _m_psrldi(self._s, (unsigned int) s);
	return self;
}

inline Vector<unsigned int, 2> operator +(const Vector<unsigned int, 2> &u) {
	return u;
}

inline Vector<unsigned int, 2> operator -(const Vector<unsigned int, 2> &u) {
	Vector<unsigned int, 2> v;
	v._s = _m_psubd(_mm_set1_pi8(0x00), u._s);
	return v;
}

inline Vector<unsigned int, 2> operator ~(const Vector<unsigned int, 2> &u) {
	Vector<unsigned int, 2> v;
	v._s = _m_pxor(u._s, _mm_set1_pi8(0xff));
	return v;
}

template <typename... Args>
inline Vector<unsigned int, 3, Pass<unsigned int, 3> >::Vector(Args... args) {
	initialize<0, unsigned int, 3, Pass<unsigned int, 3> >(*this, args...);
}

template <>
inline Vector<unsigned int, 3, Pass<unsigned int, 3> > &Vector<unsigned int, 3, Pass<unsigned int, 3> >::operator =(const Vector<unsigned int, 3, Pass<unsigned int, 3> > &v) {
	this->_s = v._s;
	return *this;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 3, Pass<unsigned int, 3> > &Vector<unsigned int, 3, Pass<unsigned int, 3> >::operator =(const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	for (size_t c = 0; c < 3; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

inline unsigned int &Vector<unsigned int, 3, Pass<unsigned int, 3> >::operator [](size_t index) {
	Pass<unsigned int, 3> perm;
	return perm(_v, index);
}

inline const unsigned int &Vector<unsigned int, 3, Pass<unsigned int, 3> >::operator [](size_t index) const {
	Pass<unsigned int, 3> perm;
	return perm(_v, index);
}

inline const unsigned int *Vector<unsigned int, 3, Pass<unsigned int, 3> >::internal() const {
	return _v;
}

inline Vector<unsigned int, 3> &operator ++(Vector<unsigned int, 3> &self) {
	self._s = _mm_add_epi32(self._s, _mm_set1_epi32(1));
	return self;
}

inline Vector<unsigned int, 3> &operator --(Vector<unsigned int, 3> &self) {
	self._s = _mm_sub_epi32(self._s, _mm_set1_epi32(1));
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 3> &operator +=(Vector<unsigned int, 3> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<unsigned int, 3> packed(v);
	self._s = _mm_add_epi32(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 3> &operator -=(Vector<unsigned int, 3> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<unsigned int, 3> packed(v);
	self._s = _mm_sub_epi32(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 3> &operator &=(Vector<unsigned int, 3> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<unsigned int, 3> packed(v);
	self._s = _mm_and_si128(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 3> &operator |=(Vector<unsigned int, 3> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<unsigned int, 3> packed(v);
	self._s = _mm_or_si128(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 3> &operator ^=(Vector<unsigned int, 3> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	Vector<unsigned int, 3> packed(v);
	self._s = _mm_xor_si128(self._s, packed._s);
	return self;
}

inline Vector<unsigned int, 3> &operator <<=(Vector<unsigned int, 3> &self, size_t s) {
	self._s = _mm_slli_epi32(self._s, (unsigned int) s);
	return self;
}

inline Vector<unsigned int, 3> &operator >>=(Vector<unsigned int, 3> &self, size_t s) {
	self._s = _mm_srli_epi32(self._s, (unsigned int) s);
	return self;
}

inline Vector<unsigned int, 3> operator +(const Vector<unsigned int, 3> &u) {
	return u;
}

inline Vector<unsigned int, 3> operator -(const Vector<unsigned int, 3> &u) {
	Vector<unsigned int, 3> v;
	v._s = _mm_sub_epi32(_mm_set1_epi8(0x00), u._s);
	return v;
}

inline Vector<unsigned int, 3> operator ~(const Vector<unsigned int, 3> &u) {
	Vector<unsigned int, 3> v;
	v._s = _mm_xor_si128(u._s, _mm_set1_epi8(0xff));
	return v;
}

template <typename... Args>
inline Vector<unsigned int, 4, Pass<unsigned int, 4> >::Vector(Args... args) {
	initialize<0, unsigned int, 4, Pass<unsigned int, 4> >(*this, args...);
}

template <>
inline Vector<unsigned int, 4, Pass<unsigned int, 4> > &Vector<unsigned int, 4, Pass<unsigned int, 4> >::operator =(const Vector<unsigned int, 4, Pass<unsigned int, 4> > &v) {
	this->_s = v._s;
	return *this;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 4, Pass<unsigned int, 4> > &Vector<unsigned int, 4, Pass<unsigned int, 4> >::operator =(const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	for (size_t c = 0; c < 4; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

inline unsigned int &Vector<unsigned int, 4, Pass<unsigned int, 4> >::operator [](size_t index) {
	Pass<unsigned int, 4> perm;
	return perm(_v, index);
}

inline const unsigned int &Vector<unsigned int, 4, Pass<unsigned int, 4> >::operator [](size_t index) const {
	Pass<unsigned int, 4> perm;
	return perm(_v, index);
}

inline const unsigned int *Vector<unsigned int, 4, Pass<unsigned int, 4> >::internal() const {
	return _v;
}

inline Vector<unsigned int, 4> &operator ++(Vector<unsigned int, 4> &self) {
	self._s = _mm_add_epi32(self._s, _mm_set1_epi32(1));
	return self;
}

inline Vector<unsigned int, 4> &operator --(Vector<unsigned int, 4> &self) {
	self._s = _mm_sub_epi32(self._s, _mm_set1_epi32(1));
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 4> &operator +=(Vector<unsigned int, 4> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<unsigned int, 4> packed(v);
	self._s = _mm_add_epi32(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 4> &operator -=(Vector<unsigned int, 4> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<unsigned int, 4> packed(v);
	self._s = _mm_sub_epi32(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 4> &operator &=(Vector<unsigned int, 4> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<unsigned int, 4> packed(v);
	self._s = _mm_and_si128(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 4> &operator |=(Vector<unsigned int, 4> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<unsigned int, 4> packed(v);
	self._s = _mm_or_si128(self._s, packed._s);
	return self;
}

template <size_t SizeV, typename PermutationV>
inline Vector<unsigned int, 4> &operator ^=(Vector<unsigned int, 4> &self, const Vector<unsigned int, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	Vector<unsigned int, 4> packed(v);
	self._s = _mm_xor_si128(self._s, packed._s);
	return self;
}

inline Vector<unsigned int, 4> &operator <<=(Vector<unsigned int, 4> &self, size_t s) {
	self._s = _mm_slli_epi32(self._s, (unsigned int) s);
	return self;
}

inline Vector<unsigned int, 4> &operator >>=(Vector<unsigned int, 4> &self, size_t s) {
	self._s = _mm_srli_epi32(self._s, (unsigned int) s);
	return self;
}

inline Vector<unsigned int, 4> operator +(const Vector<unsigned int, 4> &u) {
	return u;
}

inline Vector<unsigned int, 4> operator -(const Vector<unsigned int, 4> &u) {
	Vector<unsigned int, 4> v;
	v._s = _mm_sub_epi32(_mm_set1_epi8(0x00), u._s);
	return v;
}

inline Vector<unsigned int, 4> operator ~(const Vector<unsigned int, 4> &u) {
	Vector<unsigned int, 4> v;
	v._s = _mm_xor_si128(u._s, _mm_set1_epi8(0xff));
	return v;
}

#if 0
template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY>
inline Vector<Type, Permutation::Elements> pow(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y) {
	static_assert(Permutation::Elements == PermutationY::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = pow(x[c], y[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> exp(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = exp(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> log(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = log(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> abs(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = abs(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> floor(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = floor(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> ceil(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = ceil(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY>
inline Vector<Type, Permutation::Elements> min(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y) {
	static_assert(Permutation::Elements == PermutationY::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = min(x[c], y[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY>
inline Vector<Type, Permutation::Elements> max(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y) {
	static_assert(Permutation::Elements == PermutationY::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = max(x[c], y[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> clamp(const Vector<Type, Size, Permutation> &x, const Type &minVal, const Type &maxVal) {
	return min(max(x, Vector<Type, Permutation::Elements>(minVal)), Vector<Type, Permutation::Elements>(maxVal));
}
#endif

#if 0
template <size_t Size, typename Permutation>
inline unsigned int packHalf2x16(const Vector<float, Size, Permutation> &v) {
	static_assert(Permutation::Elements == 2, "Vector must have exactly two elements");
	Vector<float, 2> packed(v);
	__m64 bits = _mm_movepi64_pi64(_mm_castps_si128(packed._s));
	__m64 sign = _m_pand(_m_psrldi(bits, 16), _mm_set1_pi32(0xc000));
	__m64 expvalue = _m_pand(_m_psrldi(bits, 13), _mm_set1_pi32(0x3fff));
	__m64 combine = _mm_shuffle_pi16(_m_por(sign, expvalue), (0 << 0) | (2 << 2) | (0 << 4) | (2 << 6));
	return (unsigned int) _m_to_int(combine);
}

inline Vector<float, 2> unpackHalf2x16(unsigned int v) {
	__m64 bits = _mm_shuffle_pi16(_m_from_int((int) v), (0 << 0) | (2 << 2) | (1 << 4) | (3 << 6));
	__m64 sign = _m_por(bits, _mm_set1_pi32(0x8000));
	__m64 biasexp = _m_por(bits, _mm_set1_pi32(0x7c00));
	__m64 value = _m_por(bits, _mm_set1_pi32(0x03ff));
	__m64 ezero = _m_pcmpeqd(biasexp, _mm_set1_pi32(0x0000));
	__m64 emax = _m_pcmpeqd(biasexp, _mm_set1_pi32(0x7c00));
	__m64 eother = _m_pxor(_m_por(ezero, emax), _mm_set1_pi8(0xff));
	__m64 exp = _m_pxor(_m_pand(emax, _mm_set1_pi32(0x3fc00)), _m_pand(eother, _mm_set1_pi32(0x1c000)));
	__m64 combine = _m_por(_m_pslldi(sign, 16), _m_pslldi(_m_por(value, exp), 13));
	Vector<float, 2> ret;
	ret._s = _mm_castsi128_ps(_mm_movpi64_epi64(combine));
	return ret;
}
#endif

}


// Load template definitions

#define GLAM_VECTOR_TEMPLATE_IMPL
#include <glam/vector.tpp>
#undef GLAM_VECTOR_TEMPLATE_IMPL


#else //GLAM_VECTOR_SSE

// Fallback

#include <glam/vector.h>

#endif //GLAM_VECTOR_SSE

#endif //GLAM_VECTOR_SSE_H