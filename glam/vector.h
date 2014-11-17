/*
 * GLAM - GLSL Linear Algebra Math Library
 * 
 * Copyright (c) 2013-2014, Gregor Riepl <onitake@gmail.com>
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

#ifndef GLAM_VECTOR_H
#define GLAM_VECTOR_H

#include <iostream>
#include <cstddef>
#include <glam/math.h>
#include <glam/config.h>
#include <glam/exception.h>

namespace glam {

// Classes

// Passthrough member access
template <typename Type, size_t Size>
struct Pass {
	// The number of used elements of the underlying vector (equal to its size)
	static constexpr size_t Elements = Size;
	
	// Access element index of array v
	inline Type &operator ()(Type *v, size_t index) const;
	// Constant access element index of array v
	inline const Type &operator ()(const Type *v, size_t index) const;
};

// Member access with permutation
template <typename Type, size_t Size, size_t... Permutation>
struct Permutatator {
	// The number of used elements of the underlying vector (not its memory requirements/size)
	static constexpr size_t Elements = sizeof...(Permutation);
	
	// Access element index of array v
	inline Type &operator ()(Type *v, size_t index) const;
	// Constant access element index of array v
	inline const Type &operator ()(const Type *v, size_t index) const;
};

template <typename Type, size_t Size, typename Permutation = Pass<Type, Size> >
class Vector {
private:
	Type _v[Size];
	
public:
	// Catch-all variadic constructor
	// The following constructor calls are supported:
	// Vector(): Creates an empty vector
	//  All elements are 0.
	// Vector(const Type &s): Extend scalar to vector
	//  All elements are initialized to the same value.
	// Vector<ForwardIterator>(std::pair<ForwardIterator, ForwardIterator> iterator): Populate the vector from an iterator pair
	//  Elements are filled up in order of occurrence, from first to last element (ex. x, y, z, w).
	//  std::pair is used to avoid call ambiguity problems. Call std::make_pair(begin, end) to constuct your iterator pair.
	// Vector<AnyScalarOrVector...>(const AnyScalarOrVector &arg0, ...): Populate the vector from a list of vectors and/or scalars
	//  Elements are filled up in order of occurrence.
	// Vector(const Type *m): Populate the vector from a constant array
	//  Elements are filled up in order of occurrence.
	//  Note that this overload does not support range checking. Make sure that the input array has the correct size or use the ForwardIterator overload.
	// Vector<OtherType>(const Vector<OtherType> &arg0): Convert a vector of different element type
	//  The lengths of the vectors need to match.
	template <typename... Args>
	inline Vector(Args... args);

	// Assignment operator
	template <size_t SizeV, typename PermutationV>
	inline Vector<Type, Size, Permutation> &operator =(const Vector<Type, SizeV, PermutationV> &v);

	// Member access operator
	inline Type &operator [](size_t index);
	// Const member access operator
	inline const Type &operator [](size_t index) const;
};


// Types

// Float vectors
typedef Vector<float, 2> vec2;
typedef Vector<float, 3> vec3;
typedef Vector<float, 4> vec4;
// Int vectors
typedef Vector<int, 2> ivec2;
typedef Vector<int, 3> ivec3;
typedef Vector<int, 4> ivec4;
// Double vectors
typedef Vector<double, 2> dvec2;
typedef Vector<double, 3> dvec3;
typedef Vector<double, 4> dvec4;
// Boolean vectors
typedef Vector<bool, 2> bvec2;
typedef Vector<bool, 3> bvec3;
typedef Vector<bool, 4> bvec4;
// Unsigned Int vectors
typedef Vector<unsigned int, 2> uvec2;
typedef Vector<unsigned int, 3> uvec3;
typedef Vector<unsigned int, 4> uvec4;


// Macros

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

template <typename Type>
class Vector<Type, 2, Pass<Type, 2> > {
public:
	union {
		Type _v[2];
		// Produces permutations: x, y, xx, xy, yx, yy
		GLAM_PERM_MAKE2_0(Type, x, y);
		// Produces permutations: r, g, rr, rg, gr, gg
		GLAM_PERM_MAKE2_0(Type, r, g);
		// Produces permutations: s, t, ss, st, ts tt
		GLAM_PERM_MAKE2_0(Type, s, t);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<Type, 2, Pass<Type, 2> > &operator =(const Vector<Type, SizeV, PermutationV> &v);
	
	inline Type &operator [](size_t index);
	inline const Type &operator [](size_t index) const;
	
	// Direct access to internal pointer
	inline const Type *internal() const;
};

template <typename Type>
class Vector<Type, 3, Pass<Type, 3> > {
public:
	union {
		Type _v[3];
		// Produces permutations: x, y, z, xx, xy, xz, yx, yy, yz, zx, zy, zz, xxx, xxy, xxz, xyx, ...
		GLAM_PERM_MAKE3_0(Type, x, y, z);
		// Produces permutations: r, g, b, rr, rg, rb, gr, gg, gb, br, bg, bb, rrr, rrg, rrb, rgr, ...
		GLAM_PERM_MAKE3_0(Type, r, g, b);
		// Produces permutations: s, t, p, ss, st, sp, ts, tt, tp, ps, pt, pp, sss, sst, ssp, sts, ...
		GLAM_PERM_MAKE3_0(Type, s, t, p);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<Type, 3, Pass<Type, 3> > &operator =(const Vector<Type, SizeV, PermutationV> &v);
	
	inline Type &operator [](size_t index);
	inline const Type &operator [](size_t index) const;

	// Direct access to internal pointer
	inline const Type *internal() const;
};

template <typename Type>
class Vector<Type, 4, Pass<Type, 4> > {
public:
	union {
		Type _v[4];
		// Produces permutations: x, y, z, w, xx, xy, xz, xw, yx, yy, yz, yw, zx, zy, zz, zw, wx, wy, wz, ww, xxx, xxy, xxz, xyx, ..., xxxx, xxxy, ...
		GLAM_PERM_MAKE4_0(Type, x, y, z, w);
		// Produces permutations: r, g, b, a, rr, rg, rb, ra, gr, gg, gb, ga, br, bg, bb, ba, ar, ag, ab, aa, rrr, rrg, rrb, rgr, ..., rrrr, rrrg, ...
		GLAM_PERM_MAKE4_0(Type, r, g, b, a);
		// Produces permutations: s, t, p, q, ss, st, sp, sq, ts, tt, tp, tq, ps, pt, pp, pq, qs, qt, qp, qq, sss, sst, ssp, sts, ..., ssss, ssst, ...
		GLAM_PERM_MAKE4_0(Type, s, t, p, q);
	};
	
	template <typename... Args>
	inline Vector(Args... args);
	
	template <size_t SizeV, typename PermutationV>
	inline Vector<Type, 4, Pass<Type, 4> > &operator =(const Vector<Type, SizeV, PermutationV> &v);
	
	inline Type &operator [](size_t index);
	inline const Type &operator [](size_t index) const;
	
	// Direct access to internal pointer
	inline const Type *internal() const;
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


// Declarations

// Deferred constructors
template <size_t Index, typename Type, size_t Size, typename Permutation>
inline void initialize(Vector<Type, Size, Permutation> &self);
template <size_t Index, typename Type, size_t Size, typename Permutation, typename ForwardIterator>
inline void initialize(Vector<Type, Size, Permutation> &self, std::pair<ForwardIterator, ForwardIterator> iterator);
template <size_t Index, typename Type, size_t Size, typename Permutation>
inline void initialize(Vector<Type, Size, Permutation> &self, const Type *v);
template <size_t Index, typename Type, size_t Size, typename Permutation, typename TypeV, size_t SizeV, typename PermutationV>
inline void initialize(Vector<Type, Size, Permutation> &self, const Vector<TypeV, SizeV, PermutationV> &v);
template <size_t Index, typename Type, size_t Size, typename Permutation>
inline void initialize(Vector<Type, Size, Permutation> &self, const Type &s);
template <size_t Index, typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV, typename... Args>
inline void initialize(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v, Args... args);
template <size_t Index, typename Type, size_t Size, typename Permutation, typename... Args>
inline void initialize(Vector<Type, Size, Permutation> &self, const Type &s, Args... args);

// Component-wise pre-increment operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator ++(Vector<Type, Size, Permutation> &self);
// Component-wise post-increment operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> operator ++(Vector<Type, Size, Permutation> &self, int);
// Component-wise pre-decrement operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator --(Vector<Type, Size, Permutation> &self);
// Component-wise post-decrement operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> operator --(Vector<Type, Size, Permutation> &self, int);

// Component-wise self-addition operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator +=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v);
// Vector-scalar addition, the same value is added to all components
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator +=(Vector<Type, Size, Permutation> &self, const Type &v);
// Component-wise self-subtraction operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator -=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v);
// Vector-scalar subtraction, the same value is taken from all components
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator -=(Vector<Type, Size, Permutation> &self, const Type &v);
// Component-wise self-multiplication operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator *=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v);
// Vector-scalar product, all components are multiplied by the same value
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator *=(Vector<Type, Size, Permutation> &self, const Type &v);
// Component-wise self-division operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator /=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v);
// Vector-scalar division, all components are divided by the same value
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator /=(Vector<Type, Size, Permutation> &self, const Type &v);
// Component-wise self-modulus operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator %=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v);
// Vector-scalar modulo division, all components are replaced by the remainder of their division by the same value
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator %=(Vector<Type, Size, Permutation> &self, const Type &v);

// Component-wise bitwise and operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator &=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v);
// Vector-scalar bitwise and, all components are masked by the same bit pattern
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator &=(Vector<Type, Size, Permutation> &self, const Type &v);
// Component-wise bitwise or operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator |=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v);
// Vector-scalar bitwise or, all components are combined with the same bit pattern
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator |=(Vector<Type, Size, Permutation> &self, const Type &v);
// Component-wise bitwise xor operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator ^=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v);
// Vector-scalar bitwise xor, all components are flipped or kept depending on the same bit pattern
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator ^=(Vector<Type, Size, Permutation> &self, const Type &v);
// Component-wise bit shift left operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator <<=(Vector<Type, Size, Permutation> &self, size_t s);
// Component-wise bit shift right operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator >>=(Vector<Type, Size, Permutation> &self, size_t s);

// Component-wise unary plus operator (no-op if T follows standard algebra)
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator +(const Vector<Type, Size, Permutation> &u);
// Component-wise negation operator (negation)
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator -(const Vector<Type, Size, Permutation> &u);
// Component-wise ones complement operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator ~(const Vector<Type, Size, Permutation> &u);

// Component-wise addition operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator +(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise vector-scalar addition operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator +(const Vector<Type, Size, Permutation> &u, const Type &v);
// Component-wise scalar-vector addition operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator +(const Type &u, const Vector<Type, Size, Permutation> &v);
// Component-wise subtraction operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator -(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise vector-scalar subtraction operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator -(const Vector<Type, Size, Permutation> &u, const Type &v);
// Component-wise scalar-vector subtraction operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator -(const Type &u, const Vector<Type, Size, Permutation> &v);
// Component-wise multiplication operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator *(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise vector-scalar multiplication operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator *(const Vector<Type, Size, Permutation> &u, const Type &v);
// Component-wise scalar-vector multiplication operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator *(const Type &u, const Vector<Type, Size, Permutation> &v);
// Component-wise division operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator /(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise vector-scalar division operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator /(const Vector<Type, Size, Permutation> &u, const Type &v);
// Component-wise scalar-vector division operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator /(const Type &u, const Vector<Type, Size, Permutation> &v);
// Component-wise modulo operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator %(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise vector-scalar modulo operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator %(const Vector<Type, Size, Permutation> &u, const Type &v);
// Component-wise scalar-vector modulo operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator %(const Type &u, const Vector<Type, Size, Permutation> &v);

// Component-wise bitwise and operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator &(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise vector-scalar bitwise and operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator &(const Vector<Type, Size, Permutation> &u, const Type &v);
// Component-wise scalar-vector bitwise and operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator &(const Type &u, const Vector<Type, Size, Permutation> &v);
// Component-wise bitwise or operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator |(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise vector-scalar bitwise or operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator |(const Vector<Type, Size, Permutation> &u, const Type &v);
// Component-wise scalar-vector bitwise or operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator |(const Type &u, const Vector<Type, Size, Permutation> &v);
// Component-wise bitwise xor operator
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator ^(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise vector-scalar bitwise xor operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator ^(const Vector<Type, Size, Permutation> &u, const Type &v);
// Component-wise scalar-vector bitwise xor operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator ^(const Type &u, const Vector<Type, Size, Permutation> &v);
// Component-wise bit shift left operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator <<(const Vector<Type, Size, Permutation> &u, size_t s);
// Component-wise bit shift right operator
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator >>(const Vector<Type, Size, Permutation> &u, size_t s);

// Component-wise equality comparison operator (exact compare)
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> equal(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise inequality comparison operator (exact compare)
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> notEqual(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise less-than comparison operator (exact compare)
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> lessThan(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise greater-than comparison operator (exact compare)
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> greaterThan(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise less-or-equal comparison operator (exact compare)
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> lessThanEqual(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Component-wise greater-or-equal comparison operator (exact compare)
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> greaterThanEqual(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);

// Vector comparison operator (exact compare), correspondes to all(equal(u, v))
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline bool operator ==(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);
// Vector comparison operator (exact compare), correspondes to any(notEqual(u, v))
template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline bool operator !=(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v);

// Component-wise sine
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> sin(const Vector<Type, Size, Permutation> &v);
// Component-wise cosine
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> cos(const Vector<Type, Size, Permutation> &v);
// Component-wise tangent
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> tan(const Vector<Type, Size, Permutation> &v);
// Component-wise arc sine
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> asin(const Vector<Type, Size, Permutation> &v);
// Component-wise arc cosine
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> acos(const Vector<Type, Size, Permutation> &v);
// Component-wise arc tangent
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> atan(const Vector<Type, Size, Permutation> &v);
// Component-wise sinus hyperbolicus
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> sinh(const Vector<Type, Size, Permutation> &v);
// Component-wise cosinus hyperbolicus
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> cosh(const Vector<Type, Size, Permutation> &v);
// Component-wise tangens hyperbolicus
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> tanh(const Vector<Type, Size, Permutation> &v);

// Component-wise power
template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY>
inline Vector<Type, Permutation::Elements> pow(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y);
// Component-wise exponential (e^v)
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> exp(const Vector<Type, Size, Permutation> &v);
// Component-wise natural logarithm
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> log(const Vector<Type, Size, Permutation> &v);
// Component-wise power of 2
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> exp2(const Vector<Type, Size, Permutation> &v);
// Component-wise base 2 logarithm
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> log2(const Vector<Type, Size, Permutation> &v);
// Component-wise square root
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> sqrt(const Vector<Type, Size, Permutation> &v);

// Component-wise absolute value
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> abs(const Vector<Type, Size, Permutation> &v);
// Component-wise sign
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> sign(const Vector<Type, Size, Permutation> &v);
// Component-wise round towards negative infinity
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> floor(const Vector<Type, Size, Permutation> &v);
// Component-wise truncate
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> trunc(const Vector<Type, Size, Permutation> &v);
// Component-wise round to the nearest integer
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> round(const Vector<Type, Size, Permutation> &v);
// Component-wise round to the nearest even integer
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> roundEven(const Vector<Type, Size, Permutation> &v);
// Component-wise round towards positive infinity
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> ceil(const Vector<Type, Size, Permutation> &v);
// Component-wise fractional part separation
template <typename Type, size_t Size, typename Permutation, size_t SizeI, typename PermutationI>
inline Vector<Type, Permutation::Elements> modf(const Vector<Type, Size, Permutation> &v, Vector<Type, SizeI, PermutationI> &i);
// Component-wise minimum
template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY>
inline Vector<Type, Permutation::Elements> min(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y);
// Component-wise maximum
template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY>
inline Vector<Type, Permutation::Elements> max(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y);
// Saturation with single min/max values (applied to all vector components)
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> clamp(const Vector<Type, Size, Permutation> &x, const Type &minVal, const Type &maxVal);
// Component-wise linear blend (x * (1 - a) + y * a)
template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY>
inline Vector<Type, Permutation::Elements> mix(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y, const Type &a);
// Return x for each component of a that is false, y if true
template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY, size_t SizeA, typename PermutationA>
inline Vector<Type, Permutation::Elements> mix(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y, const Vector<bool, SizeA, PermutationA> &a);
// Component-wise step transition
template <typename Type, size_t Size, typename Permutation, size_t SizeX, typename PermutationX>
inline Vector<Type, Permutation::Elements> step(const Vector<Type, Size, Permutation> &edge, const Vector<Type, SizeX, PermutationX> &x);
// Component-wise step transition with common edge argument
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> step(const Type &edge, const Vector<Type, Size, Permutation> &x);
// Component-wise smooth step transition
template <typename Type, size_t Size, typename Permutation, size_t Size1, typename Permutation1, size_t SizeX, typename PermutationX>
inline Vector<Type, Permutation::Elements> smoothstep(const Vector<Type, Size, Permutation> &edge0, const Vector<Type, Size1, Permutation1> &edge1, const Vector<Type, SizeX, PermutationX> &x);
// Component-wise smooth step transition with common edge arguments
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> smoothstep(const Type &edge0, const Type &edge1, const Vector<Type, Size, Permutation> &x);
// Component-wise test for NaN
template <typename Type, size_t Size, typename Permutation>
inline Vector<bool, Permutation::Elements> isnan(const Vector<Type, Size, Permutation> &v);
// Component-wise test for infinity
template <typename Type, size_t Size, typename Permutation>
inline Vector<bool, Permutation::Elements> isinf(const Vector<Type, Size, Permutation> &v);

// Vector length (sqrt(dot(a, a)))
template <typename Type, size_t Size, typename Permutation>
inline Type length(const Vector<Type, Size, Permutation> &x);
// Distance between two vectors (length(b - a))
template <typename Type, size_t Size, typename Permutation, size_t Size1, typename Permutation1>
inline Type distance(const Vector<Type, Size, Permutation> &p0, const Vector<Type, Size1, Permutation1> &p1);
// Dot product of two vectors (a[0] * b[0] + a[1] * b[1] + ...)
template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY>
inline Type dot(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y);
// Cross product of two 3-component vectors
template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY>
inline Vector<Type, 3> cross(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y);
// Unit length vector (a / length(a))
// always_inline encouraged to keep gcc from generating IEEE correct math
template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> normalize(const Vector<Type, Size, Permutation> &x) __attribute__((always_inline));
// Check if I and Nref are facing in opposite directions and return N if yes, -N otherwise (dot(Nref, I) < 0 ? N : -N)
template <class Type>
inline Type faceforward(const Type &N, const Type &I, const Type &Nref);
// Reflect I in respect to surface normal N
template <typename Type, size_t Size, typename Permutation, size_t SizeN, typename PermutationN>
inline Vector<Type, Permutation::Elements> reflect(const Vector<Type, Size, Permutation> &I, const Vector<Type, SizeN, PermutationN> &N);
// Refract I on the surface defined by normal N using index of refraction ration eta
template <typename Type, size_t Size, typename Permutation, size_t SizeN, typename PermutationN>
inline Vector<Type, Permutation::Elements> refract(const Vector<Type, Size, Permutation> &I, const Vector<Type, SizeN, PermutationN> &N, const Type &eta);

// True if any component of x is true (or evaluates to true)
template <size_t Size, typename Permutation>
inline bool any(const Vector<bool, Size, Permutation> &v);
// True if all components of x are true (or evaluates to true)
template <size_t Size, typename Permutation>
inline bool all(const Vector<bool, Size, Permutation> &v);
// Component-wise negation of a boolean vector, this is an overloaded ! operator because
// 'not' is a reserved word in C++ with the same semantics as !, so can be used interchangeably
template <size_t Size, typename Permutation>
inline Vector<bool, Permutation::Elements> operator !(const Vector<bool, Size, Permutation> &v);

// Formatted output operator (vectors will be presented in the form "(0.1,1.50)")
template <typename Type, size_t Size, typename Permutation>
inline std::ostream &operator <<(std::ostream &stream, const Vector<Type, Size, Permutation> &v);
// Read vector components from stream, in text form, separated by whitespace
template <typename Type, size_t Size, typename Permutation>
inline std::istream &operator >>(std::istream &stream, const Vector<Type, Size, Permutation> &v);

// Convert a floating point vector into an integer vector of equal size and exactly the same bit pattern
template <size_t Size, typename Permutation>
inline Vector<int, Permutation::Elements> floatBitsToInt(const Vector<float, Size, Permutation> &v);
template <size_t Size, typename Permutation>
inline Vector<unsigned int, Permutation::Elements> floatBitsToUInt(const Vector<float, Size, Permutation> &v);
// Convert an integer vector into a floating point vector of equal size and exactly the same bit pattern
template <size_t Size, typename Permutation>
inline Vector<float, Permutation::Elements> intBitsToFloat(const Vector<int, Size, Permutation> &v);
template <size_t Size, typename Permutation>
inline Vector<float, Permutation::Elements> uintBitsToFloat(const Vector<unsigned int, Size, Permutation> &v);

// Convert two normalized floats into 16bit fixed point values, then pack them into a 32bit integer
// The first vector component will be at the lower 16 bits of the result, the second at the upper part.
// Conversion from float to fixed is equivalent to: round(clamp(c, -1, +1) * 32767.0)
template <size_t Size, typename Permutation>
inline unsigned int packSnorm2x16(const Vector<float, Size, Permutation> &v);
// Convert two packed 16bit fixed point values into a vector of normalized floats
// The first vector component is taken from the lower 16 bits, the second from the upper part.
// Conversion from fixed to float is equivalent to: clamp(f / 32767.0, -1, +1)
inline Vector<float, 2> unpackSnorm2x16(unsigned int p);
// Convert two normalized floats into 16bit fixed point values, then pack them into a 32bit integer
// The first vector component will be at the lower 16 bits of the result, the second at the upper part.
// Conversion from float to fixed is equivalent to: round(clamp(c, 0, +1) * 65535.0)
template <size_t Size, typename Permutation>
inline unsigned int packUnorm2x16(const Vector<float, Size, Permutation> &v);
// Convert two packed 16bit fixed point values into a vector of normalized floats
// The first vector component is taken from the lower 16 bits, the second from the upper part.
// Conversion from fixed to float is equivalent to: f / 65535.0
inline Vector<float, 2> unpackUnorm2x16(unsigned int p);
// Convert two 32bit floats into half (16bit) floats and pack them into a 32bit integer
// The first vector component will be at the lower 16 bits of the result, the second at the upper part.
template <size_t Size, typename Permutation>
inline unsigned int packHalf2x16(const Vector<float, Size, Permutation> &v);
// Convert two half (16bit) floats packed into a 32bit integer into a vector of 32bit floats
// The first vector component is taken from the lower 16 bits, the second from the upper part.
inline Vector<float, 2> unpackHalf2x16(unsigned int v);


// Implementation

template <typename Type, size_t Size>
inline Type &Pass<Type, Size>::operator ()(Type *v, size_t index) const {
#ifdef GLAM_RANGE_CHECK
	if (index >= Elements) {
		throw DimensionOutOfRangeException("Invalid vector element index");
	}
#endif
	return v[index];
}

template <typename Type, size_t Size>
inline const Type &Pass<Type, Size>::operator ()(const Type *v, size_t index) const {
#ifdef GLAM_RANGE_CHECK
	if (index >= Elements) {
		throw DimensionOutOfRangeException("Invalid vector element index");
	}
#endif
	return v[index];
}

template <typename Type, size_t Size, size_t... Permutation>
inline Type &Permutatator<Type, Size, Permutation...>::operator ()(Type *v, size_t index) const {
#ifdef GLAM_RANGE_CHECK
	if (index >= Elements) {
		throw DimensionOutOfRangeException("Invalid permutation element index", Elements - 1, index);
	}
#endif
	constexpr size_t p[Elements] = { Permutation... };
	if (index < Elements) {
		return v[p[index]];
	}
	return v[index];
}

template <typename Type, size_t Size, size_t... Permutation>
inline const Type &Permutatator<Type, Size, Permutation...>::operator ()(const Type *v, size_t index) const {
#ifdef GLAM_RANGE_CHECK
	if (index >= Elements) {
		throw DimensionOutOfRangeException("Invalid permutation element index", Elements - 1, index);
	}
#endif
	constexpr size_t p[Elements] = { Permutation... };
	if (index < Elements) {
		return v[p[index]];
	}
	return v[index];
}

template <typename Type>
template <typename... Args>
inline Vector<Type, 2, Pass<Type, 2> >::Vector(Args... args) {
	initialize<0, Type, 2, Pass<Type, 2> >(*this, args...);
}

template <typename Type>
template <size_t SizeV, typename PermutationV>
inline Vector<Type, 2, Pass<Type, 2> > &Vector<Type, 2, Pass<Type, 2> >::operator =(const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 2, "Vectors differ in number of elements");
	for (size_t c = 0; c < 2; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

template <typename Type>
inline Type &Vector<Type, 2, Pass<Type, 2> >::operator [](size_t index) {
	Pass<Type, 2> perm;
	return perm(_v, index);
}

template <typename Type>
inline const Type &Vector<Type, 2, Pass<Type, 2> >::operator [](size_t index) const {
	Pass<Type, 2> perm;
	return perm(_v, index);
}

template <typename Type>
inline const Type *Vector<Type, 2, Pass<Type, 2> >::internal() const {
	return _v;
}

template <typename Type>
template <typename... Args>
inline Vector<Type, 3, Pass<Type, 3> >::Vector(Args... args) {
	initialize<0, Type, 3, Pass<Type, 3> >(*this, args...);
}

template <typename Type>
template <size_t SizeV, typename PermutationV>
inline Vector<Type, 3, Pass<Type, 3> > &Vector<Type, 3, Pass<Type, 3> >::operator =(const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 3, "Vectors differ in number of elements");
	for (size_t c = 0; c < 3; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

template <typename Type>
inline Type &Vector<Type, 3, Pass<Type, 3> >::operator [](size_t index) {
	Pass<Type, 3> perm;
	return perm(_v, index);
}

template <typename Type>
inline const Type &Vector<Type, 3, Pass<Type, 3> >::operator [](size_t index) const {
	Pass<Type, 3> perm;
	return perm(_v, index);
}

template <typename Type>
inline const Type *Vector<Type, 3, Pass<Type, 3> >::internal() const {
	return _v;
}

template <typename Type>
template <typename... Args>
inline Vector<Type, 4, Pass<Type, 4> >::Vector(Args... args) {
	initialize<0, Type, 4, Pass<Type, 4> >(*this, args...);
}

template <typename Type>
template <size_t SizeV, typename PermutationV>
inline Vector<Type, 4, Pass<Type, 4> > &Vector<Type, 4, Pass<Type, 4> >::operator =(const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(PermutationV::Elements == 4, "Vectors differ in number of elements");
	for (size_t c = 0; c < 4; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

template <typename Type>
inline Type &Vector<Type, 4, Pass<Type, 4> >::operator [](size_t index) {
	Pass<Type, 4> perm;
	return perm(_v, index);
}

template <typename Type>
inline const Type &Vector<Type, 4, Pass<Type, 4> >::operator [](size_t index) const {
	Pass<Type, 4> perm;
	return perm(_v, index);
}

template <typename Type>
inline const Type *Vector<Type, 4, Pass<Type, 4> >::internal() const {
	return _v;
}

template <size_t Index, typename Type, size_t Size, typename Permutation>
inline void initialize(Vector<Type, Size, Permutation> &self) {
	static_assert(Index == 0 || Index >= Permutation::Elements, "Not enough initializers");
	if (Index == 0) {
		for (size_t p = 0; p < Permutation::Elements; p++) {
			self[p] = Type(0);
		}
	}
}

template <size_t Index, typename Type, size_t Size, typename Permutation, typename ForwardIterator>
inline void initialize(Vector<Type, Size, Permutation> &self, std::pair<ForwardIterator, ForwardIterator> iterator) {
	static_assert(Index == 0, "No other arguments are allowed when initializing a vector from an interator");
	size_t p = Index;
	for (ForwardIterator it = iterator.first; p < Permutation::Elements && it != iterator.second; p++, it++) {
		self[p] = *it;
	}
}

template <size_t Index, typename Type, size_t Size, typename Permutation>
inline void initialize(Vector<Type, Size, Permutation> &self, const Type *v) {
	static_assert(Index == 0, "No other arguments are allowed when initializing a vector from an array");
	for (size_t p = 0; p < Permutation::Elements; p++) {
		self[Index + p] = v[p];
	}
}

template <size_t Index, typename Type, size_t Size, typename Permutation, typename TypeV, size_t SizeV, typename PermutationV>
inline void initialize(Vector<Type, Size, Permutation> &self, const Vector<TypeV, SizeV, PermutationV> &v) {
	static_assert(Index == 0, "No other arguments are allowed when initializing a vector from a vector of different element type");
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	for (size_t p = 0; p < Permutation::Elements; p++) {
		self[Index + p] = Type(v[p]);
	}
}

template <size_t Index, typename Type, size_t Size, typename Permutation>
inline void initialize(Vector<Type, Size, Permutation> &self, const Type &s) {
	static_assert(Index + 1 <= Permutation::Elements, "Too many initializers");
	static_assert(Index == 0 || Index + 1 == Permutation::Elements, "Not enough initializers");
	if (Index == 0) {
		for (size_t p = 0; p < Permutation::Elements; p++) {
			self[p] = s;
		}
	} else {
		self[Index] = s;
	}
}

template <size_t Index, typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV, typename... Args>
inline void initialize(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v, Args... args) {
	// The last argument may have excess components, which are ignored
	//static_assert(Index + PermutationV::Elements <= Permutation::Elements, "Too many initializers");
	constexpr size_t m = PermutationV::Elements < Permutation::Elements - Index + 1 ? PermutationV::Elements : Permutation::Elements - Index + 1;
	for (size_t p = 0; p < m; p++) {
		self[Index + p] = v[p];
	}
	initialize<Index + PermutationV::Elements, Type, Size, Permutation>(self, args...);
}

template <size_t Index, typename Type, size_t Size, typename Permutation, typename... Args>
inline void initialize(Vector<Type, Size, Permutation> &self, const Type &s, Args... args) {
	// The condition is < and not <= because there is a dedicated terminal initializer for a single scalar
	static_assert(Index + 1 < Permutation::Elements, "Too many initializers");
	self[Index] = s;
	initialize<Index + 1, Type, Size, Permutation>(self, args...);
}


template <typename Type, size_t Size, typename Permutation>
template <typename... Args>
inline Vector<Type, Size, Permutation>::Vector(Args... args) {
	initialize<0, Type, Size, Permutation>(*this, args...);
}

template <typename Type, size_t Size, typename Permutation>
inline Type &Vector<Type, Size, Permutation>::operator [](size_t index) {
	Permutation perm;
	return perm(_v, index);
}

template <typename Type, size_t Size, typename Permutation>
inline const Type &Vector<Type, Size, Permutation>::operator [](size_t index) const {
	Permutation perm;
	return perm(_v, index);
}

template <typename Type, size_t Size, typename Permutation>
template <size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &Vector<Type, Size, Permutation>::operator =(const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	for (size_t c = 0; c < Permutation::Elements; c++) {
		(*this)[c] = v[c];
	}
	return *this;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator ++(Vector<Type, Size, Permutation> &self) {
	for (size_t c = 0; c < Permutation::Elements; c++) {
		++self[c];
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> operator ++(Vector<Type, Size, Permutation> &self, int) {
	Vector<Type, Size, Permutation> temp(self);
	++self;
	return temp;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator --(Vector<Type, Size, Permutation> &self) {
	for (size_t c = 0; c < Permutation::Elements; c++) {
		--self[c];
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> operator --(Vector<Type, Size, Permutation> &self, int) {
	Vector<Type, Size, Permutation> temp(self);
	--self;
	return temp;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator +=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	for (size_t c = 0; c < Permutation::Elements; c++) {
		self[c] += v[c];
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator -=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	for (size_t c = 0; c < Permutation::Elements; c++) {
		self[c] -= v[c];
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator *=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	for (size_t c = 0; c < Permutation::Elements; c++) {
		self[c] *= v[c];
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator /=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	for (size_t c = 0; c < Permutation::Elements; c++) {
		self[c] /= v[c];
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator %=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	for (size_t c = 0; c < Permutation::Elements; c++) {
		self[c] %= v[c];
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator &=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	for (size_t c = 0; c < Permutation::Elements; c++) {
		self[c] &= v[c];
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator |=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	for (size_t c = 0; c < Permutation::Elements; c++) {
		self[c] |= v[c];
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Size, Permutation> &operator ^=(Vector<Type, Size, Permutation> &self, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	for (size_t c = 0; c < Permutation::Elements; c++) {
		self[c] ^= v[c];
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator <<=(Vector<Type, Size, Permutation> &self, size_t s) {
	for (size_t c = 0; c < Permutation::Elements; c++) {
		self[c] <<= s;
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator >>=(Vector<Type, Size, Permutation> &self, size_t s) {
	for (size_t c = 0; c < Permutation::Elements; c++) {
		self[c] >>= s;
	}
	return self;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator +(const Vector<Type, Size, Permutation> &u) {
	Vector<Type, Permutation::Elements> v;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		v[c] = +u[c];
	}
	return v;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator -(const Vector<Type, Size, Permutation> &u) {
	Vector<Type, Size, Permutation> v;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		v[c] = -u[c];
	}
	return v;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator ~(const Vector<Type, Size, Permutation> &u) {
	Vector<Type, Permutation::Elements> v;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		v[c] = ~u[c];
	}
	return v;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator +(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> t(u);
	t += v;
	return t;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator -(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> t(u);
	t -= v;
	return t;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator *(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> t(u);
	t *= v;
	return t;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator /(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> t(u);
	t /= v;
	return t;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator %(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> t(u);
	t %= v;
	return t;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator &(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> t(u);
	t &= v;
	return t;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator |(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> t(u);
	t |= v;
	return t;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<Type, Permutation::Elements> operator ^(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> t(u);
	t ^= v;
	return t;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator <<(const Vector<Type, Size, Permutation> &u, size_t s) {
	Vector<Type, Permutation::Elements> t(u);
	t <<= s;
	return t;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator >>(const Vector<Type, Size, Permutation> &u, size_t s) {
	Vector<Type, Permutation::Elements> t(u);
	t >>= s;
	return t;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> sqrt(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = sqrt(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> sin(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = sin(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> cos(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = cos(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> tan(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = tan(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> asin(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = asin(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> acos(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = acos(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> atan(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = atan(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> sinh(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = sinh(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> cosh(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = cosh(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> tanh(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = tanh(v[c]);
	}
	return ret;
}

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

template <typename Type, size_t Size, typename Permutation>
inline Vector<bool, Permutation::Elements> isnan(const Vector<Type, Size, Permutation> &v) {
	Vector<bool, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = isnan(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<bool, Permutation::Elements> isinf(const Vector<Type, Size, Permutation> &v) {
	Vector<bool, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = isinf(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
Vector<Type, Permutation::Elements> exp2(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = exp2(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
Vector<Type, Permutation::Elements> log2(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = log2(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> sign(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = sign(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> trunc(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = trunc(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> round(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = round(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> roundEven(const Vector<Type, Size, Permutation> &v) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = roundEven(v[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeI, typename PermutationI>
inline Vector<Type, Permutation::Elements> modf(const Vector<Type, Size, Permutation> &v, Vector<Type, SizeI, PermutationI> &i) {
	static_assert(Permutation::Elements == PermutationI::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = modf(v[c], i[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> mix(const Vector<Type, Size, Permutation> &x, const Vector<Type, Size, Permutation> &y, const Type &a) {
	return x * Vector<Type, Permutation::Elements>(Type(1) - a) + y * Vector<Type, Permutation::Elements>(a);
}

template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY, size_t SizeA, typename PermutationA>
inline Vector<Type, Permutation::Elements> mix(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y, const Vector<bool, SizeA, PermutationA> &a) {
	static_assert(Permutation::Elements == PermutationY::Elements && Permutation::Elements == PermutationA::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = mix(x[c], y[c], a[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeX, typename PermutationX>
inline Vector<Type, Permutation::Elements> step(const Vector<Type, Size, Permutation> &edge, const Vector<Type, Size, Permutation> &x) {
	static_assert(Permutation::Elements == PermutationX::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = step(edge[c], x[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> step(const Type &edge, const Vector<Type, Size, Permutation> &x) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = step(edge, x[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t Size1, typename Permutation1, size_t SizeA, typename PermutationA>
inline Vector<Type, Permutation::Elements> smoothstep(const Vector<Type, Size, Permutation> &edge0, const Vector<Type, Size, Permutation> &edge1, const Vector<Type, Size, Permutation> &x) {
	static_assert(Permutation::Elements == Permutation1::Elements && Permutation::Elements == PermutationA::Elements, "Vectors differ in number of elements");
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = smoothstep(edge0[c], edge1[c], x[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> smoothstep(const Type &edge0, const Type &edge1, const Vector<Type, Size, Permutation> &x) {
	Vector<Type, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = smoothstep(edge0, edge1, x[c]);
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeI, typename PermutationI, size_t SizeN, typename PermutationN>
inline Vector<Type, Permutation::Elements> faceforward(const Vector<Type, Size, Permutation> &N, const Vector<Type, Size, Permutation> &I, const Vector<Type, Size, Permutation> &Nref) {
	if (dot(Nref, I) < Type(0)) {
		return N;
	} else {
		return -N;
	}
}

template <typename Type, size_t Size, typename Permutation>
inline Type length(const Vector<Type, Size, Permutation> &x) {
	return sqrt(dot(x, x));
}

template <typename Type, size_t Size, typename Permutation, size_t Size1, typename Permutation1>
inline Type distance(const Vector<Type, Size, Permutation> &p0, const Vector<Type, Size1, Permutation1> &p1) {
	return length(p0 - p1);
}

template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY>
inline Type dot(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y) {
	Vector<Type, Permutation::Elements> temp = x * y;
	Type ret(0);
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret += temp[c];
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeY, typename PermutationY>
inline Vector<Type, 3> cross(const Vector<Type, Size, Permutation> &x, const Vector<Type, SizeY, PermutationY> &y) {
	static_assert(Permutation::Elements == PermutationY::Elements, "Vectors differ in number of elements");
	return Vector<Type, 3>(
		x[1] * y[2] - y[1] * x[2],
		x[2] * y[0] - y[2] * x[0],
		x[0] * y[1] - y[0] * x[1]
	);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> normalize(const Vector<Type, Size, Permutation> &x) {
	return x * Vector<Type, Permutation::Elements>(inversesqrt(dot(x, x)));
}

template <typename Type, size_t Size, typename Permutation, size_t SizeN, typename PermutationN>
inline Vector<Type, Permutation::Elements> reflect(const Vector<Type, Size, Permutation> &I, const Vector<Type, SizeN, PermutationN> &N) {
	return I - Vector<Type, Permutation::Elements>(2) * dot(N, I) * N;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeN, typename PermutationN>
inline Vector<Type, Permutation::Elements> refract(const Vector<Type, Size, Permutation> &I, const Vector<Type, SizeN, PermutationN> &N, const Type &eta) {
	Type k = 1 - eta * eta * (1 - dot(N, I) * dot(N, I));
	if (k < Type(0)) {
		return Vector<Type, Permutation::Elements>();
	} else {
		return eta * I - (eta * dot(N, I) + sqrt(k)) * N;
	}
}

template <size_t Size, typename Permutation>
inline bool any(const Vector<bool, Size, Permutation> &v) {
	bool ret = false;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret = ret || v[c];
	}
	return ret;
}

template <size_t Size, typename Permutation>
inline bool all(const Vector<bool, Size, Permutation> &v) {
	bool ret = true;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret = ret && v[c];
	}
	return ret;
}

template <size_t Size, typename Permutation>
inline Vector<bool, Permutation::Elements> operator !(const Vector<bool, Size, Permutation> &v) {
	Vector <bool, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = !v[c];
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation>
inline std::ostream &operator <<(std::ostream &stream, const Vector<Type, Size, Permutation> &v) {
	stream << '(';
	bool first = true;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		if (first) {
			first = false;
		} else {
			stream << ',';
		}
		stream << v[c];
	}
	stream << ')';
	return stream;
}

template <typename Type, size_t Size, typename Permutation>
inline std::istream &operator >>(std::istream &stream, const Vector<Type, Size, Permutation> &v) {
	for (size_t c = 0; c < Permutation::Elements; c++) {
		stream >> std::skipws >> v[c];
	}
	return stream;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> equal(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<bool, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = u[c] == v[c];
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> notEqual(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<bool, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = u[c] != v[c];
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> lessThan(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<bool, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = u[c] < v[c];
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> greaterThan(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<bool, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = u[c] > v[c];
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> lessThanEqual(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<bool, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = u[c] <= v[c];
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline Vector<bool, Permutation::Elements> greaterThanEqual(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	static_assert(Permutation::Elements == PermutationV::Elements, "Vectors differ in number of elements");
	Vector<bool, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = u[c] >= v[c];
	}
	return ret;
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline bool operator ==(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	return all(equal(u, v));
}

template <typename Type, size_t Size, typename Permutation, size_t SizeV, typename PermutationV>
inline bool operator !=(const Vector<Type, Size, Permutation> &u, const Vector<Type, SizeV, PermutationV> &v) {
	return any(notEqual(u, v));
}

template <size_t Size, typename Permutation>
inline Vector<int, Permutation::Elements> floatBitsToInt(const Vector<float, Size, Permutation> &v) {
	Vector<int, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = floatBitsToInt(v[c]);
	}
	return ret;
}

template <size_t Size, typename Permutation>
inline Vector<unsigned int, Permutation::Elements> floatBitsToUInt(const Vector<float, Size, Permutation> &v) {
	Vector<unsigned int, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = floatBitsToUInt(v[c]);
	}
	return ret;
}

template <size_t Size, typename Permutation>
inline Vector<float, Permutation::Elements> intBitsToFloat(const Vector<int, Size, Permutation> &v) {
	Vector<float, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = intBitsToFloat(v[c]);
	}
	return ret;
}

template <size_t Size, typename Permutation>
inline Vector<float, Permutation::Elements> uintBitsToFloat(const Vector<unsigned int, Size, Permutation> &v) {
	Vector<float, Permutation::Elements> ret;
	for (size_t c = 0; c < Permutation::Elements; c++) {
		ret[c] = uintBitsToFloat(v[c]);
	}
	return ret;
}

template <size_t Size, typename Permutation>
inline unsigned int packSnorm2x16(const Vector<float, Size, Permutation> &v) {
	static_assert(Permutation::Elements == 2, "Vector must have exactly two elements");
	Vector<int, 2> u = round(clamp(v, -1.0f, 1.0f) * 32767.0f);
	return (u[0] & 0xffff) | (u[1] << 16);
}

inline Vector<float, 2> unpackSnorm2x16(unsigned int p) {
	// Cast to short to sign extend
	return clamp(Vector<float, 2>(short(p), short(p >> 16)) / 32767.0f, -1.0f, 1.0f);
}

template <size_t Size, typename Permutation>
inline unsigned int packUnorm2x16(const Vector<float, Size, Permutation> &v) {
	static_assert(Permutation::Elements == 2, "Vector must have exactly two elements");
	Vector<unsigned int, 2> u = round(clamp(v, 0.0f, 1.0f) * 65535.0f);
	return u[0] | (u[1] << 16);
}

inline Vector<float, 2> unpackUnorm2x16(unsigned int p) {
	return Vector<float, 2>(p & 0xffff, p >> 16) / 65535.0f;
}

template <size_t Size, typename Permutation>
inline unsigned int packHalf2x16(const Vector<float, Size, Permutation> &v) {
	static_assert(Permutation::Elements == 2, "Vector must have exactly two elements");
	Vector<unsigned int, 2> bits = floatBitsToUInt(v);
	//          | 31   | 30   | 29   | 28   | 27   | 26   | 25   | 24   | 23   | 22   | ..   | 00   |
	// binary32 | sign | exponent (bias -127)                                  | value              |
	//          | 15   |                      14   | 13   | 12   | 11   | 10   | 09   | ..   | 00   |
	// binary16 | sign |                      exponent (bias -15)              | value              |
	Vector<unsigned int, 2> sign = (bits >> 16) & 0xc000u;
	Vector<unsigned int, 2> expvalue = (bits >> 13) & 0x3fffu;
	bits = sign | expvalue;
	return bits[0] | (bits[1] << 16);
}

inline Vector<float, 2> unpackHalf2x16(unsigned int v) {
	Vector<unsigned int, 2> bits(v & 0xffffu, v >> 16);
	// Unpacking is much more complicated than packing.
	// We need to extend the exponent, which involves special handling for the special numbers.
	Vector<unsigned int, 2> sign = bits & 0x8000u;
	Vector<unsigned int, 2> biasexp = bits & 0x7c00u;
	Vector<unsigned int, 2> value = bits & 0x03ffu;
	// zero or subnormal
	Vector<bool, 2> ezero = equal(biasexp, Vector<unsigned int, 2>(0x0000, 0x0000));
	// infinity or Nan
	Vector<bool, 2> emax = equal(biasexp, Vector<unsigned int, 2>(0x7c00, 0x7c00));
	// anything else
	Vector<bool, 2> eother = !(ezero | emax);
	// bit flipping to avoid multiply and add:
	// calculate the two's complement of the booleans (false => 0x00000000, true => 0xffffffff)
	// mask the appropriate bits
	// combine
	// 0x3fc00 = (0b11111111) << 10
	// 0x1c000 = (127 - 15) << 10
	Vector<unsigned int, 2> exp = (-Vector<unsigned int, 2>(emax) & 0x3fc00u) | (-Vector<unsigned int, 2>(eother) & (biasexp + 0x1c000u));
	bits = (sign << 16) | ((value | exp) << 13);
	return uintBitsToFloat(bits);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator +=(Vector<Type, Size, Permutation> &self, const Type &v) {
	return self += Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator -=(Vector<Type, Size, Permutation> &self, const Type &v) {
	return self -= Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator *=(Vector<Type, Size, Permutation> &self, const Type &v) {
	return self *= Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator /=(Vector<Type, Size, Permutation> &self, const Type &v) {
	return self /= Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator %=(Vector<Type, Size, Permutation> &self, const Type &v) {
	return self %= Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator &=(Vector<Type, Size, Permutation> &self, const Type &v) {
	return self &= Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator |=(Vector<Type, Size, Permutation> &self, const Type &v) {
	return self |= Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Size, Permutation> &operator ^=(Vector<Type, Size, Permutation> &self, const Type &v) {
	return self ^= Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator +(const Vector<Type, Size, Permutation> &u, const Type &v) {
	return u + Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator +(const Type &u, const Vector<Type, Size, Permutation> &v) {
	return Vector<Type, Permutation::Elements>(u) + v;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator -(const Vector<Type, Size, Permutation> &u, const Type &v) {
	return u - Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator -(const Type &u, const Vector<Type, Size, Permutation> &v) {
	return Vector<Type, Permutation::Elements>(u) - v;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator *(const Vector<Type, Size, Permutation> &u, const Type &v) {
	return u * Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator *(const Type &u, const Vector<Type, Size, Permutation> &v) {
	return Vector<Type, Permutation::Elements>(u) * v;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator /(const Vector<Type, Size, Permutation> &u, const Type &v) {
	return u / Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator /(const Type &u, const Vector<Type, Size, Permutation> &v) {
	return Vector<Type, Permutation::Elements>(u) / v;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator %(const Vector<Type, Size, Permutation> &u, const Type &v) {
	return u % Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator %(const Type &u, const Vector<Type, Size, Permutation> &v) {
	return Vector<Type, Permutation::Elements>(u) % v;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator &(const Vector<Type, Size, Permutation> &u, const Type &v) {
	return u & Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator &(const Type &u, const Vector<Type, Size, Permutation> &v) {
	return Vector<Type, Permutation::Elements>(u) & v;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator |(const Vector<Type, Size, Permutation> &u, const Type &v) {
	return u | Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator |(const Type &u, const Vector<Type, Size, Permutation> &v) {
	return Vector<Type, Permutation::Elements>(u) | v;
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator ^(const Vector<Type, Size, Permutation> &u, const Type &v) {
	return u ^ Vector<Type, Permutation::Elements>(v);
}

template <typename Type, size_t Size, typename Permutation>
inline Vector<Type, Permutation::Elements> operator ^(const Type &u, const Vector<Type, Size, Permutation> &v) {
	return Vector<Type, Permutation::Elements>(u) ^ v;
}


}

#endif //GLAM_VECTOR_H
