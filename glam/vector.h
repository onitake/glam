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
#include <string>
#include <vector>
#include <list>
#include <sstream>
#include <glam/math.h>
#include <glam/config.h>
#include <glam/exception.h>

namespace glam {

// Declaration

template <class T, unsigned int C>
class Vector {
private:
	T _v[C];

	// Overloadable template functions for use by the generator macro
	inline T permutation(unsigned int x) const;
	inline Vector<T, 2> permutation(unsigned int x, unsigned int y) const;
	inline Vector<T, 3> permutation(unsigned int x, unsigned int y, unsigned int z) const;
	inline Vector<T, 4> permutation(unsigned int x, unsigned int y, unsigned int z, unsigned int w) const;
	
	// Supporting initializers to avoid trouble from delegate constructors
	template <unsigned int X>
	inline void initialize();
	template <unsigned int X, typename ForwardIterator>
	inline void initialize(std::pair<ForwardIterator, ForwardIterator> iterator);
	template <unsigned int X>
	inline void initialize(const T *v);
	template <unsigned int X, class U>
	inline void initialize(const Vector<U, C> &v);
	template <unsigned int X>
	inline void initialize(const T &s);
	template <unsigned int X, typename... Args, unsigned int P>
	inline void initialize(const Vector<T, P> &v, Args... args);
	template <unsigned int X, typename... Args>
	inline void initialize(const T &v, Args... args);
	
public:
	// Catch-all variadic constructor
	// The following constructor calls are supported:
	// Vector(): Creates an empty vector
	//  All elements are 0.
	// Vector(const T &s): Extend scalar to vector
	//  All elements are initialized to the same value.
	// Vector<ForwardIterator>(std::pair<ForwardIterator, ForwardIterator> iterator): Populate the vector from an iterator pair
	//  Elements are filled up in order of occurrence, from first to last element (ex. x, y, z, w).
	//  std::pair is used to avoid call ambiguity problems. Call std::make_pair(begin, end) to constuct your iterator pair.
	// Vector<AnyScalarOrVector...>(const AnyScalarOrVector &arg0, ...): Populate the vector from a list of vectors and/or scalars
	//  Elements are filled up in order of occurrence.
	// Vector(const T *m): Populate the vector from a constant array
	//  Elements are filled up in order of occurrence.
	//  Note that this overload does not support range checking. Make sure that the input array has the correct size or use the ForwardIterator overload.
	// Vector<OtherType>(const Vector<OtherType> &arg0): Convert a vector of different element type
	//  The lengths of the vectors need to match.
	template <typename... Args>
	inline Vector(Args... args);

	// Assignment operator
	inline Vector<T, C> &operator =(const Vector<T, C> &v);
	
	// Member access operator
	inline T &operator [](unsigned int index);
	// Const member access operator
	inline const T &operator [](unsigned int index) const;
	// Direct access to internal pointer
	inline const T *internal() const;

	// Permutation functions
	// These are in the form of function calls, not member access like in GLSL
	// Eg.: v.xxyy() instead of v.xxyy
	// All the 340 1-, 2-, 3- and 4-member permutations are generated on the fly.
	#define GLAM_PERM_MAKE1_1(type, length, name, member0) inline T name() const { static_assert(member0 < length, "Invalid member index for this vector type"); return permutation(member0); }
	#define GLAM_PERM_MAKE2_2(type, length, name, member0, member1) inline Vector<type, 2> name() const { static_assert(member0 < length && member1 < length, "Invalid member index for this vector type"); return permutation(member0, member1); }
	#define GLAM_PERM_MAKE3_3(type, length, name, member0, member1, member2) inline Vector<type, 3> name() const { static_assert(member0 < length && member1 < length && member2 < length, "Invalid member index for this vector type"); return permutation(member0, member1, member2); }
	#define GLAM_PERM_MAKE4_4(type, length, name, member0, member1, member2, member3) inline Vector<type, 4> name() const { static_assert(member0 < length && member1 < length && member2 < length && member3 < length, "Invalid member index for this vector type"); return permutation(member0, member1, member2, member3); }
	#define GLAM_PERM_MAKE4_3(type, length, prefix, member0, member1, member2) \
		GLAM_PERM_MAKE4_4(type, length, prefix##x, member0, member1, member2, 0) \
		GLAM_PERM_MAKE4_4(type, length, prefix##y, member0, member1, member2, 1) \
		GLAM_PERM_MAKE4_4(type, length, prefix##z, member0, member1, member2, 2) \
		GLAM_PERM_MAKE4_4(type, length, prefix##w, member0, member1, member2, 3) \
		GLAM_PERM_MAKE3_3(type, length, prefix, member0, member1, member2)
	#define GLAM_PERM_MAKE4_2(type, length, prefix, member0, member1) \
		GLAM_PERM_MAKE4_3(type, length, prefix##x, member0, member1, 0) \
		GLAM_PERM_MAKE4_3(type, length, prefix##y, member0, member1, 1) \
		GLAM_PERM_MAKE4_3(type, length, prefix##z, member0, member1, 2) \
		GLAM_PERM_MAKE4_3(type, length, prefix##w, member0, member1, 3) \
		GLAM_PERM_MAKE2_2(type, length, prefix, member0, member1)
	#define GLAM_PERM_MAKE4_1(type, length, prefix, member0) \
		GLAM_PERM_MAKE4_2(type, length, prefix##x, member0, 0) \
		GLAM_PERM_MAKE4_2(type, length, prefix##y, member0, 1) \
		GLAM_PERM_MAKE4_2(type, length, prefix##z, member0, 2) \
		GLAM_PERM_MAKE4_2(type, length, prefix##w, member0, 3) \
		GLAM_PERM_MAKE1_1(type, length, prefix, member0)
	#define GLAM_PERM_MAKE4_0(type, length) \
		GLAM_PERM_MAKE4_1(type, length, x, 0) \
		GLAM_PERM_MAKE4_1(type, length, y, 1) \
		GLAM_PERM_MAKE4_1(type, length, z, 2) \
		GLAM_PERM_MAKE4_1(type, length, w, 3)
	GLAM_PERM_MAKE4_0(T, C)
	#undef GLAM_PERM_MAKE1_1
	#undef GLAM_PERM_MAKE2_2
	#undef GLAM_PERM_MAKE3_3
	#undef GLAM_PERM_MAKE4_4
	#undef GLAM_PERM_MAKE4_3
	#undef GLAM_PERM_MAKE4_2
	#undef GLAM_PERM_MAKE4_1
	#undef GLAM_PERM_MAKE4_0
};

// Component-wise pre-increment operator
template <class T, unsigned int C>
inline Vector<T, C> &operator ++(Vector<T, C> &self);
// Component-wise post-increment operator
template <class T, unsigned int C>
inline Vector<T, C> operator ++(Vector<T, C> &self, int);
// Component-wise pre-decrement operator
template <class T, unsigned int C>
inline Vector<T, C> &operator --(Vector<T, C> &self);
// Component-wise post-decrement operator
template <class T, unsigned int C>
inline Vector<T, C> operator --(Vector<T, C> &self, int);

// Component-wise self-addition operator
template <class T, unsigned int C>
inline Vector<T, C> &operator +=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise self-subtraction operator
template <class T, unsigned int C>
inline Vector<T, C> &operator -=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise self-multiplication operator
template <class T, unsigned int C>
inline Vector<T, C> &operator *=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise self-division operator
template <class T, unsigned int C>
inline Vector<T, C> &operator /=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise self-modulus operator
template <class T, unsigned int C>
inline Vector<T, C> &operator %=(Vector<T, C> &self, const Vector<T, C> &v);

// Component-wise bitwise and operator
template <class T, unsigned int C>
inline Vector<T, C> &operator &=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise bitwise or operator
template <class T, unsigned int C>
inline Vector<T, C> &operator |=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise bitwise xor operator
template <class T, unsigned int C>
inline Vector<T, C> &operator ^=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise bit shift left operator
template <class T, unsigned int C>
inline Vector<T, C> &operator <<=(Vector<T, C> &self, unsigned int s);
// Component-wise bit shift right operator
template <class T, unsigned int C>
inline Vector<T, C> &operator >>=(Vector<T, C> &self, unsigned int s);

// Component-wise unary plus operator (no-op if T follows standard algebra)
template <class T, unsigned int C>
inline Vector<T, C> operator +(const Vector<T, C> &u);
// Component-wise negation operator (negation)
template <class T, unsigned int C>
inline Vector<T, C> operator -(const Vector<T, C> &u);
// Component-wise ones complement operator
template <class T, unsigned int C>
inline Vector<T, C> operator ~(const Vector<T, C> &u);

// Component-wise addition operator
template <class T, unsigned int C>
inline Vector<T, C> operator +(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise subtraction operator
template <class T, unsigned int C>
inline Vector<T, C> operator -(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise multiplication operator
template <class T, unsigned int C>
inline Vector<T, C> operator *(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise division operator
template <class T, unsigned int C>
inline Vector<T, C> operator /(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise modulus operator
template <class T, unsigned int C>
inline Vector<T, C> operator %(const Vector<T, C> &u, const Vector<T, C> &v);

// Component-wise bitwise and operator
template <class T, unsigned int C>
inline Vector<T, C> operator &(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise bitwise or operator
template <class T, unsigned int C>
inline Vector<T, C> operator |(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise bitwise xor operator
template <class T, unsigned int C>
inline Vector<T, C> operator ^(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise bit shift left operator
template <class T, unsigned int C>
inline Vector<T, C> operator <<(const Vector<T, C> &u, unsigned int s);
// Component-wise bit shift right operator
template <class T, unsigned int C>
inline Vector<T, C> operator >>(const Vector<T, C> &u, unsigned int s);

// Component-wise equality comparison operator (exact compare)
template <class T, unsigned int C>
inline Vector<bool, C> equal(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise inequality comparison operator (exact compare)
template <class T, unsigned int C>
inline Vector<bool, C> notEqual(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise less-than comparison operator (exact compare)
template <class T, unsigned int C>
inline Vector<bool, C> lessThan(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise greater-than comparison operator (exact compare)
template <class T, unsigned int C>
inline Vector<bool, C> greaterThan(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise less-or-equal comparison operator (exact compare)
template <class T, unsigned int C>
inline Vector<bool, C> lessThanEqual(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise greater-or-equal comparison operator (exact compare)
template <class T, unsigned int C>
inline Vector<bool, C> greaterThanEqual(const Vector<T, C> &u, const Vector<T, C> &v);

// Vector comparison operator (exact compare), correspondes to all(equal(u, v))
template <class T, unsigned int C>
inline bool operator ==(const Vector<T, C> &u, const Vector<T, C> &v);
// Vector comparison operator (exact compare), correspondes to any(notEqual(u, v))
template <class T, unsigned int C>
inline bool operator !=(const Vector<T, C> &u, const Vector<T, C> &v);

// Component-wise sine
template <class T, unsigned int C>
inline Vector<T, C> sin(const Vector<T, C> &v);
// Component-wise cosine
template <class T, unsigned int C>
inline Vector<T, C> cos(const Vector<T, C> &v);
// Component-wise tangent
template <class T, unsigned int C>
inline Vector<T, C> tan(const Vector<T, C> &v);
// Component-wise arc sine
template <class T, unsigned int C>
inline Vector<T, C> asin(const Vector<T, C> &v);
// Component-wise arc cosine
template <class T, unsigned int C>
inline Vector<T, C> acos(const Vector<T, C> &v);
// Component-wise arc tangent
template <class T, unsigned int C>
inline Vector<T, C> atan(const Vector<T, C> &v);
// Component-wise sinus hyperbolicus
template <class T, unsigned int C>
inline Vector<T, C> sinh(const Vector<T, C> &v);
// Component-wise cosinus hyperbolicus
template <class T, unsigned int C>
inline Vector<T, C> cosh(const Vector<T, C> &v);
// Component-wise tangens hyperbolicus
template <class T, unsigned int C>
inline Vector<T, C> tanh(const Vector<T, C> &v);

// Component-wise power
template <class T, unsigned int C>
inline Vector<T, C> pow(const Vector<T, C> &x, const Vector<T, C> &y);
// Component-wise exponential (e^v)
template <class T, unsigned int C>
inline Vector<T, C> exp(const Vector<T, C> &v);
// Component-wise natural logarithm
template <class T, unsigned int C>
inline Vector<T, C> log(const Vector<T, C> &v);
// Component-wise power of 2
template <class T, unsigned int C>
inline Vector<T, C> exp2(const Vector<T, C> &v);
// Component-wise base 2 logarithm
template <class T, unsigned int C>
inline Vector<T, C> log2(const Vector<T, C> &v);
// Component-wise square root
template <class T, unsigned int C>
inline Vector<T, C> sqrt(const Vector<T, C> &v);

// Component-wise absolute value
template <class T, unsigned int C>
inline Vector<T, C> abs(const Vector<T, C> &v);
// Component-wise sign
template <class T, unsigned int C>
inline Vector<T, C> sign(const Vector<T, C> &v);
// Component-wise round towards negative infinity
template <class T, unsigned int C>
inline Vector<T, C> floor(const Vector<T, C> &v);
// Component-wise truncate
template <class T, unsigned int C>
inline Vector<T, C> trunc(const Vector<T, C> &v);
// Component-wise round to the nearest integer
template <class T, unsigned int C>
inline Vector<T, C> round(const Vector<T, C> &v);
// Component-wise round to the nearest even integer
template <class T, unsigned int C>
inline Vector<T, C> roundEven(const Vector<T, C> &v);
// Component-wise round towards positive infinity
template <class T, unsigned int C>
inline Vector<T, C> ceil(const Vector<T, C> &v);
// Component-wise fractional part
template <class T, unsigned int C>
inline Vector<T, C> fract(const Vector<T, C> &v);
// Component-wise fractional part separation
template <class T, unsigned int C>
inline Vector<T, C> modf(const Vector<T, C> &v, Vector<T, C> &i);
// Component-wise minimum
template <class T, unsigned int C>
inline Vector<T, C> min(const Vector<T, C> &x, const Vector<T, C> &y);
// Component-wise maximum
template <class T, unsigned int C>
inline Vector<T, C> max(const Vector<T, C> &x, const Vector<T, C> &y);
// Saturation with single min/max values (applied to all vector components)
template <class T, unsigned int C>
inline Vector<T, C> clamp(const Vector<T, C> &x, const T &minVal, const T &maxVal);
// Component-wise linear blend (x * (1 - a) + y * a)
template <class T, unsigned int C>
inline Vector<T, C> mix(const Vector<T, C> &x, const Vector<T, C> &y, const T &a);
// Return x for each component of a that is false, y if true
template <class T, unsigned int C>
inline Vector<T, C> mix(const Vector<T, C> &x, const Vector<T, C> &y, const Vector<bool, C> &a);
// Component-wise step transition
template <class T, unsigned int C>
inline Vector<T, C> step(const Vector<T, C> &edge, const Vector<T, C> &x);
// Component-wise step transition with common edge argument
template <class T, unsigned int C>
inline Vector<T, C> step(const T &edge, const Vector<T, C> &x);
// Component-wise smooth step transition
template <class T, unsigned int C>
inline Vector<T, C> smoothstep(const Vector<T, C> &edge0, const Vector<T, C> &edge1, const Vector<T, C> &x);
// Component-wise smooth step transition with common edge arguments
template <class T, unsigned int C>
inline Vector<T, C> smoothstep(const T &edge0, const T &edge1, const Vector<T, C> &x);
// Component-wise test for NaN
template <class T, unsigned int C>
inline Vector<bool, C> isnan(const Vector<T, C> &v);
// Component-wise test for infinity
template <class T, unsigned int C>
inline Vector<bool, C> isinf(const Vector<T, C> &v);

// Vector length (sqrt(dot(a, a)))
template <class T, unsigned int C>
inline T length(const Vector<T, C> &x);
// Distance between two vectors (length(b - a))
template <class T, unsigned int C>
inline T distance(const Vector<T, C> &p0, const Vector<T, C> &p1);
// Dot product of two vectors (a[0] * b[0] + a[1] * b[1] + ...)
template <class T, unsigned int C>
inline T dot(const Vector<T, C> &x, const Vector<T, C> &y);
// Cross product of two 3-component vectors
template <class T>
inline Vector<T, 3> cross(const Vector<T, 3> &x, const Vector<T, 3> &y);
// Unit length vector (a / length(a))
// always_inline encouraged to keep gcc from generating crap code
template <class T, unsigned int C>
inline Vector<T, C> normalize(const Vector<T, C> &x) __attribute__((always_inline));
// Check if I and Nref are facing in opposite directions and return N if yes, -N otherwise (dot(Nref, I) < 0 ? N : -N)
template <class T>
inline T faceforward(const T &N, const T &I, const T &Nref);
// Reflect I in respect to surface normal N
template <class T, unsigned int C>
inline Vector<T, C> reflect(const Vector<T, C> &I, const Vector<T, C> &N);
// Refract I on the surface defined by normal N using index of refraction ration eta
template <class T, unsigned int C>
inline Vector<T, C> refract(const Vector<T, C> &I, const Vector<T, C> &N, const T &eta);

// True if any component of x is true (or evaluates to true)
template <unsigned int C>
inline bool any(const Vector<bool, C> &v);
// True if all components of x are true (or evaluates to true)
template <unsigned int C>
inline bool all(const Vector<bool, C> &v);
// Component-wise negation of a boolean vector, this is an overloaded ! operator because
// 'not' is a reserved word in C++ with the same semantics as !
template <unsigned int C>
inline Vector<bool, C> operator !(const Vector<bool, C> &v);

// Formatted output operator (vectors will be presented in the form "(0.1,1.50)")
template <class T, unsigned int C>
inline std::ostream &operator <<(std::ostream &stream, const Vector<T, C> &v);
// Read vector components from stream, in text form, separated by whitespace
template <class T, unsigned int C>
inline std::istream &operator >>(std::istream &stream, const Vector<T, C> &v);

// Convert a floating point vector into an integer vector of equal size and exactly the same bit pattern
template <unsigned int C>
inline Vector<int, C> floatBitsToInt(const Vector<float, C> &v);
template <unsigned int C>
inline Vector<unsigned int, C> floatBitsToUInt(const Vector<float, C> &v);
// Convert an integer vector into a floating point vector of equal size and exactly the same bit pattern
template <unsigned int C>
inline Vector<float, C> intBitsToFloat(const Vector<int, C> &v);
template <unsigned int C>
inline Vector<float, C> uintBitsToFloat(const Vector<unsigned int, C> &v);

// Convert two normalized floats into 16bit fixed point values, then pack them into a 32bit integer
// The first vector component will be at the lower 16 bits of the result, the second at the upper part.
// Conversion from float to fixed is equivalent to: round(clamp(c, -1, +1) * 32767.0)
inline unsigned int packSnorm2x16(const Vector<float, 2> &v);
// Convert two packed 16bit fixed point values into a vector of normalized floats
// The first vector component is taken from the lower 16 bits, the second from the upper part.
// Conversion from fixed to float is equivalent to: clamp(f / 32767.0, -1, +1)
inline Vector<float, 2> unpackSnorm2x16(unsigned int p);
// Convert two normalized floats into 16bit fixed point values, then pack them into a 32bit integer
// The first vector component will be at the lower 16 bits of the result, the second at the upper part.
// Conversion from float to fixed is equivalent to: round(clamp(c, 0, +1) * 65535.0)
inline unsigned int packUnorm2x16(const Vector<float, 2> &v);
// Convert two packed 16bit fixed point values into a vector of normalized floats
// The first vector component is taken from the lower 16 bits, the second from the upper part.
// Conversion from fixed to float is equivalent to: f / 65535.0
inline Vector<float, 2> unpackUnorm2x16(unsigned int p);


// GLSL types

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


// Implementation

template <class T, unsigned int C>
inline T Vector<T, C>::permutation(unsigned int x) const {
	return _v[x];
}

template <class T, unsigned int C>
inline Vector<T, 2> Vector<T, C>::permutation(unsigned int x, unsigned int y) const {
	return Vector<T, 2>(_v[x], _v[y]);
}

template <class T, unsigned int C>
inline Vector<T, 3> Vector<T, C>::permutation(unsigned int x, unsigned int y, unsigned int z) const {
	return Vector<T, 3>(_v[x], _v[y], _v[z]);
}

template <class T, unsigned int C>
inline Vector<T, 4> Vector<T, C>::permutation(unsigned int x, unsigned int y, unsigned int z, unsigned int w) const {
	return Vector<T, 4>(_v[x], _v[y], _v[z], _v[w]);
}

template <class T, unsigned int C>
template <unsigned int X>
inline void Vector<T, C>::initialize() {
	static_assert(X == 0 || X == C, "Not enough initializers");
	if (X == 0) {
		for (unsigned int p = 0; p < C; p++) {
			(*this)[p] = static_cast<T>(0);
		}
	}
}

template <class T, unsigned int C>
template <unsigned int X, typename... Args>
inline void Vector<T, C>::initialize(const T &s, Args... args) {
	// The condition is < and not <= because there is a dedicated terminal initializer for a single scalar
	static_assert(X + 1 < C, "Too many initializers");
	(*this)[X] = s;
	initialize<X + 1>(args...);
}

template <class T, unsigned int C>
template <unsigned int X>
inline void Vector<T, C>::initialize(const T &s) {
	static_assert(X + 1 <= C, "Too many initializers");
	static_assert(X == 0 || X + 1 == C, "Not enough initializers");
	if (X == 0) {
		for (unsigned int p = 0; p < C; p++) {
			(*this)[p] = s;
		}
	} else {
		(*this)[X] = s;
	}
}

template <class T, unsigned int C>
template <unsigned int X, typename... Args, unsigned int P>
inline void Vector<T, C>::initialize(const Vector<T, P> &v, Args... args) {
	static_assert(X + P <= C, "Too many initializers");
	for (unsigned int p = 0; p < P; p++) {
		(*this)[X + p] = v[p];
	}
	initialize<X + P>(args...);
}

template <class T, unsigned int C>
template <unsigned int X, typename ForwardIterator>
inline void Vector<T, C>::initialize(std::pair<ForwardIterator, ForwardIterator> iterator) {
	static_assert(X == 0, "No other arguments are allowed when initializing a vector from an interator");
	unsigned int p = X;
	for (ForwardIterator it = iterator.first; p < C && it != iterator.second; p++, it++) {
		(*this)[p] = *it;
	}
}

template <class T, unsigned int C>
template <unsigned int X>
inline void Vector<T, C>::initialize(const T *v) {
	static_assert(X == 0, "No other arguments are allowed when initializing a vector from an array");
	for (unsigned int p = 0; p < C; p++) {
		(*this)[X + p] = v[p];
	}
}

template <class T, unsigned int C>
template <unsigned int X, class U>
inline void Vector<T, C>::initialize(const Vector<U, C> &v) {
	static_assert(X == 0, "No other arguments are allowed when initializing a vector from a vector of different element type");
	for (unsigned int p = 0; p < C; p++) {
		(*this)[X + p] = static_cast<T>(v[p]);
	}
}

template <class T, unsigned int C>
template <typename... Args>
inline Vector<T, C>::Vector(Args... args) {
	initialize<0>(args...);
}

template <class T, unsigned int C>
inline T &Vector<T, C>::operator [](unsigned int index) {
#ifdef GLAM_RANGE_CHECKS
	if (index < C) {
		return _v[index];
	}
	throw DimensionOutOfRangeException("Invalid vector element index", C - 1, index);
#else
	return _v[index];
#endif
}

template <class T, unsigned int C>
inline const T &Vector<T, C>::operator [](unsigned int index) const {
#ifdef GLAM_RANGE_CHECKS
	if (index < C) {
		return _v[index];
	}
	throw DimensionOutOfRangeException("Invalid vector element index", C - 1, index);
#else
	return _v[index];
#endif
}

template <class T, unsigned int C>
inline Vector<T, C> &Vector<T, C>::operator =(const Vector<T, C> &v) {
	for (unsigned int c = 0; c < C; c++) {
		_v[c] = v[c];
	}
	return *this;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator ++(Vector<T, C> &self) {
	for (unsigned int c = 0; c < C; c++) {
		self[c]++;
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> operator ++(Vector<T, C> &self, int) {
	Vector<T, C> temp(self);
	for (unsigned int c = 0; c < C; c++) {
		self[c]++;
	}
	return temp;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator --(Vector<T, C> &self) {
	for (unsigned int c = 0; c < C; c++) {
		self[c]--;
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> operator --(Vector<T, C> &self, int) {
	Vector<T, C> temp(self);
	for (unsigned int c = 0; c < C; c++) {
		self[c]--;
	}
	return temp;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator +=(Vector<T, C> &self, const Vector<T, C> &v) {
	for (unsigned int c = 0; c < C; c++) {
		self[c] += v[c];
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator -=(Vector<T, C> &self, const Vector<T, C> &v) {
	for (unsigned int c = 0; c < C; c++) {
		self[c] -= v[c];
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator *=(Vector<T, C> &self, const Vector<T, C> &v) {
	for (unsigned int c = 0; c < C; c++) {
		self[c] *= v[c];
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator /=(Vector<T, C> &self, const Vector<T, C> &v) {
	for (unsigned int c = 0; c < C; c++) {
		self[c] /= v[c];
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator %=(Vector<T, C> &self, const Vector<T, C> &v) {
	for (unsigned int c = 0; c < C; c++) {
		self[c] %= v[c];
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator &=(Vector<T, C> &self, const Vector<T, C> &v) {
	for (unsigned int c = 0; c < C; c++) {
		self[c] &= v[c];
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator |=(Vector<T, C> &self, const Vector<T, C> &v) {
	for (unsigned int c = 0; c < C; c++) {
		self[c] |= v[c];
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator ^=(Vector<T, C> &self, const Vector<T, C> &v) {
	for (unsigned int c = 0; c < C; c++) {
		self[c] ^= v[c];
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator <<=(Vector<T, C> &self, unsigned int s) {
	for (unsigned int c = 0; c < C; c++) {
		self[c] <<= s;
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> &operator >>=(Vector<T, C> &self, unsigned int s) {
	for (unsigned int c = 0; c < C; c++) {
		self[c] >>= s;
	}
	return self;
}

template <class T, unsigned int C>
inline Vector<T, C> operator +(const Vector<T, C> &u) {
	Vector<T, C> v;
	for (unsigned int c = 0; c < C; c++) {
		v[c] = +u[c];
	}
	return v;
}

template <class T, unsigned int C>
inline Vector<T, C> operator -(const Vector<T, C> &u) {
	Vector<T, C> v;
	for (unsigned int c = 0; c < C; c++) {
		v[c] = -u[c];
	}
	return v;
}

template <class T, unsigned int C>
inline Vector<T, C> operator ~(const Vector<T, C> &u) {
	Vector<T, C> v;
	for (unsigned int c = 0; c < C; c++) {
		v[c] = ~u[c];
	}
	return v;
}

template <class T, unsigned int C>
inline Vector<T, C> operator +(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<T, C> t(u);
	t += v;
	return t;
}

template <class T, unsigned int C>
inline Vector<T, C> operator -(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<T, C> t(u);
	t -= v;
	return t;
}

template <class T, unsigned int C>
inline Vector<T, C> operator *(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<T, C> t(u);
	t *= v;
	return t;
}

template <class T, unsigned int C>
inline Vector<T, C> operator /(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<T, C> t(u);
	t /= v;
	return t;
}

template <class T, unsigned int C>
inline Vector<T, C> operator %(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<T, C> t(u);
	t %= v;
	return t;
}

template <class T, unsigned int C>
inline Vector<T, C> operator &(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<T, C> t(u);
	t &= v;
	return t;
}

template <class T, unsigned int C>
inline Vector<T, C> operator |(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<T, C> t(u);
	t |= v;
	return t;
}

template <class T, unsigned int C>
inline Vector<T, C> operator ^(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<T, C> t(u);
	t ^= v;
	return t;
}

template <class T, unsigned int C>
inline Vector<T, C> operator <<(const Vector<T, C> &u, unsigned int s) {
	Vector<T, C> t(u);
	t <<= s;
	return t;
}

template <class T, unsigned int C>
inline Vector<T, C> operator >>(const Vector<T, C> &u, unsigned int s) {
	Vector<T, C> t(u);
	t >>= s;
	return t;
}

template <class T, unsigned int C>
inline Vector<T, C> sqrt(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = sqrt(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> sin(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = sin(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> cos(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = cos(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> tan(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = tan(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> asin(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = asin(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> acos(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = acos(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> atan(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = atan(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> sinh(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = sinh(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> cosh(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = cosh(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> tanh(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = tanh(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> pow(const Vector<T, C> &x, const Vector<T, C> &y) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = pow(x[c], y[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> exp(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = exp(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> log(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = log(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> abs(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = abs(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> floor(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = floor(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> ceil(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = ceil(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> min(const Vector<T, C> &x, const Vector<T, C> &y) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = min(x[c], y[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> max(const Vector<T, C> &x, const Vector<T, C> &y) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = max(x[c], y[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<bool, C> isnan(const Vector<T, C> &v) {
	Vector<bool, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = isnan(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<bool, C> isinf(const Vector<T, C> &v) {
	Vector<bool, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = isinf(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
Vector<T, C> exp2(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = exp2(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
Vector<T, C> log2(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = log2(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> sign(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = sign(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> trunc(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = trunc(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> round(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = round(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> roundEven(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = roundEven(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> fract(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = fract(v[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> modf(const Vector<T, C> &v, Vector<T, C> &i) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = modf(v[c], i[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> clamp(const Vector<T, C> &x, const T &minVal, const T &maxVal) {
	return min(max(x, Vector<T, C>(minVal)), Vector<T, C>(maxVal));
}

template <class T, unsigned int C>
inline Vector<T, C> mix(const Vector<T, C> &x, const Vector<T, C> &y, const T &a) {
	return x * Vector<T, C>(T(1) - a) + y * Vector<T, C>(a);
}

template <class T, unsigned int C>
inline Vector<T, C> mix(const Vector<T, C> &x, const Vector<T, C> &y, const Vector<bool, C> &a) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = mix(x[c], y[c], a[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> step(const Vector<T, C> &edge, const Vector<T, C> &x) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = step(edge[c], x[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> step(const T &edge, const Vector<T, C> &x) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = step(edge, x[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> smoothstep(const Vector<T, C> &edge0, const Vector<T, C> &edge1, const Vector<T, C> &x) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = smoothstep(edge0[c], edge1[c], x[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> smoothstep(const T &edge0, const T &edge1, const Vector<T, C> &x) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = smoothstep(edge0, edge1, x[c]);
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<T, C> faceforward(const Vector<T, C> &N, const Vector<T, C> &I, const Vector<T, C> &Nref) {
	if (dot(Nref, I) < T(0)) {
		return N;
	} else {
		return -N;
	}
}

template <class T, unsigned int C>
inline T length(const Vector<T, C> &x) {
	return sqrt(dot(x, x));
}

template <class T, unsigned int C>
inline T distance(const Vector<T, C> &p0, const Vector<T, C> &p1) {
	return length(p0 - p1);
}

template <class T, unsigned int C>
inline T dot(const Vector<T, C> &x, const Vector<T, C> &y) {
	if (C == 0) {
		return 0;
	} else {
		Vector<T, C> temp = x * y;
		T ret = temp[0];
		for (unsigned int c = 1; c < C; c++) {
			ret += temp[c];
		}
		return ret;
	}
}

template <class T>
inline Vector<T, 3> cross(const Vector<T, 3> &x, const Vector<T, 3> &y) {
	return Vector<T, 3>(
		x[1] * y[2] - y[1] * x[2],
		x[2] * y[0] - y[2] * x[0],
		x[0] * y[1] - y[0] * x[1]
	);
}

template <class T, unsigned int C>
inline Vector<T, C> normalize(const Vector<T, C> &x) {
	return x * Vector<T, C>(inversesqrt(dot(x, x)));
}

template <class T, unsigned int C>
inline Vector<T, C> reflect(const Vector<T, C> &I, const Vector<T, C> &N) {
	Vector<T, C> temp(dot(N, I) * N);
	return I - Vector<T, C>(2) * dot(N, I) * N;
}

template <class T, unsigned int C>
inline Vector<T, C> refract(const Vector<T, C> &I, const Vector<T, C> &N, const T &eta) {
	T k = 1 - eta * eta * (1 - dot(N, I) * dot(N, I));
	if (k < T(0)) {
		return Vector<T, C>();
	} else {
		return eta * I - (eta * dot(N, I) + sqrt(k)) * N;
	}
}

template <unsigned int C>
inline bool any(const Vector<bool, C> &v) {
	bool ret = false;
	for (unsigned int c = 0; c < C; c++) {
		ret = ret || v[c];
	}
	return ret;
}

template <unsigned int C>
inline bool all(const Vector<bool, C> &v) {
	bool ret = true;
	for (unsigned int c = 0; c < C; c++) {
		ret = ret && v[c];
	}
	return ret;
}

template <unsigned int C>
inline Vector<bool, C> operator !(const Vector<bool, C> &v) {
	Vector <bool, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = !v[c];
	}
	return ret;
}

template <class T, unsigned int C>
inline std::ostream &operator <<(std::ostream &stream, const Vector<T, C> &v) {
	stream << "(";
	for (unsigned int c = 0; c < C - 1; c++) {
		stream << v[c] << ",";
	}
	if (C > 0) {
		stream << v[C - 1];
	}
	stream << ")";
	return stream;
}

template <class T, unsigned int C>
inline std::istream &operator >>(std::istream &stream, const Vector<T, C> &v) {
	for (unsigned int c = 0; c < C; c++) {
		stream >> std::skipws >> v[c];
	}
	return stream;
}

template <class T, unsigned int C>
inline Vector<bool, C> equal(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<bool, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = u[c] == v[c];
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<bool, C> notEqual(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<bool, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = u[c] != v[c];
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<bool, C> lessThan(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<bool, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = u[c] < v[c];
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<bool, C> greaterThan(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<bool, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = u[c] > v[c];
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<bool, C> lessThanEqual(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<bool, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = u[c] <= v[c];
	}
	return ret;
}

template <class T, unsigned int C>
inline Vector<bool, C> greaterThanEqual(const Vector<T, C> &u, const Vector<T, C> &v) {
	Vector<bool, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = u[c] >= v[c];
	}
	return ret;
}

template <class T, unsigned int C>
inline bool operator ==(const Vector<T, C> &u, const Vector<T, C> &v) {
	return all(equal(u, v));
}

template <class T, unsigned int C>
inline bool operator !=(const Vector<T, C> &u, const Vector<T, C> &v) {
	return any(notEqual(u, v));
}

template <class T, unsigned int C>
inline const T *Vector<T, C>::internal() const {
	return _v;
}

template <unsigned int C>
inline Vector<int, C> floatBitsToInt(const Vector<float, C> &v) {
	Vector<int, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = floatBitsToInt(v[c]);
	}
	return ret;
}

template <unsigned int C>
inline Vector<unsigned int, C> floatBitsToUInt(const Vector<float, C> &v) {
	Vector<unsigned int, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = floatBitsToUInt(v[c]);
	}
	return ret;
}

template <unsigned int C>
inline Vector<float, C> intBitsToFloat(const Vector<int, C> &v) {
	Vector<float, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = intBitsToFloat(v[c]);
	}
	return ret;
}

template <unsigned int C>
inline Vector<float, C> uintBitsToFloat(const Vector<unsigned int, C> &v) {
	Vector<float, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = uintBitsToFloat(v[c]);
	}
	return ret;
}

inline unsigned int packSnorm2x16(const Vector<float, 2> &v) {
	unsigned short v0 = static_cast<unsigned short>(int(round(clamp(v[0], -1.0f, +1.0f) * 32767.0f)));
	unsigned short v1 = static_cast<unsigned short>(int(round(clamp(v[1], -1.0f, +1.0f) * 32767.0f)));
	return v0 | (v1 << 16);
}

inline Vector<float, 2> unpackSnorm2x16(unsigned int p) {
	short v0 = static_cast<short>(p);
	short v1 = static_cast<short>(p >> 16);
	return Vector<float, 2>(clamp(v0 / 32767.0f, -1.0f, +1.0f), clamp(v1 / 32767.0f, -1.0f, +1.0f));
}

inline unsigned int packUnorm2x16(const Vector<float, 2> &v) {
	unsigned short v0 = static_cast<unsigned short>(round(clamp(v[0], 0.0f, +1.0f) * 65535.0f));
	unsigned short v1 = static_cast<unsigned short>(round(clamp(v[1], 0.0f, +1.0f) * 65535.0f));
	return v0 | (v1 << 16);
}

inline Vector<float, 2> unpackUnorm2x16(unsigned int p) {
	unsigned short v0 = static_cast<unsigned short>(p);
	unsigned short v1 = static_cast<unsigned short>(p >> 16);
	return Vector<float, 2>(v0 / 65535.0f, v1 / 65535.0f);
}

}

#endif //GLAM_VECTOR_H
