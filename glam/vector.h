/* Copyright (c) 2013-2014, Gregor Riepl <onitake@gmail.com>
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

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <sstream>
#include <glam/config.h>
#include <glam/exception.h>

namespace glam {

// Declaration

template <class T, unsigned int C>
class Vector {
private:
	T _v[C];

	// Overloadable template functions for use by the generator macro
	T permutation(unsigned int x) const;
	Vector<T, 2> permutation(unsigned int x, unsigned int y) const;
	Vector<T, 3> permutation(unsigned int x, unsigned int y, unsigned int z) const;
	Vector<T, 4> permutation(unsigned int x, unsigned int y, unsigned int z, unsigned int w) const;
	
	// Supporting initializers to avoid trouble from delegate constructors
	template <unsigned int X>
	void initialize();
	template <unsigned int X, typename ForwardIterator>
	void initialize(std::pair<ForwardIterator, ForwardIterator> iterator);
	template <unsigned int X>
	void initialize(const T *v);
	template <unsigned int X, class U>
	void initialize(const Vector<U, C> &v);
	template <unsigned int X>
	void initialize(const T &s);
	template <unsigned int X, typename... Args, unsigned int P>
	void initialize(const Vector<T, P> &v, Args... args);
	template <unsigned int X, typename... Args>
	void initialize(const T &v, Args... args);
	
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
	Vector(Args... args);

	// Assignment operator
	Vector<T, C> &operator =(const Vector<T, C> &v);
	
	// Member access operator
	T &operator [](unsigned int index);
	// Const member access operator
	const T &operator [](unsigned int index) const;
	// Direct access to internal pointer
	const T *internal() const;

	// Permutation functions
	// These are in the form of function calls, not member access like in GLSL
	// Eg.: v.xxyy() instead of v.xxyy
	// All the 340 1-, 2-, 3- and 4-member permutations are generated on the fly.
	#define GLAM_PERM_MAKE1_1(type, length, name, member0) T name() const { static_assert(member0 < length, "Invalid member index for this vector type"); return permutation(member0); }
	#define GLAM_PERM_MAKE2_2(type, length, name, member0, member1) Vector<type, 2> name() const { static_assert(member0 < length && member1 < length, "Invalid member index for this vector type"); return permutation(member0, member1); }
	#define GLAM_PERM_MAKE3_3(type, length, name, member0, member1, member2) Vector<type, 3> name() const { static_assert(member0 < length && member1 < length && member2 < length, "Invalid member index for this vector type"); return permutation(member0, member1, member2); }
	#define GLAM_PERM_MAKE4_4(type, length, name, member0, member1, member2, member3) Vector<type, 4> name() const { static_assert(member0 < length && member1 < length && member2 < length && member3 < length, "Invalid member index for this vector type"); return permutation(member0, member1, member2, member3); }
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
Vector<T, C> &operator ++(Vector<T, C> &self);
// Component-wise post-increment operator
template <class T, unsigned int C>
Vector<T, C> operator ++(Vector<T, C> &self, int);
// Component-wise pre-decrement operator
template <class T, unsigned int C>
Vector<T, C> &operator --(Vector<T, C> &self);
// Component-wise post-decrement operator
template <class T, unsigned int C>
Vector<T, C> operator --(Vector<T, C> &self, int);

// Component-wise self-addition operator
template <class T, unsigned int C>
Vector<T, C> &operator +=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise self-subtraction operator
template <class T, unsigned int C>
Vector<T, C> &operator -=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise self-multiplication operator
template <class T, unsigned int C>
Vector<T, C> &operator *=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise self-division operator
template <class T, unsigned int C>
Vector<T, C> &operator /=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise self-modulus operator
template <class T, unsigned int C>
Vector<T, C> &operator %=(Vector<T, C> &self, const Vector<T, C> &v);

// Component-wise bitwise and operator
template <class T, unsigned int C>
Vector<T, C> &operator &=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise bitwise or operator
template <class T, unsigned int C>
Vector<T, C> &operator |=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise bitwise xor operator
template <class T, unsigned int C>
Vector<T, C> &operator ^=(Vector<T, C> &self, const Vector<T, C> &v);
// Component-wise bit shift left operator
template <class T, unsigned int C>
Vector<T, C> &operator <<=(Vector<T, C> &self, unsigned int s);
// Component-wise bit shift right operator
template <class T, unsigned int C>
Vector<T, C> &operator >>=(Vector<T, C> &self, unsigned int s);

// Component-wise unary plus operator (no-op if T follows standard algebra)
template <class T, unsigned int C>
Vector<T, C> operator +(const Vector<T, C> &u);
// Component-wise negation operator (negation)
template <class T, unsigned int C>
Vector<T, C> operator -(const Vector<T, C> &u);
// Component-wise ones complement operator
template <class T, unsigned int C>
Vector<T, C> operator ~(const Vector<T, C> &u);

// Component-wise addition operator
template <class T, unsigned int C>
Vector<T, C> operator +(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise subtraction operator
template <class T, unsigned int C>
Vector<T, C> operator -(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise multiplication operator
template <class T, unsigned int C>
Vector<T, C> operator *(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise division operator
template <class T, unsigned int C>
Vector<T, C> operator /(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise modulus operator
template <class T, unsigned int C>
Vector<T, C> operator %(const Vector<T, C> &u, const Vector<T, C> &v);

// Component-wise bitwise and operator
template <class T, unsigned int C>
Vector<T, C> operator &(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise bitwise or operator
template <class T, unsigned int C>
Vector<T, C> operator |(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise bitwise xor operator
template <class T, unsigned int C>
Vector<T, C> operator ^(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise bit shift left operator
template <class T, unsigned int C>
Vector<T, C> operator <<(const Vector<T, C> &u, unsigned int s);
// Component-wise bit shift right operator
template <class T, unsigned int C>
Vector<T, C> operator >>(const Vector<T, C> &u, unsigned int s);

// Component-wise equality comparison operator (exact compare)
template <class T, unsigned int C>
Vector<bool, C> equal(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise inequality comparison operator (exact compare)
template <class T, unsigned int C>
Vector<bool, C> notEqual(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise less-than comparison operator (exact compare)
template <class T, unsigned int C>
Vector<bool, C> lessThan(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise greater-than comparison operator (exact compare)
template <class T, unsigned int C>
Vector<bool, C> greaterThan(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise less-or-equal comparison operator (exact compare)
template <class T, unsigned int C>
Vector<bool, C> lessThanEqual(const Vector<T, C> &u, const Vector<T, C> &v);
// Component-wise greater-or-equal comparison operator (exact compare)
template <class T, unsigned int C>
Vector<bool, C> greaterThanEqual(const Vector<T, C> &u, const Vector<T, C> &v);

// Vector comparison operator (exact compare), correspondes to all(equal(u, v))
template <class T, unsigned int C>
bool operator ==(const Vector<T, C> &u, const Vector<T, C> &v);
// Vector comparison operator (exact compare), correspondes to any(notEqual(u, v))
template <class T, unsigned int C>
bool operator !=(const Vector<T, C> &u, const Vector<T, C> &v);

// Convert degrees to radians (works for scalars and vectors, component-wise)
template <class T>
T radians(const T &degrees);
// Convert radians to degrees (works for scalars and vectors, component-wise)
template <class T>
T degrees(const T &radians);
// Sclar sine
template <class T>
T sin(const T &x);
// Component-wise sine
template <class T, unsigned int C>
Vector<T, C> sin(const Vector<T, C> &v);
// Scalar cosine
template <class T>
T cos(const T &x);
// Component-wise cosine
template <class T, unsigned int C>
Vector<T, C> cos(const Vector<T, C> &v);
// Scalar tangent
template <class T>
T tan(const T &x);
// Component-wise tangent
template <class T, unsigned int C>
Vector<T, C> tan(const Vector<T, C> &v);
// Scalar arc sine
template <class T>
T asin(const T &x);
// Component-wise arc sine
template <class T, unsigned int C>
Vector<T, C> asin(const Vector<T, C> &v);
// Scalar arc cosine
template <class T>
T acos(const T &x);
// Component-wise arc cosine
template <class T, unsigned int C>
Vector<T, C> acos(const Vector<T, C> &v);
// Scalar arc tangent
template <class T>
T atan(const T &x);
// Component-wise arc tangent
template <class T, unsigned int C>
Vector<T, C> atan(const Vector<T, C> &v);
// Scalar arc tangent of y/xs (works for scalars and vectors, component-wise)
template <class T>
T atan(const T &y, const T &x);
// Scalar sinus hyperbolicus
template <class T>
T sinh(const T &x);
// Component-wise sinus hyperbolicus
template <class T, unsigned int C>
Vector<T, C> sinh(const Vector<T, C> &v);
// Scalar cosinus hyperbolicus
template <class T>
T cosh(const T &x);
// Component-wise cosinus hyperbolicus
template <class T, unsigned int C>
Vector<T, C> cosh(const Vector<T, C> &v);
// Scalar tangens hyperbolicus
template <class T>
T tanh(const T &x);
// Component-wise tangens hyperbolicus
template <class T, unsigned int C>
Vector<T, C> tanh(const Vector<T, C> &v);
// Inverse of sinus hyperbolicus (works for scalars and vectors, component-wise)
template <class T>
T asinh(const T &x);
// Inverse of cosinus hyperbolicus (works for scalars and vectors, component-wise)
template <class T>
T acosh(const T &x);
// Inverse of tangens hyperbolicus (works for scalars and vectors, component-wise)
template <class T>
T atanh(const T &x);

// x raised to the power of y
template <class T>
T pow(const T &x, const T &y);
// Component-wise power
template <class T, unsigned int C>
Vector<T, C> pow(const Vector<T, C> &x, const Vector<T, C> &y);
// e raised to the power of x
template <class T>
T exp(const T &x);
// Component-wise exponential (e^v)
template <class T, unsigned int C>
Vector<T, C> exp(const Vector<T, C> &v);
// Natual logarithm (base e)
template <class T>
T log(const T &x);
// Component-wise natural logarithm
template <class T, unsigned int C>
Vector<T, C> log(const Vector<T, C> &v);
// 2 raised to the power of x
template <class T>
T exp2(const T &x);
// Component-wise power of 2
template <class T, unsigned int C>
Vector<T, C> exp2(const Vector<T, C> &v);
// Base 2 logarithm
template <class T>
T log2(const T &x);
// Component-wise base 2 logarithm
template <class T, unsigned int C>
Vector<T, C> log2(const Vector<T, C> &v);
// Square root
template <class T>
T sqrt(const T &x);
// Component-wise square root
template <class T, unsigned int C>
Vector<T, C> sqrt(const Vector<T, C> &v);
// Inverse squareroot
template <class T>
T inversesqrt(const T &x);

// Absolute value |x|
template <class T>
T abs(const T &x);
// Component-wise absolute value
template <class T, unsigned int C>
Vector<T, C> abs(const Vector<T, C> &v);
// Sign (1 if positive, -1 if negative, 0 if 0)
template <class T>
T sign(const T &x);
// Component-wise sign
template <class T, unsigned int C>
Vector<T, C> sign(const Vector<T, C> &v);
// Round towards negative infinity
template <class T>
T floor(const T &x);
// Component-wise round towards negative infinity
template <class T, unsigned int C>
Vector<T, C> floor(const Vector<T, C> &v);
// Truncate (round towards 0)
template <class T>
T trunc(const T &x);
// Component-wise truncate
template <class T, unsigned int C>
Vector<T, C> trunc(const Vector<T, C> &v);
// Round to the nearest integer
template <class T>
T round(const T &x);
// Component-wise round to the nearest integer
template <class T, unsigned int C>
Vector<T, C> round(const Vector<T, C> &v);
// Round to the nearest even integer
template <class T>
T roundEven(const T &x);
// Component-wise round to the nearest even integer
template <class T, unsigned int C>
Vector<T, C> roundEven(const Vector<T, C> &v);
// Round towards positive infinity
template <class T>
T ceil(const T &x);
// Component-wise round towards positive infinity
template <class T, unsigned int C>
Vector<T, C> ceil(const Vector<T, C> &v);
// Fractional part (x - trunc(x))
template <class T>
T fract(const T &x);
// Component-wise fractional part
template <class T, unsigned int C>
Vector<T, C> fract(const Vector<T, C> &v);
// Modulus (x - y * floor (x/y)), works for both scalar and vector arguments
template <class T>
T mod(const T &x, const T &y);
// Separate x into its fractional and integer parts. The fraction is returned, the integral assigned to i.
template <class T>
T modf(const T &x, T &i);
// Component-wise fractional part separation
template <class T, unsigned int C>
Vector<T, C> modf(const Vector<T, C> &v, Vector<T, C> &i);
// Return the smaller value of x and y
template <class T>
T &min(T &x, T &y);
template <class T>
const T &min(const T &x, const T &y);
// Component-wise minimum
template <class T, unsigned int C>
Vector<T, C> min(const Vector<T, C> &x, const Vector<T, C> &y);
// Return the larger value of x and y
template <class T>
T &max(T &x, T &y);
template <class T>
const T &max(const T &x, const T &y);
// Component-wise maximum
template <class T, unsigned int C>
Vector<T, C> max(const Vector<T, C> &x, const Vector<T, C> &y);
// Saturation (min(max(x, minVal), maxVal)), works for both scalar and vector arguments
template <class T>
T clamp(const T &x, const T &minVal, const T &maxVal);
// Saturation with single min/max values (applied to all vector components)
template <class T, unsigned int C>
T clamp(const Vector<T, C> &x, const T &minVal, const T &maxVal);
// Linear blend (x * (1 - a) + y * a)
template <class T>
T mix(const T &x, const T &y, const T &a);
// Component-wise linear blend (x * (1 - a) + y * a)
template <class T, unsigned int C>
Vector<T, C> mix(const Vector<T, C> &x, const Vector<T, C> &y, const T &a);
// Return x if a is false, y if true
template <class T>
T mix(const T &x, const T &y, bool a);
// Return x for each component of a that is false, y if true
template <class T, unsigned int C>
Vector<T, C> mix(const Vector<T, C> &x, const Vector<T, C> &y, const Vector<bool, C> &a);
// Step transition (0 if x < edge, 1 otherwise)
template <class T>
T step(const T &edge, const T &x);
// Component-wise step transition
template <class T, unsigned int C>
Vector<T, C> step(const Vector<T, C> &edge, const Vector<T, C> &x);
// Component-wise step transition with common edge argument
template <class T, unsigned int C>
Vector<T, C> step(const T &edge, const Vector<T, C> &x);
// Smooth step transition (0 if x <= edge0, 1 if x >= edge1, interpolated with a Hermite term in between)
template <class T>
T smoothstep(const T &edge0, const T &edge1, const T &x);
// Component-wise smooth step transition
template <class T, unsigned int C>
Vector<T, C> smoothstep(const Vector<T, C> &edge0, const Vector<T, C> &edge1, const Vector<T, C> &x);
// Component-wise smooth step transition with common edge arguments
template <class T, unsigned int C>
Vector<T, C> smoothstep(const T &edge0, const T &edge1, const Vector<T, C> &x);
// Return true if x is a NaN
template <class T>
bool isnan(const T &x);
// Component-wise test for NaN
template <class T, unsigned int C>
Vector<bool, C> isnan(const Vector<T, C> &v);
// Return true if x is positive or negative infinity
template <class T>
bool isinf(const T &x);
// Component-wise test for infinity
template <class T, unsigned int C>
Vector<bool, C> isinf(const Vector<T, C> &v);

// Vector length (sqrt(dot(a, a)))
template <class T, unsigned int C>
T length(const Vector<T, C> &x);
// Distance between two vectors (length(b - a))
template <class T, unsigned int C>
T distance(const Vector<T, C> &p0, const Vector<T, C> &p1);
// Dot product of two vectors (a[0] * b[0] + a[1] * b[1] + ...)
template <class T, unsigned int C>
T dot(const Vector<T, C> &x, const Vector<T, C> &y);
// Cross product of two 3-component vectors
template <class T>
Vector<T, 3> cross(const Vector<T, 3> &x, const Vector<T, 3> &y);
// Unit length vector (a / length(a))
// always_inline encouraged to keep gcc from generating crap code
template <class T, unsigned int C>
Vector<T, C> normalize(const Vector<T, C> &x) __attribute__((always_inline));
// Check if I and Nref are facing in opposite directions and return N if yes, -N otherwise (dot(Nref, I) < 0 ? N : -N)
template <class T>
T faceforward(const T &N, const T &I, const T &Nref);
// Reflect I in respect to surface normal N
template <class T, unsigned int C>
Vector<T, C> reflect(const Vector<T, C> &I, const Vector<T, C> &N);
// Refract I on the surface defined by normal N using index of refraction ration eta
template <class T, unsigned int C>
Vector<T, C> refract(const Vector<T, C> &I, const Vector<T, C> &N, const T &eta);

// True if any component of x is true (or evaluates to true)
template <unsigned int C>
bool any(const Vector<bool, C> &v);
// True if all components of x are true (or evaluates to true)
template <unsigned int C>
bool all(const Vector<bool, C> &v);
// Component-wise negation of a boolean vector, this is an overloaded ! operator because
// not is a reserved word in C++ with the same semantics as !
template <unsigned int C>
Vector<bool, C> operator !(const Vector<bool, C> &v);

// Formatted output operator (vectors will be presented in the form "(0.1,1.50)")
template <class T, unsigned int C>
std::ostream &operator <<(std::ostream &stream, const Vector<T, C> &v);
// Read vector components from stream, in text form, separated by whitespace
template <class T, unsigned int C>
std::istream &operator >>(std::istream &stream, const Vector<T, C> &v);


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
	if (index < C) {
		return _v[index];
	}
	throw DimensionOutOfRangeException("Invalid vector element index", C - 1, index);
}

template <class T, unsigned int C>
inline const T &Vector<T, C>::operator [](unsigned int index) const {
	if (index < C) {
		return _v[index];
	}
	throw DimensionOutOfRangeException("Invalid vector element index", C - 1, index);
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

template <class T>
T atan(const T &y, const T &x) {
	return atan(y / x);
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

template <class T>
inline T radians(const T &degrees) {
	return T(M_PI / 180.0) * degrees;
}

template <class T>
inline T degrees(const T &radians) {
	return T(180.0 / M_PI) * radians;
}

template <class T>
inline T asinh(const T &x) {
	return log(x + sqrt(x * x + T(1)));
}

template <class T>
inline T acosh(const T &x) {
	return log(x + sqrt(x * x - T(1)));
}

template <class T>
inline T atanh(const T &x) {
	return T(0.5) * log((T(1) + x) / (T(1) - x));
}

template <class T>
inline T exp2(const T &x) {
	return std::exp2(x);
}

template <class T, unsigned int C>
Vector<T, C> exp2(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = exp2(v[c]);
	}
	return ret;
}

template <class T>
inline T log2(const T &x) {
	return std::log2(x);
}

template <class T, unsigned int C>
Vector<T, C> log2(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = log2(v[c]);
	}
	return ret;
}

template <class T>
inline T inversesqrt(const T &x) {
	return T(1) / sqrt(x);
}

template <class T>
inline T sign(const T &x) {
#if 0
	if (x == T(0)) return T(0);
	return std::copysign(T(1), x);
#else
	if (x < T(0)) return T(-1);
	if (x > T(0)) return T(1);
	return T(0);
#endif
}

template <class T, unsigned int C>
inline Vector<T, C> sign(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = sign(v[c]);
	}
	return ret;
}

template <class T>
inline T trunc(const T &x) {
	return std::trunc(x);
}

template <class T, unsigned int C>
Vector<T, C> trunc(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = trunc(v[c]);
	}
	return ret;
}

template <class T>
inline T round(const T &x) {
	return std::round(x);
}

template <class T, unsigned int C>
Vector<T, C> round(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = round(v[c]);
	}
	return ret;
}

template <class T>
inline T roundEven(const T &x) {
	return std::rint(x);
}

template <class T, unsigned int C>
Vector<T, C> roundEven(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = roundEven(v[c]);
	}
	return ret;
}

template <class T>
inline T fract(const T &x) {
	return x - trunc(x);
}

template <class T, unsigned int C>
Vector<T, C> fract(const Vector<T, C> &v) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = fract(v[c]);
	}
	return ret;
}

template <class T>
inline T mod(const T &x, const T &y) {
	return x - y * floor (x / y);
}

template <class T>
inline T modf(const T &x, T &i) {
	return std::modf(x, &i);
}

template <class T, unsigned int C>
Vector<T, C> modf(const Vector<T, C> &v, Vector<T, C> &i) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = modf(v[c], i[c]);
	}
	return ret;
}

template <class T>
inline T sqrt(const T &x) {
	return std::sqrt(x);
}

template <class T>
inline T sin(const T &x) {
	return std::sin(x);
}

template <class T>
inline T cos(const T &x) {
	return std::cos(x);
}

template <class T>
inline T tan(const T &x) {
	return std::tan(x);
}

template <class T>
inline T asin(const T &x) {
	return std::asin(x);
}

template <class T>
inline T acos(const T &x) {
	return std::acos(x);
}

template <class T>
inline T atan(const T &x) {
	return std::atan(x);
}

template <class T>
inline T sinh(const T &x) {
	return std::sinh(x);
}

template <class T>
inline T cosh(const T &x) {
	return std::cosh(x);
}

template <class T>
inline T tanh(const T &x) {
	return std::tanh(x);
}

template <class T>
inline T pow(const T &x, const T &y) {
	return std::pow(x, y);
}

template <class T>
inline T exp(const T &x) {
	return std::exp(x);
}

template <class T>
inline T log(const T &x) {
	return std::log(x);
}

template <class T>
inline T abs(const T &x) {
	return std::abs(x);
}

template <class T>
inline T floor(const T &x) {
	return std::floor(x);
}

template <class T>
inline T ceil(const T &x) {
	return std::ceil(x);
}

template <class T>
inline T &min(T &x, T &y) {
	if (x < y) return x;
	return y;
}

template <class T>
inline const T &min(const T &x, const T &y) {
	if (x < y) return x;
	return y;
}

template <class T>
inline T &max(T &x, T &y) {
	if (x > y) return x;
	return y;
}

template <class T>
inline const T &max(const T &x, const T &y) {
	if (x > y) return x;
	return y;
}

template <class T>
inline bool isnan(const T &x) {
	return std::isnan(x);
}

template <class T>
inline bool isinf(const T &x) {
	return std::isinf(x);
}

template <class T>
inline T clamp(const T &x, const T &minVal, const T &maxVal) {
	return min(max(x, minVal), maxVal);
}

template <class T, unsigned int C>
T clamp(const Vector<T, C> &x, const T &minVal, const T &maxVal) {
	return min(max(x, Vector<T, C>(minVal)), Vector<T, C>(maxVal));
}

template <class T>
inline T mix(const T &x, const T &y, const T &a) {
	return x * (T(1) - a) + y * a;
}

template <class T, unsigned int C>
Vector<T, C> mix(const Vector<T, C> &x, const Vector<T, C> &y, const T &a) {
	return x * Vector<T, C>(T(1) - a) + y * Vector<T, C>(a);
}

template <class T>
T mix(const T &x, const T &y, bool a) {
	return a ? y : x;
}

template <class T, unsigned int C>
Vector<T, C> mix(const Vector<T, C> &x, const Vector<T, C> &y, const Vector<bool, C> &a) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = mix(x[c], y[c], a[c]);
	}
	return ret;
}

template <class T>
T step(const T &edge, const T &x) {
	return x < edge ? T(0) : T(1);
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
Vector<T, C> step(const T &edge, const Vector<T, C> &x) {
	Vector<T, C> ret;
	for (unsigned int c = 0; c < C; c++) {
		ret[c] = step(edge, x[c]);
	}
	return ret;
}

template <class T>
T smoothstep(const T &edge0, const T &edge1, const T &x) {
	if (x <= edge0) {
		return T(0);
	} else if (x >= edge1) {
		return T(1);
	} else {
		T t = clamp((x - edge0) / (edge1 - edge0), T(0), T(1));
		return t * t * (3 - 2 * t);
	}
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
Vector<T, C> smoothstep(const T &edge0, const T &edge1, const Vector<T, C> &x) {
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
bool any(const Vector<bool, C> &v) {
	bool ret = false;
	for (unsigned int c = 0; c < C; c++) {
		ret = ret || v[c];
	}
	return ret;
}

template <unsigned int C>
bool all(const Vector<bool, C> &v) {
	bool ret = true;
	for (unsigned int c = 0; c < C; c++) {
		ret = ret && v[c];
	}
	return ret;
}

template <unsigned int C>
Vector<bool, C> operator !(const Vector<bool, C> &v) {
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
const T *Vector<T, C>::internal() const {
	return _v;
}

}

#endif //GLAM_VECTOR_H
