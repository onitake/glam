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
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHType HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUType NOType LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENType SHALL THE COPYRIGHType HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUType NOType LIMITED TO, PROCUREMENType OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICType LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUType OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef GLAM_MATH_H
#define GLAM_MATH_H

#include <cmath>
#include <type_traits>
#include <glam/config.h>

namespace glam {

// Mathematical constants
// The generic implementation derives its values from double (or long double, if available) precision rational numbers.
// If your type needs more precision or different values, provide a template specialization.
#ifdef GLAM_MATH_CONST
template <typename Type>
struct CONSTANTS {
#ifdef GLAM_HAS_LONG_DOUBLE
	static constexpr Type PI = Type(3.141592653589793238462643383279502884L);
	static constexpr Type E = Type(2.718281828459045235360287471352662498L);
#else
	static constexpr Type PI = Type(3.14159265358979323846);
	static constexpr Type E = Type(2.7182818284590452354);
#endif
};
#endif


// Scalar/vector template discrimination helper
// Triggers SFINAE if Type is not a scalar type (or a compile-time error if no compatible overload is available)
template<typename Type>
using Scalar = typename std::enable_if<std::is_scalar<Type>::value, Type>::type;

// Convert degrees to radians (works for scalars and vectors, component-wise)
template <typename Type>
inline Type radians(const Type &degrees);
// Convert radians to degrees (works for scalars and vectors, component-wise)
template <typename Type>
inline Type degrees(const Type &radians);
// Sclar sine
template <typename Type>
inline Scalar<Type> sin(const Type &x);
// Scalar cosine
template <typename Type>
inline Scalar<Type> cos(const Type &x);
// Scalar tangent
template <typename Type>
inline Scalar<Type> tan(const Type &x);
// Scalar arc sine
template <typename Type>
inline Scalar<Type> asin(const Type &x);
// Scalar arc cosine
template <typename Type>
inline Scalar<Type> acos(const Type &x);
// Scalar arc tangent
template <typename Type>
inline Scalar<Type> atan(const Type &x);
// Scalar arc tangent of y/xs (works for scalars and vectors, component-wise)
template <typename Type>
inline Type atan(const Type &y, const Type &x);
// Scalar sinus hyperbolicus
template <typename Type>
inline Scalar<Type> sinh(const Type &x);
// Scalar cosinus hyperbolicus
template <typename Type>
inline Scalar<Type> cosh(const Type &x);
// Scalar tangens hyperbolicus
template <typename Type>
inline Scalar<Type> tanh(const Type &x);
// Inverse of sinus hyperbolicus (works for scalars and vectors, component-wise)
template <typename Type>
inline Type asinh(const Type &x);
// Inverse of cosinus hyperbolicus (works for scalars and vectors, component-wise)
template <typename Type>
inline Type acosh(const Type &x);
// Inverse of tangens hyperbolicus (works for scalars and vectors, component-wise)
template <typename Type>
inline Type atanh(const Type &x);

// x raised to the power of y
template <typename Type>
inline Scalar<Type> pow(const Type &x, const Type &y);
// e raised to the power of x
template <typename Type>
inline Scalar<Type> exp(const Type &x);
// Natual logarithm (base e)
template <typename Type>
inline Scalar<Type> log(const Type &x);
// 2 raised to the power of x
template <typename Type>
inline Scalar<Type> exp2(const Type &x);
// Base 2 logarithm
template <typename Type>
inline Scalar<Type> log2(const Type &x);
// Square root
template <typename Type>
inline Scalar<Type> sqrt(const Type &x);
// Inverse squareroot (scalar and component-wise)
template <typename Type>
inline Type inversesqrt(const Type &x);

// Absolute value |x|
template <typename Type>
inline Scalar<Type> abs(const Type &x);
// Sign (1 if positive, -1 if negative, 0 if 0)
template <typename Type>
inline Scalar<Type> sign(const Type &x);
// Round towards negative infinity
template <typename Type>
inline Scalar<Type> floor(const Type &x);
// Truncate (round towards 0)
template <typename Type>
inline Scalar<Type> trunc(const Type &x);
// Round to the nearest integer
template <typename Type>
inline Scalar<Type> round(const Type &x);
// Round to the nearest even integer
template <typename Type>
inline Scalar<Type> roundEven(const Type &x);
// Round towards positive infinity
template <typename Type>
inline Scalar<Type> ceil(const Type &x);
// Fractional part (x - trunc(x)), works for both scalars and vectors.
template <typename Type>
inline Type fract(const Type &x);
// Modulus (x - y * floor (x/y)), works for both scalar and vector arguments
template <typename Type>
inline Type mod(const Type &x, const Type &y);
// Separate x into its fractional and integer parts. The fraction is returned, the integral assigned to i.
template <typename Type>
inline Scalar<Type> modf(const Type &x, Type &i);
// Return the smaller value of x and y
template <typename Type>
inline Scalar<Type> min(const Type &x, const Type &y);
// Return the larger value of x and y
template <typename Type>
inline Scalar<Type> max(const Type &x, const Type &y);
// Saturation (min(max(x, minVal), maxVal)), works for both scalar and vector arguments
template <typename Type>
inline Type clamp(const Type &x, const Type &minVal, const Type &maxVal);
// Linear blend (x * (1 - a) + y * a)
template <typename Type>
inline Type mix(const Type &x, const Type &y, const Type &a);
// Return x if a is false, y if true
template <typename Type>
inline Scalar<Type> mix(const Type &x, const Type &y, bool a);
// Step transition (0 if x < edge, 1 otherwise)
template <typename Type>
inline Scalar<Type> step(const Type &edge, const Type &x);
// Smooth step transition (0 if x <= edge0, 1 if x >= edge1, interpolated with a Hermite term in between)
template <typename Type>
inline Scalar<Type> smoothstep(const Type &edge0, const Type &edge1, const Type &x);
// Return true if x is a NaN
template <typename Type>
inline bool isnan(const Scalar<Type> &x);
// Return true if x is positive or negative infinity
template <typename Type>
inline bool isinf(const Scalar<Type> &x);

// Convert a floating point number into an integer of equal size and exactly the same bit pattern
inline int floatBitsToInt(float v);
inline unsigned int floatBitsToUInt(float v);
// Convert an integer into a floating point number of equal size and exactly the same bit pattern
inline float intBitsToFloat(int v);
inline float uintBitsToFloat(unsigned int v);


// GLSL types

typedef unsigned int uint;


// Implementation

template <typename Type>
inline Type atan(const Type &y, const Type &x) {
	return atan(y / x);
}

template <typename Type>
inline Type radians(const Type &degrees) {
#ifdef GLAM_MATH_CONST
	return degrees * CONSTANTS<Type>::PI / Type(180);
#else
#ifdef GLAM_HAS_LONG_DOUBLE
	return degrees * Type(3.141592653589793238462643383279502884L) / Type(180.0);
#else
	return degrees * Type(3.14159265358979323846) / Type(180.0);
#endif
#endif
}

template <typename Type>
inline Type degrees(const Type &radians) {
#ifdef GLAM_MATH_CONST
	return radians * Type(180) / CONSTANTS<Type>::PI;
#else
#ifdef GLAM_HAS_LONG_DOUBLE
	return radians * Type(180) / Type(3.141592653589793238462643383279502884L);
#else
	return radians * Type(180) / Type(3.14159265358979323846);
#endif
#endif
}

template <typename Type>
inline Type asinh(const Type &x) {
	return log(x + sqrt(x * x + Type(1)));
}

template <typename Type>
inline Type acosh(const Type &x) {
	return log(x + sqrt(x * x - Type(1)));
}

template <typename Type>
inline Type atanh(const Type &x) {
	return Type(0.5) * log((Type(1) + x) / (Type(1) - x));
}

template <typename Type>
inline Scalar<Type> exp2(const Type &x) {
	return std::exp2(x);
}

template <typename Type>
inline Scalar<Type> log2(const Type &x) {
	return std::log2(x);
}

template <typename Type>
inline Type inversesqrt(const Type &x) {
	return Type(1) / sqrt(x);
}

template <typename Type>
inline Scalar<Type> sign(const Type &x) {
#if 0
	if (x == Type(0)) return Type(0);
	return std::copysign(Type(1), x);
#else
	if (x < Type(0)) return Type(-1);
	if (x > Type(0)) return Type(1);
	return Type(0);
#endif
}

template <typename Type>
inline Scalar<Type> trunc(const Type &x) {
	return std::trunc(x);
}

template <typename Type>
inline Scalar<Type> round(const Type &x) {
	return std::round(x);
}

template <typename Type>
inline Scalar<Type> roundEven(const Type &x) {
	return std::rint(x);
}

template <typename Type>
inline Type fract(const Type &x) {
	return x - trunc(x);
}

template <typename Type>
inline Type mod(const Type &x, const Type &y) {
	return x - y * floor (x / y);
}

template <typename Type>
inline Scalar<Type> modf(const Type &x, Type &i) {
	return std::modf(x, &i);
}

template <typename Type>
inline Scalar<Type> sqrt(const Type &x) {
	return std::sqrt(x);
}

template <typename Type>
inline Scalar<Type> sin(const Type &x) {
	return std::sin(x);
}

template <typename Type>
inline Scalar<Type> cos(const Type &x) {
	return std::cos(x);
}

template <typename Type>
inline Scalar<Type> tan(const Type &x) {
	return std::tan(x);
}

template <typename Type>
inline Scalar<Type> asin(const Type &x) {
	return std::asin(x);
}

template <typename Type>
inline Scalar<Type> acos(const Type &x) {
	return std::acos(x);
}

template <typename Type>
inline Scalar<Type> atan(const Type &x) {
	return std::atan(x);
}

template <typename Type>
inline Scalar<Type> sinh(const Type &x) {
	return std::sinh(x);
}

template <typename Type>
inline Scalar<Type> cosh(const Type &x) {
	return std::cosh(x);
}

template <typename Type>
inline Scalar<Type> tanh(const Type &x) {
	return std::tanh(x);
}

template <typename Type>
inline Scalar<Type> pow(const Type &x, const Type &y) {
	return std::pow(x, y);
}

template <typename Type>
inline Scalar<Type> exp(const Type &x) {
	return std::exp(x);
}

template <typename Type>
inline Scalar<Type> log(const Type &x) {
	return std::log(x);
}

template <typename Type>
inline Scalar<Type> abs(const Type &x) {
	return std::abs(x);
}

template <typename Type>
inline Scalar<Type> floor(const Type &x) {
	return std::floor(x);
}

template <typename Type>
inline Scalar<Type> ceil(const Type &x) {
	return std::ceil(x);
}

template <typename Type>
inline Scalar<Type> min(const Type &x, const Type &y) {
	if (x < y) return x;
	return y;
}

template <typename Type>
inline Scalar<Type> max(const Type &x, const Type &y) {
	if (x > y) return x;
	return y;
}

template <typename Type>
inline bool isnan(const Type &x) {
	return std::isnan(x);
}

template <typename Type>
inline bool isinf(const Type &x) {
	return std::isinf(x);
}

template <typename Type>
inline Type clamp(const Type &x, const Type &minVal, const Type &maxVal) {
	return min(max(x, minVal), maxVal);
}

template <typename Type>
inline Type mix(const Type &x, const Type &y, const Type &a) {
	return x * (Type(1) - a) + y * a;
}

template <typename Type>
inline Scalar<Type> mix(const Type &x, const Type &y, bool a) {
	return a ? y : x;
}

template <typename Type>
inline Scalar<Type> step(const Type &edge, const Type &x) {
	return x < edge ? Type(0) : Type(1);
}

template <typename Type>
inline Scalar<Type> smoothstep(const Type &edge0, const Type &edge1, const Type &x) {
	if (x <= edge0) {
		return Type(0);
	} else if (x >= edge1) {
		return Type(1);
	} else {
		Type t = clamp((x - edge0) / (edge1 - edge0), Type(0), Type(1));
		return t * t * (Type(3) - Type(2) * t);
	}
}

inline int floatBitsToInt(float v) {
	union {
		float f;
		int i;
	} c;
	c.f = v;
	return c.i;
}

inline unsigned int floatBitsToUInt(float v) {
	union {
		float f;
		unsigned int i;
	} c;
	c.f = v;
	return c.i;
}

inline float intBitsToFloat(int v) {
	union {
		float f;
		int i;
	} c;
	c.i = v;
	return c.f;
}

inline float uintBitsToFloat(unsigned int v) {
	union {
		float f;
		unsigned int i;
	} c;
	c.i = v;
	return c.f;
}

}

#endif //GLAM_MATH_H
