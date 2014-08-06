/* Copyright (c) 2014, Gregor Riepl <onitake@gmail.com>
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

#ifndef GLAM_MATH_H
#define GLAM_MATH_H

#include <cmath>
#include <glam/config.h>

namespace glam {

// Convert degrees to radians (works for scalars and vectors, component-wise)
template <class T>
inline T radians(const T &degrees);
// Convert radians to degrees (works for scalars and vectors, component-wise)
template <class T>
inline T degrees(const T &radians);
// Sclar sine
template <class T>
inline T sin(const T &x);
// Scalar cosine
template <class T>
inline T cos(const T &x);
// Scalar tangent
template <class T>
inline T tan(const T &x);
// Scalar arc sine
template <class T>
inline T asin(const T &x);
// Scalar arc cosine
template <class T>
inline T acos(const T &x);
// Scalar arc tangent
template <class T>
inline T atan(const T &x);
// Scalar arc tangent of y/xs (works for scalars and vectors, component-wise)
template <class T>
inline T atan(const T &y, const T &x);
// Scalar sinus hyperbolicus
template <class T>
inline T sinh(const T &x);
// Scalar cosinus hyperbolicus
template <class T>
inline T cosh(const T &x);
// Scalar tangens hyperbolicus
template <class T>
inline T tanh(const T &x);
// Inverse of sinus hyperbolicus (works for scalars and vectors, component-wise)
template <class T>
inline T asinh(const T &x);
// Inverse of cosinus hyperbolicus (works for scalars and vectors, component-wise)
template <class T>
inline T acosh(const T &x);
// Inverse of tangens hyperbolicus (works for scalars and vectors, component-wise)
template <class T>
inline T atanh(const T &x);

// x raised to the power of y
template <class T>
inline T pow(const T &x, const T &y);
// e raised to the power of x
template <class T>
inline T exp(const T &x);
// Natual logarithm (base e)
template <class T>
inline T log(const T &x);
// 2 raised to the power of x
template <class T>
inline T exp2(const T &x);
// Base 2 logarithm
template <class T>
inline T log2(const T &x);
// Square root
template <class T>
inline T sqrt(const T &x);
// Inverse squareroot (scalar and component-wise)
template <class T>
inline T inversesqrt(const T &x);

// Absolute value |x|
template <class T>
inline T abs(const T &x);
// Sign (1 if positive, -1 if negative, 0 if 0)
template <class T>
inline T sign(const T &x);
// Round towards negative infinity
template <class T>
inline T floor(const T &x);
// Truncate (round towards 0)
template <class T>
inline T trunc(const T &x);
// Round to the nearest integer
template <class T>
inline T round(const T &x);
// Round to the nearest even integer
template <class T>
inline T roundEven(const T &x);
// Round towards positive infinity
template <class T>
inline T ceil(const T &x);
// Fractional part (x - trunc(x))
template <class T>
inline T fract(const T &x);
// Modulus (x - y * floor (x/y)), works for both scalar and vector arguments
template <class T>
inline T mod(const T &x, const T &y);
// Separate x into its fractional and integer parts. The fraction is returned, the integral assigned to i.
template <class T>
inline T modf(const T &x, T &i);
// Return the smaller value of x and y
template <class T>
inline T &min(T &x, T &y);
template <class T>
inline const T &min(const T &x, const T &y);
// Return the larger value of x and y
template <class T>
inline T &max(T &x, T &y);
template <class T>
inline const T &max(const T &x, const T &y);
// Saturation (min(max(x, minVal), maxVal)), works for both scalar and vector arguments
template <class T>
inline T clamp(const T &x, const T &minVal, const T &maxVal);
// Linear blend (x * (1 - a) + y * a)
template <class T>
inline T mix(const T &x, const T &y, const T &a);
// Return x if a is false, y if true
template <class T>
inline T mix(const T &x, const T &y, bool a);
// Step transition (0 if x < edge, 1 otherwise)
template <class T>
inline T step(const T &edge, const T &x);
// Smooth step transition (0 if x <= edge0, 1 if x >= edge1, interpolated with a Hermite term in between)
template <class T>
inline T smoothstep(const T &edge0, const T &edge1, const T &x);
// Return true if x is a NaN
template <class T>
inline bool isnan(const T &x);
// Return true if x is positive or negative infinity
template <class T>
inline bool isinf(const T &x);


// Implementation

template <class T>
T atan(const T &y, const T &x) {
	return atan(y / x);
}

template <class T>
inline T radians(const T &degrees) {
	return T(3.14159265358979323846 / 180.0) * degrees;
}

template <class T>
inline T degrees(const T &radians) {
	return T(180.0 / 3.14159265358979323846) * radians;
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

template <class T>
inline T log2(const T &x) {
	return std::log2(x);
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

template <class T>
inline T trunc(const T &x) {
	return std::trunc(x);
}

template <class T>
inline T round(const T &x) {
	return std::round(x);
}

template <class T>
inline T roundEven(const T &x) {
	return std::rint(x);
}

template <class T>
inline T fract(const T &x) {
	return x - trunc(x);
}

template <class T>
inline T mod(const T &x, const T &y) {
	return x - y * floor (x / y);
}

template <class T>
inline T modf(const T &x, T &i) {
	return std::modf(x, &i);
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

template <class T>
inline T mix(const T &x, const T &y, const T &a) {
	return x * (T(1) - a) + y * a;
}

template <class T>
inline T mix(const T &x, const T &y, bool a) {
	return a ? y : x;
}

template <class T>
inline T step(const T &edge, const T &x) {
	return x < edge ? T(0) : T(1);
}

template <class T>
inline T smoothstep(const T &edge0, const T &edge1, const T &x) {
	if (x <= edge0) {
		return T(0);
	} else if (x >= edge1) {
		return T(1);
	} else {
		T t = clamp((x - edge0) / (edge1 - edge0), T(0), T(1));
		return t * t * (T(3) - T(2) * t);
	}
}

}

#endif //GLAM_MATH_H
